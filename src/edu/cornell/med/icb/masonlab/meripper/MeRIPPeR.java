package edu.cornell.med.icb.masonlab.meripper;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.IntervalTree;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import edu.cornell.med.icb.masonlab.jenotator.activity.GetChromosomeSizes;
import edu.cornell.med.icb.masonlab.meripper.util.FishersTestJunctionThread;
import edu.cornell.med.icb.masonlab.meripper.util.FishersTestThread;
import edu.cornell.med.icb.masonlab.meripper.util.Interval;
import edu.cornell.med.icb.masonlab.meripper.util.ReadCounter;
import edu.cornell.med.icb.masonlab.meripper.util.ReadCounterJunctions;
import edu.cornell.med.icb.masonlab.meripper.util.WindowFilterJunctionsThread;
import edu.cornell.med.icb.masonlab.meripper.util.WindowFilterThread;

public class MeRIPPeR {
	private static int WINDOW_SIZE = 25;
	private static int STEP_SIZE = WINDOW_SIZE;
	private static Map<String, Integer> GENOME_SIZES;
	private static int NUM_THREADS = 1;
	private static double ALPHA = 0.05;
	private static String P_ADJUST = "BenjaminiHochberg";
	private static int MIN_WINDOW_SIZE = 100;
	private static boolean JUNCTIONS = false;

	public static void main(String[] args) {
		CommandLineParser parser = new GnuParser();
		Options options = new Options();
		buildOptions(options);
		
		try {
			CommandLine cmd = parser.parse(options, args, true);
			processCmd(cmd);
			ExecutorService threadPool = Executors.newFixedThreadPool(NUM_THREADS);
			PrintStream windowPrintStream = new PrintStream(new File(cmd.getOptionValue("output")));
			
			// window maps
			Map<String, Map<Integer, Double>> significantWindows = new HashMap<String, Map<Integer, Double>>();
			Map<String, IntervalTree<Double>> junctions = new HashMap<String, IntervalTree<Double>>();
			Map<String, Map<Integer, Set<Integer>>> junctions_map = new HashMap<String, Map<Integer, Set<Integer>>>();
			
			// counts & fishers
			windowCountsFishers(cmd, threadPool, significantWindows, junctions, junctions_map);
			
			// p-value adjustment
			pValueAdjustment(significantWindows);
			
			// print out the final windows
			printFilteredWindows(threadPool, significantWindows, junctions, junctions_map, windowPrintStream);
			
			// shutdown the threadpool
			threadPool.shutdown();
			while(!threadPool.isShutdown()) {
				Thread.sleep(1000);
			}
			
			windowPrintStream.close();
		} catch(ParseException e) {
			e.printStackTrace();
			System.err.println("Improper arguments");
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp("MeRIPPER", options, true);
			System.exit(1); // so that make stops running
		} catch (IOException e) {
			e.printStackTrace();
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp("MeRIPPER", options, true);
			System.exit(1); // so that make stops running
		} catch (Throwable t) {
			t.printStackTrace();
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp("MeRIPPER", options, true);
			System.exit(1); // so that make stops running
		}
	}
	
	public static void processCmd(CommandLine cmd) throws IOException, ParseException {
		String genomeSizesFilename = cmd.getOptionValue("genome-sizes");
		
		if(cmd.hasOption("alpha")) {
			ALPHA = Double.parseDouble(cmd.getOptionValue("alpha"));
		}
		
		if(cmd.hasOption("p.adjust")) {
			P_ADJUST = cmd.getOptionValue("p.adjust");
			
			if(!P_ADJUST.equalsIgnoreCase("none") && !P_ADJUST.equalsIgnoreCase("Bonferroni") && !P_ADJUST.equalsIgnoreCase("BenjaminiHochberg")) {
				throw new IllegalArgumentException("Unsupported p-value adjustment method: '" + P_ADJUST + "'");
			}
		}
		
		if(cmd.hasOption("num-threads")) {
			NUM_THREADS = Integer.parseInt(cmd.getOptionValue("num-threads"));
		}
		
		if(cmd.hasOption("window-size")) {
			WINDOW_SIZE = Integer.parseInt(cmd.getOptionValue("window-size"));
		}
		
		if(cmd.hasOption("step-size")) {
			STEP_SIZE = Integer.parseInt(cmd.getOptionValue("step-size"));
		}
		
		if(STEP_SIZE > WINDOW_SIZE) {
			System.err.println("WARNING: computing windows with STEP_SIZE=" + STEP_SIZE + 
					" > WINDOW_SIZE=" + WINDOW_SIZE + ". Windows are not overlapping or book-ended!");
		}
		
		if(cmd.hasOption("min-window")) {
			MIN_WINDOW_SIZE = Integer.parseInt(cmd.getOptionValue("min-window"));
		}
		
		GENOME_SIZES = GetChromosomeSizes.get(genomeSizesFilename);
		
		System.out.println("WINDOW_SIZE: " + WINDOW_SIZE);
		System.out.println("STEP_SIZE: " + STEP_SIZE);
	}

	@SuppressWarnings("rawtypes")
	public static void windowCountsFishers(CommandLine cmd, ExecutorService threadPool, Map<String, Map<Integer, Double>> significantWindows,
			Map<String, IntervalTree<Double>> junctions, Map<String, Map<Integer, Set<Integer>>> junctions_map) throws IOException {
		System.out.println("MeRIPPeR: Window Counter & Fisher's Test...started");
		
		/*
		 * Readers and Output
		 */
		final SamReader sampleReader, controlReader;
		int sampleReadCounts, controlReadCounts;
		String[] sampleFilenames = cmd.getOptionValues("merip");
		String[] controlFilenames = cmd.getOptionValues("control");
		
		/*
		 * Counters
		 */
		Map<String, Map<Integer, Integer>> sample_counters = new HashMap<String, Map<Integer, Integer>>();
		Map<String, Map<Integer, Integer>> control_counters = new HashMap<String, Map<Integer, Integer>>();
		for(String s : GENOME_SIZES.keySet()) {
			sample_counters.put(s, new HashMap<Integer, Integer>());
			control_counters.put(s,  new HashMap<Integer, Integer>());
			significantWindows.put(s,  new HashMap<Integer, Double>());
			junctions.put(s,  new IntervalTree<Double>());
			junctions_map.put(s, new HashMap<Integer, Set<Integer>>());
		}
		
		/*
		 * Junction-specific counters and maps
		 */
		Map<String, Map<IntervalTree.Node<Double>, Integer>> junctions_merip_counter = null, 
				junctions_control_counter = null;
		if(cmd.hasOption("genes")) {
			JUNCTIONS = true;
			addJunctionsFromGenes(junctions, junctions_map, cmd.getOptionValue("genes"));
		}
		
		if(cmd.hasOption("junctions")) {
			JUNCTIONS = true;
			int minimum_coverage = 5;
			if(cmd.hasOption("junctions-min-coverage")) {
				minimum_coverage = Integer.parseInt(cmd.getOptionValue("junctions-min-coverage"));
			}
			
			addJunctionsFromFile(junctions, junctions_map, cmd.getOptionValue("junctions"), minimum_coverage);
		}

		if(JUNCTIONS) {
			junctions_merip_counter = new HashMap<String, Map<IntervalTree.Node<Double>, Integer>>();
			junctions_control_counter = new HashMap<String, Map<IntervalTree.Node<Double>, Integer>>();
			for(String chr : junctions.keySet()) {
				junctions_merip_counter.put(chr,  new HashMap<IntervalTree.Node<Double>, Integer>());
				junctions_control_counter.put(chr,  new HashMap<IntervalTree.Node<Double>, Integer>());
				
				for(Iterator<IntervalTree.Node<Double>> i = junctions.get(chr).iterator(); i.hasNext();) {
					IntervalTree.Node<Double> node = i.next();
					junctions_merip_counter.get(chr).put(node,  0);
					junctions_control_counter.get(chr).put(node,  0);
				}
			}
		}
		
		/* TODO: Implement other formats
		switch(format) {
			case BED:
				sampleReader = new Bed3IntervalReader(sampleFilename);
				controlReader = new Bed3IntervalReader(controlFilename);
				break;
			case SAM:
				sampleReader = new SamOrBamFileIntervalReader(sampleFilename);
				controlReader = new SamOrBamFileIntervalReader(controlFilename);
				break;
			case UNSUPPORTED:
			default:
				throw new IllegalArgumentException("Unsupported input file type");
		}
		*/
		
		sampleReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(sampleFilenames[0]));
		controlReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(controlFilenames[0]));
		
		System.err.println("Reading in sample reads from: " + sampleFilenames[0]);
		System.err.println("Reading in control reads from: " + controlFilenames[0]);
		
		Future<Integer> sampleReaderFuture, controlReaderFuture;
		if(junctions_merip_counter != null) {
			sampleReaderFuture = threadPool.submit(new ReadCounterJunctions(sampleReader, sample_counters, 
					junctions, junctions_merip_counter, WINDOW_SIZE, STEP_SIZE));
			controlReaderFuture = threadPool.submit(new ReadCounterJunctions(controlReader, control_counters,
					junctions, junctions_control_counter, WINDOW_SIZE, STEP_SIZE));
		} else {
			sampleReaderFuture = threadPool.submit(new ReadCounter(sampleReader, sample_counters, 
					WINDOW_SIZE, STEP_SIZE));
			controlReaderFuture = threadPool.submit(new ReadCounter(controlReader, control_counters,
					WINDOW_SIZE, STEP_SIZE));
		}
		
		try {
			// try to wait until the two threads finish
			long startTime = System.currentTimeMillis();
			sampleReadCounts = sampleReaderFuture.get();
			controlReadCounts = controlReaderFuture.get();
			long readTime = System.currentTimeMillis();
			
			System.err.println("Read " + sampleReadCounts + " from " + sampleFilenames[0]);
			System.err.println("Read " + controlReadCounts + " from " + controlFilenames[0]);
			System.err.println("Reading took " + ((readTime - startTime) / 1000) + " seconds");
			System.err.println("Read counting completed.");
			
			System.err.println("Starting Fishers Exact Tests.");
			// now submit fishers threads
			List<Future<Integer>> futures = new ArrayList<Future<Integer>>();
			for(String chr : GENOME_SIZES.keySet()) {
				futures.add(threadPool.submit(new FishersTestThread(sample_counters.get(chr), control_counters.get(chr), 
						sampleReadCounts, controlReadCounts, significantWindows.get(chr))));
				
				if(junctions.get(chr).size() > 0) {
					futures.add(threadPool.submit(new FishersTestJunctionThread(junctions_merip_counter.get(chr), 
							junctions_control_counter.get(chr), sampleReadCounts, controlReadCounts, junctions.get(chr))));
				}
			}
			
			// now "get" the Futures and cause a wait
			for(Future future : futures) {
				future.get();
			}
			
			long endTime = System.currentTimeMillis();
			System.err.println("Fishers Exact Test tasks completed in " + ((endTime - readTime) / 1000) + " seconds.");
		} catch (InterruptedException e) {
			e.printStackTrace();
		} catch (ExecutionException e) {
			e.printStackTrace();
		} finally {
			sampleReader.close();
			controlReader.close();
			System.err.println("MeRIPPeR PeakFinder: Window Counter & Fisher's Test...finished.");
		}
	}
	
	private static void addJunctionsFromGenes(Map<String, IntervalTree<Double>> junctions,
			Map<String, Map<Integer, Set<Integer>>>junctions_map,
			String genesFilename) throws FileNotFoundException {
		Scanner reader = new Scanner(new File(genesFilename));
		while(reader.hasNextLine()) {
			String[] line = reader.nextLine().split("\t");
			String chr = line[0];
			
			if(junctions.containsKey(chr)) {
				int txStart = Integer.parseInt(line[1]);
				// int txEnd = Integer.parseInt(line[2]);
				int nExons = Integer.parseInt(line[9]);
				
				String[] exonLengths = line[10].split(",");
				String[] exonStarts = line[11].split(",");
				
				for(int i = 0; i < nExons - 1; i++) {
					int start = txStart + Integer.parseInt(exonStarts[i]) + Integer.parseInt(exonLengths[i]) - 1;
					int end = txStart + Integer.parseInt(exonStarts[i + 1]);
					addJunction(junctions.get(chr), junctions_map.get(chr), start, end);
				}
			}
		}
		reader.close();
	}

	private static void addJunctionsFromFile(Map<String, IntervalTree<Double>> junctions,
			Map<String, Map<Integer, Set<Integer>>>junctions_map,
			String junctionsFilename, int minimumCoverage) throws FileNotFoundException {
		Scanner reader = new Scanner(new File(junctionsFilename));
		while(reader.hasNextLine()) {
			String[] line = reader.nextLine().split("\t");
			String chr = line[0];
			int start = Integer.parseInt(line[1]) - 2;
			int end = Integer.parseInt(line[2]);
			int cov = Integer.parseInt(line[6]);
			
			Map<Integer, Set<Integer>> chr_map = junctions_map.get(chr);
			if(cov >= minimumCoverage && chr_map != null) {
				addJunction(junctions.get(chr), chr_map, start, end);
			}
		}
		reader.close();
	}
	
	private static void addJunction(IntervalTree<Double> junctions_chr,
			Map<Integer, Set<Integer>> junctions_map_chr, int start, int end) {
		Set<Integer> set = junctions_map_chr.get(start);
		if(set == null) {
			set = new HashSet<Integer>();
			junctions_map_chr.put(start, set);
		}
		
		if(!set.contains(end)) {
			set.add(end);
			junctions_chr.put(start - WINDOW_SIZE,  start, 1.0);
			junctions_chr.put(end, end + WINDOW_SIZE - 1,  1.0);
		}
	}
	
	private static void pValueAdjustment(Map<String, Map<Integer, Double>> significantWindows) {
		
		long time0 = System.currentTimeMillis();

		// Perform p-value adjustment
		if(P_ADJUST  != null && !P_ADJUST.equalsIgnoreCase("none")) {
			long time1 = System.currentTimeMillis();
			System.err.println("Adjusting p-values using method: " + P_ADJUST);
			int N = 0;
			for(Map.Entry<String, Integer> entry : GENOME_SIZES.entrySet()) {
				N += (int) Math.ceil(entry.getValue() / STEP_SIZE);
			}
			System.err.println("Using N = " + N);
			
			if(P_ADJUST.equalsIgnoreCase("Bonferroni")) {
				ALPHA /= N;
			} else if(P_ADJUST.equalsIgnoreCase("BenjaminiHochberg")) {
				List<Double> pvalues = new ArrayList<Double>();
				for(String chr : GENOME_SIZES.keySet()) {
					for(double d : significantWindows.get(chr).values()) {
						pvalues.add(d);
					}
				}
				
				Collections.sort(pvalues);
				
				int k = 0;
				while(pvalues.get(k) <= ALPHA * k / N && k < pvalues.size() - 1) {
					k++;
				}
				
				ALPHA = pvalues.get(k - 1);
			}

			time1 = System.currentTimeMillis();
			System.err.println("Adjusting p-values completed in " + ((time1-time0) / 1000) + " seconds.");
			System.err.println("New alpha = " + ALPHA);
		}
	}

	private static void printFilteredWindows(ExecutorService threadPool, Map<String, Map<Integer, Double>> significantWindows,
			Map<String, IntervalTree<Double>> junctions, Map<String, Map<Integer, Set<Integer>>> junctions_map, 
			PrintStream out) throws InterruptedException, ExecutionException {
		
		List<String> chrs = new ArrayList<String>(significantWindows.keySet());
		Collections.sort(chrs);
		
		Map<String, List<Interval>> filteredWindows = new HashMap<String, List<Interval>>();
		List<Future<Integer>> futures = new ArrayList<Future<Integer>>();
		for(String chr : chrs) {
			filteredWindows.put(chr,  new ArrayList<Interval>());
			if(junctions.get(chr).size() > 0) {
				futures.add(threadPool.submit(new WindowFilterJunctionsThread(significantWindows.get(chr), 
						junctions.get(chr), junctions_map.get(chr), WINDOW_SIZE, GENOME_SIZES.get(chr), MIN_WINDOW_SIZE,
						filteredWindows.get(chr), ALPHA)));
			} else {
				futures.add(threadPool.submit(new WindowFilterThread(significantWindows.get(chr), 
						WINDOW_SIZE, GENOME_SIZES.get(chr), MIN_WINDOW_SIZE,
						filteredWindows.get(chr), ALPHA)));
			}
		}
		
		// cause a wait
		for(Future<Integer> f : futures) {
			f.get();
		}
		
		for(String chr : chrs) {
			for(Interval i : filteredWindows.get(chr)) {
				out.println(chr + "\t" + i.start + "\t" + i.end);
			}
		}
	}
	
	@SuppressWarnings("static-access")
	public static void buildOptions(Options options) {
		Option sample = OptionBuilder.withArgName("MeRIP-sample-file")
								  .isRequired(true)
								  .withLongOpt("merip")
								  .withDescription("MeRIP sample file")
								  .hasArg()
								  .withValueSeparator(',')
								  .create('m');
		options.addOption(sample);
		
		Option control = OptionBuilder.withArgName("control-file")
				  .isRequired(true)
				  .withLongOpt("control")
				  .withDescription("control file")
				  .hasArg()
				  .withValueSeparator(',')
				  .create('c');
		options.addOption(control);
		
		/*
		Option format = OptionBuilder.withArgName("file format")
							.withDescription("file format")
							.withLongOpt("format")
							.hasArg()
							.isRequired()
							.create('f');
		options.addOption(format);
		*/
		
		Option output = OptionBuilder.withArgName("output")
				.hasArg()
				.isRequired()
				.withLongOpt("output")
				.withDescription("output file")
				.create('o');
		options.addOption(output);
		
		Option genome = OptionBuilder.withArgName("genome-sizes")
				.hasArg()
				.isRequired()
				.withLongOpt("genome-sizes")
				.withDescription("Genome Chromosome Sizes")
				.create('g');
		options.addOption(genome);
		
		Option windowSize = OptionBuilder.withArgName("window size")
				.hasArg()
				.isRequired(false)
				.withLongOpt("window-size")
				.withDescription("Window Size")
				.create('w');
		options.addOption(windowSize);
		
		Option stepSize = OptionBuilder.withArgName("step size")
				.hasArg()
				.isRequired(false)
				.withLongOpt("step-size")
				.withDescription("Window Step Size")
				.create('s');
		options.addOption(stepSize);
		
		Option threads = OptionBuilder.withArgName("# threads")
							.hasArg()
							.isRequired(false)
							.withLongOpt("num-threads")
							.withDescription("# of threads")
							.create('t');
		options.addOption(threads);
		
		Option annotation = OptionBuilder.withArgName("Bed 12 genes file")
				.hasArg()
				.isRequired(false)
				.withLongOpt("genes")
				.withDescription("Genes annotation (RefSeq)")
				.create('r');
		options.addOption(annotation);
		
		Option junctions = OptionBuilder.withArgName("STAR Junctions file")
				.hasArg()
				.isRequired(false)
				.withLongOpt("junctions")
				.withDescription("splice junctions from STAR")
				.create('j');
		options.addOption(junctions);
		
		Option junctions_min = OptionBuilder.withArgName("Minimum coverage")
				.hasArg()
				.isRequired(false)
				.withLongOpt("junctions-min-coverage")
				.create('k');
		
		options.addOption(junctions_min);
		
		Option alpha = OptionBuilder.withArgName("threshold")
				  .withLongOpt("alpha")
				  .withDescription("p-value alpha threshold")
				  .hasArg()
				  .create('a');
		options.addOption(alpha);
		
		Option adjust = OptionBuilder.withArgName("BenjaminiHochberg|Bonferroni")
				  .withLongOpt("p.adjust")
				  .withDescription("p-value adjustment method")
				  .hasArg()
				  .create('p');
		options.addOption(adjust);
		
		Option minimum = OptionBuilder.withArgName("minimum window size")
				  .withLongOpt("min-window")
				  .withDescription("Filter out windows smaller than the minimum size")
				  .hasArg()
				  .create('n');
		options.addOption(minimum);

	}
}
