package edu.cornell.med.icb.masonlab.meripper.peakfinder;

import htsjdk.samtools.SamReader;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import edu.cornell.med.icb.masonlab.jenotator.activity.GetChromosomeSizes;
import edu.cornell.med.icb.masonlab.meripper.util.FishersTestThread;
import edu.cornell.med.icb.masonlab.meripper.util.ReadCounter;

public class WindowCounter {
	public static int WINDOW_SIZE = 25;
	
	public enum FileFormat {
		BED,
		SAM,
		UNSUPPORTED;
		
		public static FileFormat parseFileFormat (String s) {
			if(s.equalsIgnoreCase("bed")) {
				return BED;
			} else if(s.equalsIgnoreCase("sam") || s.equalsIgnoreCase("bam")) {
				return SAM;
			} else {
				return UNSUPPORTED;
			}
		}
	}
	
	public static void main(String[] args) {
		CommandLineParser parser = new GnuParser();
		Options options = new Options();
		buildOptions(options);
		String sampleFilename, controlFilename, outputFilename, genomeSizesFilename;
		FileFormat format;
		int nthreads = 1;
		
		try {
			CommandLine cmd = parser.parse(options, args, true);
			sampleFilename = cmd.getOptionValue("merip");
			controlFilename = cmd.getOptionValue("control");
			outputFilename = cmd.getOptionValue("output");
			format = FileFormat.parseFileFormat(cmd.getOptionValue("format"));
			genomeSizesFilename = cmd.getOptionValue("genome-sizes");
			
			if(cmd.hasOption("num-threads")) {
				nthreads = Integer.parseInt(cmd.getOptionValue("num-threads"));
			}
			
			if(cmd.hasOption("window-size")) {
				WINDOW_SIZE = Integer.parseInt(cmd.getOptionValue("window-size"));
			}
			
			run(genomeSizesFilename, sampleFilename, controlFilename, outputFilename, format, nthreads);
		} catch(ParseException e) {
			e.printStackTrace();
			System.err.println("Improper arguments");
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp("MeRIPPER.PeakFinder", options, true);
			System.exit(1); // so that make stops running
		} catch (IOException e) {
			e.printStackTrace();
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp("MeRIPPER.PeakFinder", options, true);
			System.exit(1); // so that make stops running
		} catch (Throwable t) {
			t.printStackTrace();
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp("MeRIPPER.PeakFinder", options, true);
			System.exit(1); // so that make stops running
		}
	}
	
	@SuppressWarnings("rawtypes")
	public static void run(String genomeSizesFilename, String sampleFilename, String controlFilename, String outputFilename, FileFormat format, int nthreads) throws IOException {
		System.out.println("MeRIPPeR PeakFinder...");
		ExecutorService threadPool = Executors.newFixedThreadPool(nthreads);
		Map<String, int[]> sample_counters = new HashMap<String, int[]>();
		Map<String, int[]> control_counters = new HashMap<String, int[]>();
		Map<String, Integer> sizes = GetChromosomeSizes.get(genomeSizesFilename);
		final SamReader sampleReader, controlReader;
		PrintStream windowPrintStream = new PrintStream(new File(outputFilename));
		int sampleReadCounts, controlReadCounts;
		
		// initialize counters
		for(Map.Entry<String, Integer> entry : sizes.entrySet()) {
			sample_counters.put(entry.getKey(), new int[(int) Math.ceil(entry.getValue() / WINDOW_SIZE)]);
			control_counters.put(entry.getKey(), new int[(int) Math.ceil(entry.getValue() / WINDOW_SIZE)]);
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
		
		sampleReader = new SAMFileReader(new File(sampleFilename));
		controlReader = new SAMFileReader(new File(controlFilename));
		sampleReader.setValidationStringency(ValidationStringency.SILENT);
		controlReader.setValidationStringency(ValidationStringency.SILENT);
		
		System.err.println("Reading in sample reads from: " + sampleFilename);
		System.err.println("Reading in control reads from: " + controlFilename);
		Future<Integer> sampleReaderFuture = threadPool.submit(new ReadCounter(sampleReader, sample_counters, WINDOW_SIZE));
		Future<Integer> controlReaderFuture = threadPool.submit(new ReadCounter(controlReader, control_counters, WINDOW_SIZE));
		
		try {
			// try to wait until the two threads finish
			long startTime = System.currentTimeMillis();
			sampleReadCounts = sampleReaderFuture.get();
			controlReadCounts = controlReaderFuture.get();
			long readTime = System.currentTimeMillis();
			
			System.err.println("Read " + sampleReadCounts + " from " + sampleFilename);
			System.err.println("Read " + controlReadCounts + " from " + controlFilename);
			System.err.println("Reading took " + ((readTime - startTime) / 1000) + " seconds");
			System.err.println("Read counting completed.");
			
			System.err.println("Starting Fishers Exact Tests.");
			// now submit fishers threads
			List<Future> futures = new ArrayList<Future>();
			for(String chr : sizes.keySet()) {
				futures.add(threadPool.submit(new FishersTestThread(chr, sample_counters.get(chr), control_counters.get(chr), sampleReadCounts, controlReadCounts, windowPrintStream)));
			}
			
			// now "get" the Futures and cause a wait
			for(Future future : futures) {
				future.get();
			}
			threadPool.shutdown();
			while(!threadPool.isShutdown()) {
				Thread.sleep(1000);
			}
			
			long endTime = System.currentTimeMillis();
			System.err.println("Fishers Exact Test tasks completed in " + ((endTime - readTime) / 1000) + " seconds.");
			
		} catch (InterruptedException e) {
			e.printStackTrace();
		} catch (ExecutionException e) {
			e.printStackTrace();
		} finally {
			windowPrintStream.close();
			System.err.println("All tasks completed.");
		}
	}
	
	@SuppressWarnings("static-access")
	public static void buildOptions(Options options) {
		Option sample = OptionBuilder.withArgName("MeRIP-sample-file")
								  .isRequired(true)
								  .withLongOpt("merip")
								  .withDescription("MeRIP sample file")
								  .hasArg()
								  .create('m');
		options.addOption(sample);
		
		Option control = OptionBuilder.withArgName("control-file")
				  .isRequired(true)
				  .withLongOpt("control")
				  .withDescription("control file")
				  .hasArg()
				  .create('c');
		options.addOption(control);
		
		Option format = OptionBuilder.withArgName("file format")
							.withDescription("file format")
							.withLongOpt("format")
							.hasArg()
							.isRequired()
							.create('f');
		options.addOption(format);
		
		Option output = OptionBuilder.withArgName("output")
				.hasArg()
				.isRequired()
				.withLongOpt("output")
				.withDescription("output file")
				.create('o');
		options.addOption(output);
		
		Option genome = OptionBuilder.withArgName("genome.sizes")
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
		
		Option threads = OptionBuilder.withArgName("# threads")
							.hasArg()
							.isRequired(false)
							.withLongOpt("num-threads")
							.withDescription("# of threads")
							.create('t');
		options.addOption(threads);
	}
}
