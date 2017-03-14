package edu.cornell.med.icb.masonlab.meripper.peakfinder;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import edu.cornell.med.icb.masonlab.jenotator.activity.GetChromosomeSizes;
import edu.cornell.med.icb.masonlab.jenotator.io.input.Bed6IntervalReader;
import edu.cornell.med.icb.masonlab.jenotator.model.interval.Bed6Interval;
import edu.cornell.med.icb.masonlab.meripper.util.PAdjust;

public class PValueWindowFilter {
	public static void main(String[] args) {
		Options options = new Options();
		buildOptions(options);
		CommandLineParser parser = new GnuParser();
		
		String inputFilename, outputFilename, genomeSizesFilename, pAdjust = null;
		double alpha = 0.05;
		int window_size = 25;
		
		try {
			CommandLine cmd = parser.parse(options, args, true);
			inputFilename = cmd.getOptionValue("input");
			outputFilename = cmd.getOptionValue("output");
			genomeSizesFilename = cmd.getOptionValue("genome-sizes");
			Map<String, Integer> sizes = GetChromosomeSizes.get(genomeSizesFilename);
			
			if(cmd.hasOption("alpha")) {
				alpha = Double.parseDouble(cmd.getOptionValue("alpha"));
			}
			
			if(cmd.hasOption("p.adjust")) {
				pAdjust = cmd.getOptionValue("p.adjust");
			}
			
			if(cmd.hasOption("window-size")) {
				window_size = Integer.parseInt(cmd.getOptionValue("window-size"));
			}
			
			Bed6IntervalReader reader = new Bed6IntervalReader(inputFilename);
			PrintStream out = new PrintStream(new File(outputFilename));
			List<Bed6Interval> windowsList = new ArrayList<Bed6Interval>();
			
			System.err.println("Reading in windows from file: " + inputFilename);
			long time0s = System.currentTimeMillis();
			while(reader.hasNext()) {
				windowsList.add(reader.next());
			}
			long time0f = System.currentTimeMillis();
			System.err.println("Reading of windows completed in: " + ((time0f-time0s)/1000) + " seconds.");
			Bed6Interval[] windows = (Bed6Interval[]) windowsList.toArray(new Bed6Interval[1]);
			
			// Perform p-value adjustment
			if(pAdjust != null && !pAdjust.equals("none")) {
				long time1 = System.currentTimeMillis();
				System.err.println("Adjusting p-values using method: " + pAdjust);
				int N = 0;
				for(Map.Entry<String, Integer> entry : sizes.entrySet()) {
					N += (int) Math.ceil(entry.getValue() / window_size);
				}
				System.err.println("Using N = " + N);
				PAdjust.adjust(pAdjust, alpha, windows, N);
				long time2 = System.currentTimeMillis();

				time1 = System.currentTimeMillis();
				System.err.println("Adjusting p-values completed in " + ((time2-time1) / 1000) + " seconds.");
				System.err.println("New alpha = " + alpha);
			}
			
			int counter = 0, counterA = 0;
			
			System.err.println("Printing out significant windows to file: " + outputFilename);
			long time1 = System.currentTimeMillis();
			for(Bed6Interval window : windows) {
				if(window.getScore() <= alpha) {
					out.println(window);
					counter++;
				}
				counterA++;
			}
			long time2 = System.currentTimeMillis();
			out.close();
			System.err.println("Printing of significant windows completed in " + ((time2-time1)/1000) + " seconds.");
			System.err.println(counter + " / " + counterA + " windows significant.");
		} catch (Throwable t) {
			t.printStackTrace();
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp("MeRIPPER.PValueWindowFilter", options, true);
			System.exit(1);
		}
	}
	
	@SuppressWarnings("static-access")
	public static void buildOptions(Options options) {
		Option input = OptionBuilder.withArgName("input")
								  .isRequired(true)
								  .withLongOpt("input")
								  .withDescription("input filename")
								  .hasArg()
								  .create('i');
		options.addOption(input);
		
		Option output = OptionBuilder.withArgName("output")
				  .isRequired(true)
				  .withLongOpt("output")
				  .withDescription("output filename")
				  .hasArg()
				  .create('o');
		options.addOption(output);
		
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
	}
}
