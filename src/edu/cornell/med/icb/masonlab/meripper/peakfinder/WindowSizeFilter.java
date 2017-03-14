package edu.cornell.med.icb.masonlab.meripper.peakfinder;

import java.io.File;
import java.io.PrintStream;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import edu.cornell.med.icb.masonlab.jenotator.io.input.Bed3IntervalReader;
import edu.cornell.med.icb.masonlab.jenotator.model.interval.Bed3Interval;

public class WindowSizeFilter {
	public static void main(String[] args) {
		long time0 = System.currentTimeMillis();
		Options options = new Options();
		buildOptions(options);
		CommandLineParser parser = new GnuParser();
		
		String inputFilename, outputFilename;
		int min = 100;
		
		try {
			CommandLine cmd = parser.parse(options, args, true);
			inputFilename = cmd.getOptionValue("input");
			outputFilename = cmd.getOptionValue("output");
			
			if(cmd.hasOption("min")) {
				min = Integer.parseInt(cmd.getOptionValue("min"));
			}
			
			Bed3IntervalReader reader = new Bed3IntervalReader(inputFilename);
			PrintStream out = new PrintStream(new File(outputFilename));
			
			// Read in windows
			System.out.println("Reading in windows from file: " + inputFilename);
			System.out.println("Outputting filtered windows to file: " + outputFilename);
			
			int countA = 0, countF = 0;
			while(reader.hasNext()) {
				Bed3Interval interval = reader.next();
				if(interval.getLength() >= min) {
					out.println(interval);
					countF++;
				}
				countA++;
			}
			long time1 = System.currentTimeMillis();
			System.out.println("Reading in windows completed in " + ((time1-time0) / 1000) + " seconds.");
			System.out.println(countF + "/" + countA + " windows passed minimum threshold.");
		} catch (Throwable t) {
			t.printStackTrace();
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp("MeRIPPER.WindowSizeFilter", options, true);
			System.exit(1); // so that make stops running
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
		
		Option output = OptionBuilder.withArgName("output file")
				  .isRequired(true)
				  .withLongOpt("output")
				  .withDescription("output filename")
				  .hasArg()
				  .create('o');
		options.addOption(output);
		
		Option min = OptionBuilder.withArgName("window size [bp]")
				  .withLongOpt("min")
				  .withDescription("minimum window size")
				  .hasArg()
				  .create('s');
		options.addOption(min);
	}
}
