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

public class WindowSplitter {
	public static void main(String[] args) {
		long time0 = System.currentTimeMillis();
		Options options = new Options();
		buildOptions(options);
		CommandLineParser parser = new GnuParser();
		
		String inputFilename, outputFilename;
		int max = 200;
		
		try {
			CommandLine cmd = parser.parse(options, args, true);
			inputFilename = cmd.getOptionValue("input");
			outputFilename = cmd.getOptionValue("output");
			
			if(cmd.hasOption("max")) {
				max = Integer.parseInt(cmd.getOptionValue("max"));
			}
			
			Bed3IntervalReader reader = new Bed3IntervalReader(inputFilename);
			PrintStream out = new PrintStream(new File(outputFilename));
			
			// Read in windows
			System.out.println("Reading in windows from file: " + inputFilename);
			System.out.println("Outputting split windows to file: " + outputFilename);
			
			int countA = 0, countFS = 0;
			while(reader.hasNext()) {
				Bed3Interval interval = reader.next();
				
				if(interval.getLength() > max) {
					int s = interval.getStart();
					int e = interval.getEnd();
					int l = interval.getLength();
					int N = (int) Math.ceil(1.0 * l / max);
					int incr = (int) Math.ceil(1.0 * l / N);
					for(long i = s; i < e; i += incr) {
						out.println(interval.getChromosome() + "\t" + i + "\t" +
								Math.min(i+incr, e));
						countFS++;
					}
				} else {
					out.println(interval);
					countFS++;
				}
				countA++;
			}
			long time1 = System.currentTimeMillis();
			System.out.println("Reading in windows completed in " + ((time1-time0) / 1000) + " seconds.");
			System.out.println(countA + " windows split into " + countFS + " windows.");
		} catch (Throwable t) {
			t.printStackTrace();
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp("MeRIPPER.WindowSplitter", options, true);
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
		
		Option max = OptionBuilder.withArgName("window size [bp]")
				  .withLongOpt("max")
				  .withDescription("max")
				  .hasArg()
				  .create('l');
		options.addOption(max);
	}
}
