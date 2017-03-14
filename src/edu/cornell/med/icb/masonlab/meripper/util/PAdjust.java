package edu.cornell.med.icb.masonlab.meripper.util;

import java.util.Arrays;
import java.util.Comparator;

import edu.cornell.med.icb.masonlab.jenotator.model.interval.Bed6Interval;

public class PAdjust {
	public static void adjust(String method, double alpha, Bed6Interval[] windows, int N) {		
		if(method.equalsIgnoreCase("Bonferroni")) {
			for(Bed6Interval window : windows) {
				window.setScore(window.getScore() / N);
			}
		} else if(method.equalsIgnoreCase("BenjaminiHochberg")) {
			// now sort the pvalue indexes
			Arrays.sort(windows, new Comparator<Bed6Interval> () {
				@Override
				public int compare(Bed6Interval o1, Bed6Interval o2) {
					return Double.compare(o1.getScore(), o2.getScore());
				}});

			int i = 1;
			for(Bed6Interval window : windows) {
				window.setScore(1.0d * N * window.getScore() / i);
				i++;
			}
			
			// now resort
			Arrays.sort(windows);
		} else {
			throw new IllegalArgumentException("P-value adjustment method '" + method + "' not found.");
		}
	}
}
