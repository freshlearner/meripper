package edu.cornell.med.icb.masonlab.meripper.util;

import java.util.Map;
import java.util.concurrent.Callable;

import org.rfoundation.R.library.stats.FisherTest;

public class FishersTestThread implements Callable<Integer> {
	private final Map<Integer, Integer> sample;
	private final Map<Integer, Integer> control;
	private final int sampleReadCounts;
	private final int controlReadCounts;
	private final Map<Integer, Double> significantWindows;
	
	public FishersTestThread(Map<Integer, Integer> sample, Map<Integer, Integer> control, int sampleReadCounts, int controlReadCounts, 
			Map<Integer, Double> significantWindows) {
		this.sample = sample;
		this.control = control;
		this.sampleReadCounts = sampleReadCounts;
		this.controlReadCounts = controlReadCounts;
		this.significantWindows = significantWindows;
	}

	@Override
	public Integer call() {
		int[][] fisherstable = new int[2][2];
		
		for(int window : sample.keySet()) {
			int sample_count = 0;
			int control_count = 0;
			if(sample.containsKey(window)) {
				sample_count = sample.get(window);
			}
			if(control.containsKey(window)) {
				control_count = control.get(window);
			}
			
			if(sample_count > 0 && (1.0 * sample_count / sampleReadCounts >= 1.0 * control_count / controlReadCounts)) {
				fisherstable[0][0] = sample_count;
				fisherstable[0][1] = sampleReadCounts - sample_count;
				fisherstable[1][0] = control_count;
				fisherstable[1][1] = controlReadCounts - control_count;
				
				double pvalue = FisherTest.test(fisherstable);
				if(pvalue <= 0.05) {
					/* NOTE: ASSUME ONLY ONE THREAD WILL ACCESS THIS AT A TIME */
					significantWindows.put(window, pvalue);
				}
			}
		}
		
		return new Integer(0);
	}
}
