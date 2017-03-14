package edu.cornell.med.icb.masonlab.meripper.util;

import htsjdk.samtools.util.IntervalTree;

import java.util.Map;
import java.util.concurrent.Callable;

import org.rfoundation.R.library.stats.FisherTest;

public class FishersTestJunctionThread implements Callable<Integer> {
	private int sampleReadCounts;
	private int controlReadCounts;
	private final IntervalTree<Double> junctions;
	protected final Map<IntervalTree.Node<Double>, Integer> sample, control;

	public FishersTestJunctionThread(Map<IntervalTree.Node<Double>, Integer> sample, Map<IntervalTree.Node<Double>, Integer> control,
			int sampleReadCounts, int controlReadCounts, IntervalTree<Double> junctions) {
		this.sample = sample;
		this.control = control;
		this.sampleReadCounts = sampleReadCounts;
		this.controlReadCounts = controlReadCounts;
		this.junctions = junctions;
	}

	@Override
	public Integer call() {
		int[][] fisherstable = new int[2][2];
		
		for(IntervalTree.Node<Double> node : junctions) {
			int sample_count = sample.get(node);
			int control_count = control.get(node);
						
			if(sample_count > 0 && (1.0 * sample_count / sampleReadCounts >= 1.0 * control_count / controlReadCounts)) {
				fisherstable[0][0] = sample_count;
				fisherstable[0][1] = sampleReadCounts - sample_count;
				fisherstable[1][0] = control_count;
				fisherstable[1][1] = controlReadCounts - control_count;
				
				double pvalue = FisherTest.test(fisherstable);
				node.setValue(pvalue);
			}
		}
		
		return new Integer(0);
	}
}
