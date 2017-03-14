package edu.cornell.med.icb.masonlab.meripper.util;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.IntervalTree;

import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;

public class ReadCounterJunctions implements Callable<Integer> {
	protected final SamReader reader;
	protected final Map<String, Map<Integer, Integer>> counters;
	protected final SAMRecordIterator iterator;
	protected final int WINDOW_SIZE;
	protected final int STEP_SIZE;
	
	protected final Map<String, IntervalTree<Double>> junctions;
	protected final Map<String, Map<IntervalTree.Node<Double>, Integer>> junction_counters;
	
	public ReadCounterJunctions(final SamReader reader, Map<String, Map<Integer, Integer>> counters, Map<String, IntervalTree<Double>> junctions,
			Map<String, Map<IntervalTree.Node<Double>, Integer>> junction_counters, int window_size, int step_size) {
		this.counters = counters;
		this.reader = reader;
		this.iterator = reader.iterator();
		this.junctions = junctions;
		this.junction_counters = junction_counters;
		this.WINDOW_SIZE = window_size;
		this.STEP_SIZE = step_size;
	}

	@Override
	public Integer call() {
		int count = 0;
		
		while(iterator.hasNext()) {
			SAMRecord record = iterator.next();
			if(!record.getReadUnmappedFlag() && counters.containsKey(record.getReferenceName())) {
				String chr = record.getReferenceName();
				Set<Integer> windows = new HashSet<Integer>();
				Set<IntervalTree.Node<Double>> overlaps = new HashSet<IntervalTree.Node<Double>>();
				for(AlignmentBlock block : record.getAlignmentBlocks()) {
					int start = Math.max(0, ((block.getReferenceStart() - 1) - (WINDOW_SIZE - STEP_SIZE)) / STEP_SIZE);
					int end = (block.getReferenceStart() - 1 + block.getLength()) / STEP_SIZE;
					
					for(int i = start; i <= end; i++) {
						windows.add(i * STEP_SIZE);
					}
					
					for(Iterator<IntervalTree.Node<Double>> i = junctions.get(chr).overlappers(
								block.getReferenceStart() - 1, 
								block.getReferenceStart() - 1 + block.getLength() - 1);
							i.hasNext();) {
						overlaps.add(i.next());
					}
				}
				
				// now increment the counters
				// this prevents a read that's spliced twice onto the same exon from being double counted on that exon
				for(int window : windows) {
					if(counters.get(chr).containsKey(window)) {
						counters.get(chr).put(window, counters.get(chr).get(window) + 1);
					} else {
						counters.get(chr).put(window,  1);
					}
				}
				
				for(IntervalTree.Node<Double> node : overlaps) {
					junction_counters.get(chr).put(node,
							junction_counters.get(chr).get(node) + 1);
				}
				
				count++;
			}
		}
		
		return count;
	}
}