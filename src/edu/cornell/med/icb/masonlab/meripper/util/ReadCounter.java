package edu.cornell.med.icb.masonlab.meripper.util;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;

import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;

public class ReadCounter implements Callable<Integer> {
	protected final SamReader reader;
	protected final Map<String, Map<Integer, Integer>> counters;
	protected final SAMRecordIterator iterator;
	protected final int WINDOW_SIZE;
	protected final int STEP_SIZE;

	
	public ReadCounter(final SamReader reader, Map<String, Map<Integer, Integer>> counters, int window_size, int step_size) {
		this.counters = counters;
		this.reader = reader;
		this.iterator = reader.iterator();
		this.WINDOW_SIZE = window_size;
		this.STEP_SIZE = step_size;
	}

	@Override
	public Integer call() {
		int count = 0;
		
		Set<Integer> windows = new HashSet<Integer>();
		while(iterator.hasNext()) {
			SAMRecord record = iterator.next();
			String chr = record.getReferenceName();
			if(!record.getReadUnmappedFlag() && counters.containsKey(chr)) {
				for(AlignmentBlock block : record.getAlignmentBlocks()) {
					int start = Math.max(0, ((block.getReferenceStart() - 1) - (WINDOW_SIZE - STEP_SIZE)) / STEP_SIZE);
					int end = (block.getReferenceStart() - 1 + block.getLength()) / STEP_SIZE;
					
					for(int i = start; i <= end; i++) {
						windows.add(i * STEP_SIZE);
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
				
				windows.clear();
				
				count++;
			}
		}
		
		return count;
	}
}