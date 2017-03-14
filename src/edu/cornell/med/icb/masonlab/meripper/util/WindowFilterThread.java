package edu.cornell.med.icb.masonlab.meripper.util;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;

public class WindowFilterThread implements Callable<Integer>{
	private final Map<Integer, Double> pvalues;
	private final int CHR_MAX;
	private final int WINDOW_MIN;
	private final int WINDOW_SIZE;
	private final double ALPHA;
	private final List<Interval> final_list;

	public WindowFilterThread(Map<Integer, Double> pvalues, int window_size, int chr_max, int window_min, List<Interval> list, double alpha) {
		this.pvalues = pvalues;
		this.CHR_MAX = chr_max;
		this.WINDOW_MIN = window_min;
		this.WINDOW_SIZE = window_size;
		this.ALPHA = alpha;
		this.final_list = list;
	}
	
	@Override
	public Integer call() throws Exception {
		if(pvalues.size() > 0) {
			List<Integer> windows_indexes = new ArrayList<Integer>(pvalues.keySet());
			for(Iterator<Integer> i = windows_indexes.iterator(); i.hasNext(); ) {
				if(pvalues.get(i.next()) > ALPHA) {
					i.remove();
				}
			}
			Collections.sort(windows_indexes);
			
			if(windows_indexes.size() > 0) {
				int prev_start = windows_indexes.get(0);
				int prev_end = Math.min(prev_start + WINDOW_SIZE, CHR_MAX);
				
				for(int j = 1; j < windows_indexes.size(); j++) {
					int i = windows_indexes.get(j);
					int start = i;
					int end = Math.min(start + WINDOW_SIZE, CHR_MAX);
					
					if(prev_end >= start) {
						// merged case
						prev_end = Math.max(prev_end, end);
					} else {
						if(prev_end - prev_start >= WINDOW_MIN) {
							final_list.add(new Interval(prev_start, prev_end));
						}
						
						prev_start = start;
						prev_end = end;
					}
				}
				
				if(prev_end - prev_start >= WINDOW_MIN) {
					final_list.add(new Interval(prev_start, prev_end));
				}
			}
		}
		
		return 0;
	}

}
