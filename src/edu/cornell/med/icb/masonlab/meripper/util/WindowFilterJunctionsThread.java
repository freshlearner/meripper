package edu.cornell.med.icb.masonlab.meripper.util;

import htsjdk.samtools.util.IntervalTree;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;

public class WindowFilterJunctionsThread implements Callable<Integer>{
	private final Map<Integer, Double> pvalues;
	private final IntervalTree<Double> junctions;
	private final int CHR_MAX;
	private final int WINDOW_MIN;
	private final int WINDOW_SIZE;
	private final double ALPHA;
	private final List<Interval> final_list;
	private final Map<Integer, Set<Integer>> junctions_map;

	public WindowFilterJunctionsThread(Map<Integer, Double> pvalues, IntervalTree<Double> junctions,
			Map<Integer, Set<Integer>> junctions_map,
			int window_size, int chr_max, int window_min, List<Interval> list, double alpha) {
		this.pvalues = pvalues;
		this.junctions = junctions;
		this.CHR_MAX = chr_max;
		this.WINDOW_MIN = window_min;
		this.WINDOW_SIZE = window_size;
		this.final_list = list;
		this.ALPHA = alpha;
		this.junctions_map = junctions_map;
	}
	
	@Override
	public Integer call() throws Exception {
		IntervalTree<Integer> tree = new IntervalTree<Integer>();
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
					// this window has a significant p-value
					int start = i;
					int end = Math.min(start + WINDOW_SIZE, CHR_MAX);
					
					if(prev_end >= start) {
						// merged case
						prev_end = Math.max(prev_end, end);
					} else {
						if(prev_end - prev_start >= WINDOW_MIN) {
							tree.put(prev_start, prev_end - 1, 1);
						} else {
							tree.put(prev_start, prev_end - 1, 0);
						}
						
						prev_start = start;
						prev_end = end;
					}
				}
				
				if(prev_end - prev_start >= WINDOW_MIN) {
					tree.put(prev_start, prev_end - 1, 1);
				} else {
					tree.put(prev_start, prev_end - 1, 0);
				}
			}
		}
		
		for(IntervalTree.Node<Double> node : junctions) {
			if(node.getValue() <= ALPHA) {
				tree.put(node.getStart(), node.getEnd(), 0);
			}
		}
		
		for(IntervalTree.Node<Double> node : junctions) {
			if(node.getValue() <= ALPHA) {
				List<IntervalTree.Node<Integer>> list = new ArrayList<IntervalTree.Node<Integer>>();
				int prev_start = node.getStart();
				int prev_end = node.getEnd();
				Iterator<IntervalTree.Node<Integer>> i = tree.overlappers(prev_start, prev_end);
				while(i.hasNext()) {
					IntervalTree.Node<Integer> interval = i.next();
					list.add(interval);
					prev_start = Math.min(interval.getStart(), prev_start);
					prev_end = Math.max(interval.getEnd(), prev_end);
				}
				
				if(prev_end - prev_start >= WINDOW_MIN) {
					for(IntervalTree.Node<Integer> k : list) {
						k.setValue(1);
					}
				}
			}
		}
		
		for(int i : junctions_map.keySet()) {
			for(int j : junctions_map.get(i)) {
				List<IntervalTree.Node<Integer>> list = new ArrayList<IntervalTree.Node<Integer>>();
				int prev_start = -1, prev_end = -1, length = 0;
				Iterator<IntervalTree.Node<Integer>> it = tree.overlappers(i - WINDOW_SIZE, i - 1);
				if(it.hasNext()) {
					IntervalTree.Node<Integer> node = it.next();
					list.add(node);
					prev_start = node.getStart();
					prev_end = node.getEnd();
					
					while(it.hasNext()) {
						node = it.next();
						list.add(node);
						prev_start = Math.min(prev_start,  node.getStart());
						prev_end = Math.max(prev_end, node.getEnd());
					}
					
					length += prev_end - prev_start;
				}
				
				it = tree.overlappers(j, j + WINDOW_SIZE - 1);
				if(it.hasNext()) {
					IntervalTree.Node<Integer> node = it.next();
					prev_start = node.getStart();
					prev_end = node.getEnd();
					
					while(it.hasNext()) {
						node = it.next();
						list.add(node);
						prev_start = Math.min(prev_start,  node.getStart());
						prev_end = Math.max(prev_end, node.getEnd());
					}
					
					length += prev_end - prev_start;
				}
				
				if(length >= WINDOW_MIN) {
					for(IntervalTree.Node<Integer> node : list) {
						node.setValue(1);
					}
				}
			}
		}
		
		Iterator<IntervalTree.Node<Integer>> i = tree.iterator();
		if(i.hasNext()) {
			IntervalTree.Node<Integer> interval = i.next();
			
			while(interval.getValue() < 1 && i.hasNext()) {
				interval = i.next();
			}
			
			if(interval.getValue() > 0) {
				int prev_start = interval.getStart();
				int prev_end = interval.getEnd();
				
				while(i.hasNext()) {
					interval = i.next();
					if(interval.getValue() > 0) {
						if(prev_end >= interval.getStart()) {
							// merged case
							prev_start = Math.min(prev_start, interval.getStart());
							prev_end = Math.max(prev_end, interval.getEnd());
						} else {
							final_list.add(new Interval(prev_start, prev_end + 1));
							
							prev_start = interval.getStart();
							prev_end = interval.getEnd();
						}
					}
				}

				final_list.add(new Interval(prev_start, prev_end + 1));
			}
		}
		
		return 0;
	}

}
