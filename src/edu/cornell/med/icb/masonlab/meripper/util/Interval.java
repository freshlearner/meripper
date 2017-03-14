package edu.cornell.med.icb.masonlab.meripper.util;

public class Interval implements Comparable<Interval> {
	public final int start;
	public final int end;
	
	public Interval(int start, int end) {
		this.start = start;
		this.end = end;
	}

	@Override
	public int compareTo(Interval j) {		
		return this.start - j.start;
	}
}
