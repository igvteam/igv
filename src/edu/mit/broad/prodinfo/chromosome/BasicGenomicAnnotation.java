/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package edu.mit.broad.prodinfo.chromosome;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Stack;

//import edu.mit.broad.prodinfo.gene.GeneAnnotation;
import edu.mit.broad.prodinfo.genomicplot.GenomicAnnotation;
import edu.mit.broad.prodinfo.genomicplot.TwoSubjectAnnotation;
import edu.mit.broad.prodinfo.sequence.Sequence;

public class BasicGenomicAnnotation implements GenomicAnnotation {
	private int start;
	private int end;
	private String name;
	//private String displayName;
	private String orientation;
	private String chromosome;
	private Sequence sequence;
	private int fivePrimerBuffer;
	private int threePrimerBuffer;
	private String id;
	private double score;
	
	//This is a hack to speed up sequencial search and operations that involve a list
	int lastBreakInIndex = 0;
	
	protected BasicGenomicAnnotation() {
		super();
	}
	public BasicGenomicAnnotation(String name) {
		super();
		this.name = name;
		//this.displayName = name;
	}
	
	public BasicGenomicAnnotation(String name, String chr, int start, int end) {
		this(name);
		setChromosome(chr);
		setStart(start);
		setEnd(end);
	}
	
	public BasicGenomicAnnotation(GenomicAnnotation ga) {
		start = ga.getStart();
		end   = ga.getEnd();
		name  = ga.getName();
		orientation = ga.getOrientation();
		score = ga.getScore();
		chromosome = ga.getChromosome() != null ? ga.getChromosome().intern() : null;
		sequence = ga.getSequence();
		fivePrimerBuffer = ga.getFivePrimeBases();
		threePrimerBuffer = ga.getThreePrimeBases();
		id = ga.getId();
	}
	
	/**
	 * A basic genomic annotation may be built from a raw string array
	 * where the first entry is the chromosome, the second the start and
	 * the third the end.
	 */
	public BasicGenomicAnnotation(String [] info) {
		super();
		//? info[0] : info[0] + "_" + info[1] + "-" + info[2])
		if(info[0].contains(":") ) {
			String [] firstSplit = info[0].split(":");
			setChromosome(firstSplit[0]);
			String [] secondSplit = firstSplit[1].split("-");
			setStart(Integer.parseInt(secondSplit[0]));
			setEnd(Integer.parseInt(secondSplit[1]));
		} else {
			setChromosome(info[0].substring(3));
			setStart(Integer.parseInt(info[1]));
			setEnd(Integer.parseInt(info[2]));
			if(info.length > 3) {
				setId(info[3]);
				setName(info[3]);
			} else {
				setName(info[0] + "_" + info[1] + "-" + info[2]);
			}
		}
	}
	
	public boolean equals(Object o) {
		if(o == null || !this.getClass().equals(o.getClass())) {
			return false;
		}
		
		BasicGenomicAnnotation other = (BasicGenomicAnnotation) o;
		
		if(chromosome != null && !chromosome.equals(other.getChromosome())) {
			return false;
		}
		
		if(name != null && !name.equals(other.getName())) {
			return false;
		}
		
		return start == other.getStart() && end == other.getEnd(); 
	}
	
	public int hashCode() {
		int hash = 0;

		if(chromosome != null) {
			char [] chrChars = chromosome.toCharArray();
			for(int i = 0; i < chrChars.length; i++) {
				hash += chrChars[i];
			}
			hash = hash << 31;
		}
		
		if (name != null) {
			hash += name.hashCode() << 30;
		}
		
		hash += (start + end)/2;
		return hash;
	}
	
	public long getMiddle() {
		return Math.round( (getStart() + getEnd())/(float)2);
	}

	public int getLength() { return getEnd() - getStart(); }
	
	public String toString() {
		return getLocationString();
	}
	
	public String getChromosomeString() {
		if(getChromosome() == null) {
			return null;
		}

		String str = "";
		if(getChromosome().length() < 3 ) {
			str = "chr";
		}
		
		return str + getChromosome();
	}
	public String getLocationString() {
		return getChromosomeString()+":"+getStart()+"-"+getEnd() ;
	}
	protected void setStart(String data) {
		this.start = Integer.parseInt(data);
	}
	
	public void setStart(int start) {
		this.start = start;
	}
	
	protected void setEnd(String data) {
		this.end = Integer.parseInt(data);
	}
	
	public void setEnd(int end) {
		this.end = end;
	}
	
	public int getStart() {
		return start;
	}

	public int getEnd() {
		return end;
	}
	
	public double getScore() { return score;}
	public void setScore(double score) { this.score = score;}
	
	public void setBoundariesFromAnnoations(List<? extends GenomicAnnotation> annotations) {
		Iterator<? extends GenomicAnnotation> it = annotations.iterator();
		while(it.hasNext()) {
			GenomicAnnotation annot = it.next();
			start = start > 0 ? Math.min(start, annot.getStart()) : annot.getStart();
			end   = Math.max(end, annot.getEnd());
		}
	}

	/**
	 * Change the current instance to represent its intersection
	 * with the provided annotation if they overlap
	 * @param other
	 */
	public void takeIntersection(GenomicAnnotation other) {
		if(getStart() < other.getEnd() && getEnd() > other.getStart()) {
			setStart(Math.max(getStart(),other.getStart()));
			setEnd(Math.min(getEnd(), other.getEnd()));
		}
	}
	
	/**
	 * public void 
	 */
	
	public int getDistanceTo(GenomicAnnotation other) {
		int dist = 0;
		if(!getChromosome().equals(other.getChromosome())) {
			dist = 1000000000;
		}else if(!overlaps(other)) {
			if(getStart() > other.getEnd()) {
				dist = getStart() - other.getEnd();
			} else {
				dist = other.getStart() - getEnd();
			}
		}
		return dist;
	}
	
	public String getId() {
		return id;
	}
	public void setId(String id) {
		this.id = id.intern();
	}
	
	/**
	 * Change the current instance to represent its union
	 * with the provided annotation
	 */
	public void takeUnion(GenomicAnnotation other) {
		if(overlaps(other)) {
			setStart(Math.min(getStart(),other.getStart()));
			setEnd(Math.max(getEnd(), other.getEnd()));
		}
	}
	
	/**
	 * Change the current instance to represent the difference
	 * of the instance with the other annotation
	 * @param other
	 */
	public void reduceToDifference(GenomicAnnotation other) {
		if(overlaps(other)) {
			//System.out.println("Reduced annotation! " + this);
			if(getStart() < other.getStart()) {
				setEnd(Math.min(getEnd(), other.getStart() - 1));
			} else {
				setStart(Math.max(getStart(), other.getEnd() + 1));
			}
		}
	}
	
	/**
	 * Fragments the current annotation if it overlaps the provided one.
	 * it returns an empty list if this annotation does not ovelap the
	 * one passed in.
	 * 
	 * @param GenomicAnnotation annotation to disect the current one
	 * @return List<GenomicAnnotation> of the disected annotations, an empty list is returned if 
	 * 			the annotations do not overlap.
	 */
	public List<GenomicAnnotation> disect(GenomicAnnotation a) {
		List<GenomicAnnotation> disection = new ArrayList<GenomicAnnotation>();
		if(overlaps(a)) {
			disection.addAll(this.minus(a));
			//System.out.println("minus resulted in " + disection);
			GenomicAnnotation copy = new BasicGenomicAnnotation(this);
			copy.takeIntersection(a);
			//System.out.println("and intersection in " + copy);
			disection.add(copy);
			
			Collections.sort(disection);
		}
		//System.out.println("and disction is " + disection);
		return disection;
	}
	
	
	/**
	 * Fragments the current annotation if it overlaps the provided ones.
	 * It returns a list with one component (this annotation) if no annotation 
	 * in the provided list overlaps the discted annotaion.
	 * 
	 * @param List<GenomicAnnotation> <b>sorted</b> annotations with which to disect the current one.
	 * @return List<GenomicAnnotation> of the disected annotations, a list just this annotation is returned if 
	 * 			the annotations do not overlap.
	 */
	public List<GenomicAnnotation> disect(List<? extends GenomicAnnotation> disectors) {
		List<GenomicAnnotation> disection = new ArrayList<GenomicAnnotation>();
		disection.add(this);
		Iterator<? extends GenomicAnnotation> annotIt = disectors.iterator();
		while(annotIt.hasNext()) {
			GenomicAnnotation a= annotIt.next();
			if(a.getStart() > getEnd()) {
				break;
			}else if(!overlaps(a)) {
				continue;
			}
			List<GenomicAnnotation> newDisection = new ArrayList<GenomicAnnotation>(disection.size());
			Iterator<GenomicAnnotation> oldDisectionIt = disection.iterator();
			while(oldDisectionIt.hasNext()) {
				GenomicAnnotation disectionComponent = oldDisectionIt.next();
				if(disectionComponent.overlaps(a)) {
					newDisection.addAll(disectionComponent.disect(a));
				} else {
					newDisection.add(disectionComponent);
				}
			} 
			disection = newDisection;
		}
		return disection;
	}
	
	public void stitchTo(GenomicAnnotation other) {
		setStart(Math.min(getStart(), other.getStart()));
		setEnd(Math.max(getEnd(), other.getEnd()));
	}
	
	public boolean isFlankedBy(TwoSubjectAnnotation twoSubjectAnnotation, int buffer) {
		GenomicAnnotation left = new BasicGenomicAnnotation("left");
		left.setStart(getStart() - buffer);
		left.setEnd(getStart() + buffer);
		
		GenomicAnnotation right = new BasicGenomicAnnotation("right");
		right.setStart(getEnd() - buffer);
		right.setEnd(getEnd() + buffer);
		
		
		GenomicAnnotation A = new BasicGenomicAnnotation("A");
		A.setStart(twoSubjectAnnotation.getAStart());
		A.setEnd(twoSubjectAnnotation.getAEnd());
		
		GenomicAnnotation B = new BasicGenomicAnnotation("B");
		B.setStart(twoSubjectAnnotation.getBStart());
		B.setEnd(twoSubjectAnnotation.getBEnd());
		
		return (A.overlaps(left) && B.overlaps(right)) || (A.overlaps(right) && B.overlaps(left));
	}
	
	public List<GenomicAnnotation> minus(List<? extends GenomicAnnotation> others, boolean listIsDisjoint, int startAtIdx) {
		ArrayList<GenomicAnnotation> dif = new ArrayList<GenomicAnnotation>();
		dif.add(this);
		Collections.sort(others);

		Iterator<? extends GenomicAnnotation> it = null;
		if(!listIsDisjoint) {
			List<GenomicAnnotation> stitchedOthers = stitchList(others, 0);
			it = stitchedOthers.iterator();
		} else {
			it = others.listIterator(startAtIdx);
		}

		lastBreakInIndex = startAtIdx;
		while(it.hasNext()) {
			GenomicAnnotation a = it.next();
			if(a.getStart() > getEnd()) {
				break;
			}
			lastBreakInIndex++;
			ArrayList<GenomicAnnotation> newList = new ArrayList<GenomicAnnotation>(dif.size());
			Iterator<GenomicAnnotation> oldListIt = dif.iterator();
			while(oldListIt.hasNext()) {
				GenomicAnnotation oldRegion = oldListIt.next();
				newList.addAll(oldRegion.minus(a));
			}
			dif = newList;
		}
		
		//Collections.sort(dif, new PositionComparator());
		
		return dif;
	}
	
	public List<GenomicAnnotation> minus(List<? extends GenomicAnnotation> others) {
		return minus(others, false, 0);
	}
	
	public List<GenomicAnnotation> minus(GenomicAnnotation other) {
		List<GenomicAnnotation> dif = new ArrayList<GenomicAnnotation>();
		if(!overlaps(other)) {
			dif.add(this);
		} else {
			GenomicAnnotation intersection = new BasicGenomicAnnotation(this);
			intersection.takeIntersection(other);
			GenomicAnnotation first = new BasicGenomicAnnotation(this);
			first.setEnd(intersection.getStart() - 1);
			if(first.getLength() > 0) {
				dif.add(first);
			}
			GenomicAnnotation second = new BasicGenomicAnnotation(this);
			second.setStart(intersection.getEnd() + 1);
			if(second.getLength() > 0) {
				dif.add(second);
			}
		}
		return dif;
	}
	/**
	 * 
	 * @param sortedList
	 * @param maxDistanceToStitch
	 * @return
	 */	
	public static <T extends GenomicAnnotation> List<T> stitchList(List<? extends T> sortedList, int maxDistanceToStitch) {

		Stack<T> result = new Stack<T>();

		if(sortedList.size() == 0 ) {
			return result;
		}
		
		Iterator<? extends T> it = sortedList.iterator();
		result.push(it.next());
		while(it.hasNext()) {
			T next = it.next();
			T curr = result.pop();
			if(curr.overlaps(next,maxDistanceToStitch)) {
				curr.stitchTo(next);
				curr.setName(curr.getName()+"-"+next.getName());
				result.push(curr);
			} else {
				result.push(curr);
				result.push(next);
			}
		}
		
		return result;
	}
	
	/**
	 * 
	 * @param other GenomicAnnotation to check of overlap
	 * @return true if the current instance overlaps with the other one.
	 */
	public boolean overlaps(GenomicAnnotation other) {
		//System.out.print("Does this<"+getName()+"_"+getStart()+"-"+getEnd()+"> overlap <"+other.getName()+"_"+other.getStart()+"-"+other.getEnd()+">");
		//System.out.println(" !<"+(getStart() < other.getEnd() && getEnd() > other.getStart())+">!");
		return overlaps(other, 0);
	}
	
	/**
	 * 
	 * @param other - other genomic annotation
	 * @param buffer if the overlap is within buffer they will be considered overlapping 
	 *        even if they do not overlap within their original boundaries.
	 * @return true if they overlap in this extended definition
	 */
	public boolean overlaps(GenomicAnnotation other, int buffer) {

		return ((other.getChromosome() == null) || other.getChromosome().equals(getChromosome())) && (getStart() < other.getEnd() + buffer) && (getEnd() > other.getStart() - buffer);
	}
	
	public boolean contains(GenomicAnnotation other) {
		//System.out.println("Do this annotation " + toString() + " contains " + other );
		return ((other.getChromosome() == null) || other.getChromosome().equals(getChromosome())) && getStart() <= other.getStart() && getEnd() >= other.getEnd();
	}

	public String getName() {
		return name;
	}
	
	public void setName(String name) {
		this.name = name;
	}
	
	public int lastIdxWhereListSearchBroke() { return lastBreakInIndex;}
	
	public void setReversedOrientation(boolean isInReversedOrientation) {
		this.orientation = isInReversedOrientation ? "-" : "=";
	}

	public boolean inReversedOrientation() {
		return "-".equals(orientation);
	}
	public Sequence getSequence() {
		return sequence;
	}
	public void setSequence(Sequence sequence) {
		this.sequence = sequence;
	}
	/*
	public String getDisplayName() {
		return displayName;
	}
	public void setDisplayName(String displayName) {
		this.displayName = displayName;
	}
	*/
	public boolean isOrientation() {
		return "+".equals(orientation);
	}
	public void setOrientation(boolean orientation) {
		setOrientation(orientation ? "+" : "-");
	}
	public void setOrientation(String orientationString) {
		this.orientation = orientationString.intern();
	}
	public String getChromosome() {
		return chromosome;
	}
	public void setChromosome(String chromosome) {
		if(chromosome != null) {
			String chr = chromosome;
			if(chromosome.startsWith("chr")) {
				chr = chromosome.replace("chr", "");
			}
			this.chromosome = chr.intern();
		}
	}
	public void setThreePrimeBuffer(int bufferSize) {
		setEnd(getEnd() + bufferSize);
		this.threePrimerBuffer = bufferSize;
	}
	
	public void setFivePrimeBuffer(int bufferSize) {
		setStart(getStart() - bufferSize);
		this.fivePrimerBuffer = bufferSize;
	}
	
	public int getFivePrimeBases() {
		return fivePrimerBuffer;
	}


	public int getThreePrimeBases() {
		return threePrimerBuffer;
	}
	
	public String getOrientation() { return inReversedOrientation() ? "-" : "+";}
	


	/** 
	 * Default compareTo method: If annotations are in the same chromosome
	 * their start/end locations are compared otherwise their chromosomes are.
	 */
	public int compareTo(GenomicAnnotation arg0) {
		if(getChromosome() != null && arg0.getChromosome() != null && !getChromosome().equals(arg0.getChromosome())) {
			return getChromosome().compareTo(arg0.getChromosome()) * 1000000000;
		}
		return getStart() != arg0.getStart() ? getStart() - arg0.getStart() : getEnd() - arg0.getEnd();
	}
	public int getOrientedEnd() {
		return inReversedOrientation() ? getStart() : getEnd();
	}
	public int getOrientedStart() {
		// TODO Auto-generated method stub
		return inReversedOrientation() ? getEnd() : getStart();
	}
	public void addBlock(String name, int start, int end) {
		//Do nothing if does not support blocks, override if desired.
	}
	public List<? extends GenomicAnnotation> getBlocks() {
		return null;
	}
	public boolean mayHaveBlocks() {
		return false;
	}

}
