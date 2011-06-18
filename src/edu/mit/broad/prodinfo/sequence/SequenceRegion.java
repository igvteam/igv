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

package edu.mit.broad.prodinfo.sequence;

import java.util.List;

import edu.mit.broad.prodinfo.chromosome.BasicGenomicAnnotation;
import edu.mit.broad.prodinfo.genomicplot.GenomicAnnotation;
import edu.mit.broad.prodinfo.genomicplot.TwoSubjectAnnotation;

public class SequenceRegion extends Sequence implements GenomicAnnotation {
	BasicGenomicAnnotation annotation;
	private String containingSequenceId;
	public static final int INF = 1000000000;
	
	public SequenceRegion(String containingSequenceId) {
		super(null);
		this.containingSequenceId = containingSequenceId;
		annotation = new BasicGenomicAnnotation(containingSequenceId);
		annotation.setSequence(this);
		annotation.setEnd(INF);
	}
	public int getRegionEnd() {
		return annotation.getEnd() == INF && super.getSequenceBases() != null && getSequenceBases().length() > 0 
			? getSequenceBases().length() + annotation.getStart() - 1
			: annotation.getEnd();
	}
	public void setRegionEnd(int regionEnd) {
		annotation.setEnd(regionEnd);
	}
	public int getRegionStart() {
		return annotation.getStart();
	}
	public void setRegionStart(int regionStart) {
		annotation.setStart(regionStart);
	}
	
	public WindowSlider getSlider(int windowSize, int overlap) {
		return WindowSlider.getSlider(this, windowSize, overlap);
	}
	
	public int getLength() {
		if(getEnd() == INF && getSequenceBases() == null && getSequenceBases().length() == 0) {
			return INF;
		} else if(getEnd() < INF) {
			return getEnd() - getStart();
		} else {
			return super.getLength();
		}
	}
	
	public String getOrientation() { return inReversedOrientation() ? "-" : "+"; }
	public void setReversedOrientation(boolean reversed) { annotation.setOrientation(!reversed); }
	
	public String getContainingSequenceId() {
		return containingSequenceId;
	}
	
	public String getId() {
		return super.getId() == null ? containingSequenceId + "_" + getStart() + "-" + getEnd() : super.getId();
	}
	
	public int getStart() {
		return annotation.getStart();
	}
	
	public int getEnd() {
		return annotation.getEnd();
	}
	
	public double getScore() {
		return annotation.getScore();
	}
	
	public void setScore(double score) {
		annotation.setScore(score);
	}
	
	
	public String getName() {
		return getId();
	}
	
	public void setId(String id) {
		super.setId(id);
		annotation.setName(id);
	}
	
	public void setName(String name) {
		setId(name);
	}
	
	public boolean inReversedOrientation() {
		return annotation.inReversedOrientation();
	}
	
	public long getMiddle() {
		return annotation.getMiddle();
	}
	
	public void setStart(int start) {
		setRegionStart(start);
	}
	
	public void setEnd(int end) {
		setRegionEnd(end);
	}
	public Sequence getSequence() {
		return this;
	}
	public void setSequence(Sequence seq) {
		super.setSequenceBases(seq.getSequenceBases());
	}
	
	public String getChromosome() {
		return annotation.getChromosome();
	}
	
	public void setChromosome(String chr) {
		annotation.setChromosome(chr);
	}
	public boolean overlaps(GenomicAnnotation other, int buffer) {
		return annotation.overlaps(other, buffer);
	}
	public boolean overlaps(GenomicAnnotation other) {
		return annotation.overlaps(other);
	}
	
	public List<GenomicAnnotation> minus(GenomicAnnotation other) {
		return annotation.minus(other);
	}
	
	public List<GenomicAnnotation> minus(List<? extends GenomicAnnotation> others) {
		return annotation.minus(others);
	}
	
	public void takeIntersection(GenomicAnnotation other) {
		annotation.takeIntersection(other);
	}
	
	public void takeUnion(GenomicAnnotation other) {
		annotation.takeUnion(other);
	}
	public void stitchTo(GenomicAnnotation other) {
		annotation.stitchTo(other);
	}
	
	public boolean isFlankedBy(TwoSubjectAnnotation twoSubjectAnnotation, int buffer) {
		return annotation.isFlankedBy(twoSubjectAnnotation, buffer);
	}
	public int getFivePrimeBases() {
		return annotation.getFivePrimeBases();
	}
	public int getThreePrimeBases() {
		return annotation.getThreePrimeBases();
	}
	
	public String toString() {
		return getContainingSequenceId() + ":" + getRegionStart() + "-" + getRegionEnd();
	}
	
	public String getLocationString() { return annotation.getLocationString();}
	public String getChromosomeString() { return annotation.getChromosomeString(); }
	
	public SequenceRegion extractRegionBasedOnGC(float targetGC, int size, int buffer) {
		SequenceRegion theRegion = super.extractRegionBasedOnGC(targetGC, size, buffer);
		if(theRegion != null) {
			theRegion.setRegionStart(getRegionStart() + theRegion.getRegionStart());
			theRegion.setRegionEnd(getRegionStart() + theRegion.getRegionEnd());
			theRegion.setChromosome(annotation.getChromosome());
		}
		return theRegion;
	}
	public boolean contains(GenomicAnnotation other) {
		return annotation.contains(other);
	}
	public int getDistanceTo(GenomicAnnotation other) {
		return annotation.getDistanceTo(other);
	}
	public void setOrientation(String orientation) {
		annotation.setOrientation(orientation);
		
	}
	public int compareTo(GenomicAnnotation arg0) {
		return annotation.compareTo(arg0);
	}
	public List<GenomicAnnotation> disect(GenomicAnnotation a) {
		return annotation.disect(a);
	}
	public List<GenomicAnnotation> disect(List<? extends GenomicAnnotation> disectors) {
		return annotation.disect(disectors);
	}
	
	public SequenceRegion getRegion(int start, int end) {
		SequenceRegion region = super.getRegion(start, end);
		if(annotation.getChromosome() != null) {
			region.setChromosome(annotation.getChromosome());
		}
		return region;
	}
	public int getOrientedEnd() {
		return annotation.getOrientedEnd();
	}
	public int getOrientedStart() {
		return annotation.getOrientedStart();
	}
	public void addBlock(String name, int start, int end) {
		//DO nothing as there is nothing to do here.
		
	}
	public List<? extends GenomicAnnotation> getBlocks() {
		return null;
	}
	public boolean mayHaveBlocks() {
		return false;
	}
	
}
