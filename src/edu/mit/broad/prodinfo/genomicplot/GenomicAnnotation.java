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

package edu.mit.broad.prodinfo.genomicplot;

import java.util.List;

import edu.mit.broad.prodinfo.sequence.Sequence;

public interface GenomicAnnotation extends Cloneable, Comparable<GenomicAnnotation> ,Feature {
	int getStart();
	int getEnd();
	void setScore(double score);
	void setName(String name);
	//String getDisplayName();
	boolean inReversedOrientation();
	void setReversedOrientation(boolean isInReversedOrientation);
	String getOrientation();
	long getMiddle();
	
	int getLength();
	
	void setStart(int start);
	void setEnd(int end);
	Sequence getSequence();
	void setSequence(Sequence seq);
	String getChromosome();
	void setChromosome(String chr);
	
	public void setOrientation(String orientationString);
	
	/**
	 * Calculates the distance to the another genomic annotation.
	 * @return 0 if the annotations overlap or the minimum between the edges otherwise.
	 */
	public int getDistanceTo(GenomicAnnotation other);
	
	/**
	 * Genomic annotation integer identification if such exists.
	 */
	public String getId();
	public void setId(String id);
	
	/**
	 * 
	 * @param other - other genomic annotation
	 * @param buffer if the overlap is within buffer they will be considered overlapping 
	 *        even if they do not overlap within their original boundaries.
	 * @return true if they overlap in this extended definition
	 */
	public boolean overlaps(GenomicAnnotation other, int buffer);
	/**
	 * 
	 * @param other GenomicAnnotation to check of overlap
	 * @return true if the current instance overlaps with the other one.
	 */
	public boolean overlaps(GenomicAnnotation other);
	
	/**
	 * Change the current instance to represent its intersection
	 * with the provided annotation 
	 * @param other
	 */
	public void takeIntersection(GenomicAnnotation other);
	
	/**
	 * Change the current instance to represent its union
	 * with the provided annotation 
	 * @param GenomicAnnotation other
	 */
	public void takeUnion(GenomicAnnotation other);
	
	/**
	 * Enlarges this annotation by 
	 * @param other
	 */
	public void stitchTo(GenomicAnnotation other);
	
	/**
	 * Checks whether a twoSubjectAnnotation flanks the instance
	 */
	boolean isFlankedBy(TwoSubjectAnnotation twoSubjectAnnotation, int buffer);
	
	
	/**
	 * @return int - the five prime buffer (i.e. bases not belonging to the annotation proper).
	 */
	public int getFivePrimeBases();

	/**
	 * @return int - the three prime buffer (i.e. bases not belonging to the annotation proper).
	 */
	public int getThreePrimeBases();
	
	public boolean contains(GenomicAnnotation other);
	
	/**
	 * Returns the difference (all regions in this genomic annotation not in the given list)
	 * between this genomic annotation and the given list
	 * 
	 * @param others - the annotations to take the difference against
	 * @return A List of genomic annotations of the resulting difference.
	 */
	public List<GenomicAnnotation> minus(List<? extends GenomicAnnotation> others); 
	
	/**
	 * Returns the difference (all regions in this genomic annotation not in the given other annotation)
	 * between this genomic annotation and the one
	 * 
	 * @param other - the annotations to take the difference against
	 * @return A List of genomic annotations of the resulting difference.
	 */
	public List<GenomicAnnotation> minus(GenomicAnnotation other);
	
	/**
	 * Fragments the current annotation if it overlaps the provided one.
	 * it returns an empty list if this annotation does not ovelap the
	 * one passed in.
	 * 
	 * @param GenomicAnnotation annotation to disect the current one
	 * @return List<GenomicAnnotation> of the disected annotations, an empty list is returned if 
	 * 			the annotations do not overlap.
	 */
	List<GenomicAnnotation> disect(GenomicAnnotation a);
	
	
	/**
	 * Fragments the current annotation if it overlaps the provided ones.
	 * It returns a list with one component (this annotation) if no annotation 
	 * in the provided list overlaps the discted annotaion.
	 * 
	 * @param List<GenomicAnnotation> <b>sorted</b> annotations with which to disect the current one.
	 * @return List<GenomicAnnotation> of the disected annotations, a list just this annotation is returned if 
	 * 			the annotations do not overlap.
	 */
	List<GenomicAnnotation> disect(List<? extends GenomicAnnotation> disectors);
	
	
	/**
	 * Gets the start of the annotation considering its orientation 
	 * @return getStart() if the feature is in direct orientation or getEnd() otherwise
	 */
	int getOrientedStart();
	
	/**
	 * Gets the end of the annotation considering its orientation 
	 * @return getEnd() if the feature is in direct orientation or getStart() otherwise
	 */
	int getOrientedEnd();
	
	/**
	 * Returns a string of the form chr<chromosome>:start-end
	 * @return
	 */
	public String getLocationString();

	/**
	 * Returns the ussual way to write a chromosome: chrSymbol if a chromosome or the full name of the scaffold 
	 * 
	 */
	public String getChromosomeString();
	
	/**
	 * Returns true if the annotation (like a full BED or a RefSeq) has sub annotations like exons or blocks
	 */
	public boolean mayHaveBlocks();
	
	/**
	 * Returns a list of blocks if the annotations has any (@see containsBlocks)
	 */
	public List<? extends GenomicAnnotation> getBlocks();
	
	/**
	 * Adds a block to the annotation if it supports blocks
	 */
	public void addBlock(String name, int start, int end);
}
