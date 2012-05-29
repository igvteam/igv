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

package edu.mit.broad.prodinfo.multiplealignment;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Stack;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import Jama.Matrix;
import edu.mit.broad.prodinfo.chromosome.BasicGenomicAnnotation;
import edu.mit.broad.prodinfo.genomicplot.GenomicAnnotation;
import edu.mit.broad.prodinfo.sequence.Sequence;
import edu.mit.broad.prodinfo.sequence.SequenceRegion;

public class MultipleAlignment  {
	private float score;
	private LinkedHashMap<String, AlignedSequence> alignments;

	private MultipleAlignmentIO ioHelper;
	public static final short A_CODE = Sequence.SHORT_ENCODED_A;
	public static final short C_CODE = Sequence.SHORT_ENCODED_C;
	public static final short G_CODE = Sequence.SHORT_ENCODED_G;
	public static final short T_CODE = Sequence.SHORT_ENCODED_T;
	public static final short GAP_CODE = Sequence.SHORT_ENCODED_GAP;

	protected static int UNGAPPED_ALPHABET_SIZE = 4;

	private String referenceId;
	private boolean isRefGapped;

	public MultipleAlignment() {
		alignments = new LinkedHashMap<String, AlignedSequence>();
	}

	/**
	 * Constructs a multiple alignment object out of piecemeal multiple alignment blocks with, possible,
	 * different number of aligned sequences.
	 * @param exonAlnBlocks The blocks to combined into a single multiple alignment block
	 * @param string reference sequence id
	 * @param seqs array of sequence ids, the end alignment will include all even in the cases where a sequence alignment is gap only.
	 */
	public MultipleAlignment(List<? extends MultipleAlignment> alnBlocks, String ref, String[] seqs) {
		//TODO: Implement.
	}

	public void encode() {
		Iterator<AlignedSequence> seqIt = alignments.values().iterator();
		while(seqIt.hasNext()) {
			AlignedSequence as = seqIt.next();
			as.encodeSequenceIgnoreCase();
			as.unloadSequence();
		}

	}

	public void encodeAsMatrix() {
		Iterator<AlignedSequence> seqIt = alignments.values().iterator();
		while(seqIt.hasNext()) {
			AlignedSequence as = seqIt.next();
			as.encodeSequenceAsVector();
			as.unloadSequence();
		}
	}

	public void reverse() {
		Iterator<AlignedSequence> seqIt = alignments.values().iterator();
		while(seqIt.hasNext()) {
			AlignedSequence as = seqIt.next();
			as.reverse();
		}

	}
	public int length() {
		return alignments.values().iterator().next().getLength();
	}

	/**
	 * Assumes that encode has been called and each aligned sequence is encoded
	 * @param int - alignment column
	 */
	public Map<String, Short> getColumn(int i) {
		LinkedHashMap<String, Short> col = new LinkedHashMap<String, Short>(getAlignedSequenceIds().size());
		Iterator<String> seqIdIt = getAlignedSequenceIds().iterator();
		while(seqIdIt.hasNext()) {
			String seqId = seqIdIt.next();
			col.put(seqId, alignments.get(seqId).getEncodedSequence()[i - getReferenceStart()]);
		}
		return col;
	}

	public Map<String, Matrix> getColumnsAsVector(int start, int number) throws ArrayIndexOutOfBoundsException{
		LinkedHashMap<String, Matrix> cols = new LinkedHashMap<String, Matrix>(getAlignedSequenceIds().size());
		Iterator<String> seqIdIt = getAlignedSequenceIds().iterator();
		start = start - getReferenceStart();
		//System.out.println("start: " + start + " ref-start " + getReferenceStart());
		while(seqIdIt.hasNext()) {
			String seqId = seqIdIt.next();
			short [] alignedSeqSeq = alignments.get(seqId).getEncodedSequence();
			Matrix alignmentMatrix = alignments.get(seqId).getVectorEncodedSequence();

			Matrix seqRegion = new Matrix(UNGAPPED_ALPHABET_SIZE, number);
			for(int j = start; j < start + number; j++) {
				if(alignedSeqSeq != null ) {
					if(alignedSeqSeq[j] < UNGAPPED_ALPHABET_SIZE) { // gapped sequence have a 0 column.
						seqRegion.set(alignedSeqSeq[j], j - start, 1);
					}
				} else if(alignmentMatrix != null) {
					for(int i = 0; i < alignmentMatrix.getRowDimension(); i++) {
						seqRegion.set(i, j - start , alignmentMatrix.get(i,j));
					}
				}
			}
			cols.put(seqId, seqRegion);
		}
		return cols;
	}

	public Map<String, Matrix> getAsMatrix() {
		LinkedHashMap<String, Matrix> matrixAlignment = new  LinkedHashMap<String, Matrix>(getAlignedSequenceIds().size());
		Iterator<String> seqIdIt = getAlignedSequenceIds().iterator();
		while(seqIdIt.hasNext()) {
			String seqId = seqIdIt.next();
			matrixAlignment.put(seqId, getAlignedSequence(seqId).getVectorEncodedSequence());
		}
		return matrixAlignment;
	}

	public void addShortEncodedColumn(Map<String, Short> col) {
		Iterator<String> it = col.keySet().iterator();
		while(it.hasNext()) {
			String seqId = it.next();
			AlignedSequence seq = getAlignedSequence(seqId);
			if(seq == null) {
				seq = new AlignedSequence(seqId);
				seq.setId(seqId);
				addSequence(seq);
			}
			seq.appendToSequence(decodeShort(col.get(seqId)));
		}

	}

	public void addShortEncodedRegion(Map<String, short[]> region) {
		if(region == null || region.isEmpty()) {
			return;
		}

		int cols = region.values().iterator().next().length;
		for(int i = 0; i < cols; i++) {
			Iterator<String> it = region.keySet().iterator();
			while(it.hasNext()) {
				String seqId = it.next();
				AlignedSequence seq = getAlignedSequence(seqId);
				if(seq == null) {
					seq = new AlignedSequence(seqId);
					seq.setId(seqId);
					addSequence(seq);
				}
				seq.appendToSequence(decodeShort(region.get(seqId)[i]));
			}
		}

	}


	public char decodeShort(Short nucleotide) {
		char decoded = ' ';
		switch (nucleotide) {
		case A_CODE : decoded = 'A'; break;
		case C_CODE : decoded = 'C'; break;
		case G_CODE : decoded = 'G'; break;
		case T_CODE : decoded = 'T'; break;
		default : decoded = '-';
		}
		return decoded;
	}

	public void addSequence(AlignedSequence seq) {
		if(seq == null) {
			throw new IllegalArgumentException("Trying to add a null sequence");
		}
		if(alignments.keySet().contains(seq.getContainingSequenceId())) {
			throw new IllegalArgumentException("Sequence " + seq.getContainingSequenceId() + " has already been added");
		}

		alignments.put(seq.getContainingSequenceId(), seq);
	}

	public void addSequences(List<AlignedSequence> sequences) {
		Iterator<AlignedSequence> it = sequences.iterator();
		while(it.hasNext()) {
			addSequence(it.next());
		}
	}

	public boolean isEmpty() { return alignments.isEmpty();}

	public void setIOHelper(MultipleAlignmentIO maio) { this.ioHelper = maio;}
	protected MultipleAlignmentIO getIOHelper() { return ioHelper;}

	public void write(BufferedWriter bw) throws IOException {
		if(ioHelper == null) {
			throw new IllegalStateException("The MultipleAlignmentIO (ioHelper) is null set it first");
		}
		ioHelper.write(bw, this);
	}

	public void write(BufferedWriter bw, List<String> order) throws IOException {
		if(ioHelper == null) {
			throw new IllegalStateException("The MultipleAlignmentIO (ioHelper) is null set it first");
		}
		ioHelper.write(bw, this, order);
	}

	public MultipleAlignment sampleColumns(int colNum, int numConsecutiveCols) {
		Random r = new Random();
		LinkedHashMap<String, List<Character>> sequenceCharacters = new LinkedHashMap<String, List<Character>>(getAlignedSequenceIds().size());

		Iterator<String> seqIdIt = getAlignedSequenceIds().iterator();

		MultipleAlignment sampled = new MultipleAlignment();

		while(seqIdIt.hasNext()) {
			String seq = seqIdIt.next();
			AlignedSequence alignedSeq = getAlignedSequence(seq);
			char [] seqBasesArray = alignedSeq.getSequenceBases().toCharArray();
			List<Character> seqChars = new ArrayList<Character>(seqBasesArray.length);

			sequenceCharacters.put(seq, seqChars);
			for(int i = 0; i < seqBasesArray.length; i++) {
				seqChars.add(seqBasesArray[i]);
			}
			//alignedSeq.unloadSequence();
			AlignedSequence newSeq = new AlignedSequence(alignedSeq);
			newSeq.unloadSequence();
			sampled.addSequence(newSeq);
		}

		int length = sequenceCharacters.values().iterator().next().size();
		sampled.setReferenceId(getReferenceId());

		for(int i = 0; i < colNum; i++) {
			int col = r.nextInt(length - numConsecutiveCols + 1);
			//System.out.println("Sampled col " + col + " out of " + length + " cols");
			seqIdIt = getAlignedSequenceIds().iterator();
			while(seqIdIt.hasNext()) {
				String seqId = seqIdIt.next();
				AlignedSequence seq = sampled.getAlignedSequence(seqId);
				List<Character> seqChars = sequenceCharacters.get(seqId);
				for(int j = 0; j < numConsecutiveCols; j++) {
					seq.appendToSequence(seqChars.get(col + j));
				}
			}
		}
		//System.out.println(getReferenceId() + " bases: " + sampled.getAlignedSequence(getReferenceId()).getSequenceBases());
		return sampled;
	}

	public void permuteColumns() {
		Random r = new Random();
		LinkedHashMap<String, List<Character>> sequenceCharacters = new LinkedHashMap<String, List<Character>>(getAlignedSequenceIds().size());

		Iterator<String> seqIdIt = getAlignedSequenceIds().iterator();

		while(seqIdIt.hasNext()) {
			String seq = seqIdIt.next();
			AlignedSequence alignedSeq = getAlignedSequence(seq);
			char [] seqBasesArray = alignedSeq.getSequenceBases().toCharArray();
			List<Character> seqChars = new ArrayList<Character>(seqBasesArray.length);

			sequenceCharacters.put(seq, seqChars);
			for(int i = 0; i < seqBasesArray.length; i++) {
				seqChars.add(seqBasesArray[i]);
			}
			alignedSeq.unloadSequence();
		}

		int length = sequenceCharacters.values().iterator().next().size();
		while (length > 0) {
			int pos = r.nextInt(length);
			seqIdIt = alignments.keySet().iterator();

			while(seqIdIt.hasNext()) {
				String seqId = seqIdIt.next();
				AlignedSequence seq = alignments.get(seqId);
				List<Character> seqChars = sequenceCharacters.get(seqId);
				seq.appendToSequence(seqChars.remove(pos));
			}
			length--;
		}
	}

	public float getScore() { return score;}
	public void setScore(float score) { this.score = score;}

	public AlignedSequence getAlignedSequence(String containingSequenceId) { return alignments.get(containingSequenceId); }

	public List<AlignedSequence> getAlignedSequences() { return new ArrayList<AlignedSequence>(alignments.values());}

	public List<String> getAlignedSequenceIds() { return new ArrayList<String>(alignments.keySet());}

	public boolean overlaps(SequenceRegion region) {
		boolean overlap = false;
		AlignedSequence seq = alignments.get(region.getContainingSequenceId());
		if (seq != null) {
			overlap = seq.overlaps(region);
		}
		return overlap;
	}

	public String getReferenceId() {
		String refId = null;
		if(referenceId != null) {
			refId = referenceId;
		} else if (!alignments.isEmpty()){
			refId =  alignments.keySet().iterator().next() ;
		}
		return refId;
	}

	public void setReferenceId(String referenceId) {
		this.referenceId = referenceId.intern();
	}

	public int getReferenceStart() {
		AlignedSequence refChunk = getReference();
		return refChunk == null ? -1 : refChunk.getStart();
	}

	public AlignedSequence getReference() {
		AlignedSequence reference = null;
		if(referenceId != null) {
			reference = getAlignedSequence(referenceId);
		} else if(alignments.size() > 0) {
			reference = alignments.values().iterator().next();
		}

		return reference;
	}

	public int getReferenceEnd() {
		AlignedSequence refChunk = getReference();
		return refChunk == null ? -1 : refChunk.getEnd();
	}

	public boolean isRefGapped() {
		return isRefGapped;
	}

	public void setRefGapped(boolean isRefGapped) {
		this.isRefGapped = isRefGapped;
	}

	protected Iterator<AlignedSequence> getAlignedSequences(String refId, int refStart, int refEnd, boolean reverse) {
		ArrayList<AlignedSequence> seqs = new ArrayList<AlignedSequence>(alignments.size());

		return seqs.iterator();
	}

	public MultipleAlignment getConcatenatedSubAlignment(String refId, List<? extends GenomicAnnotation> annotations) {
		MultipleAlignment result = new MultipleAlignment();
		result.setIOHelper(this.ioHelper);
		Iterator<? extends GenomicAnnotation> it = annotations.iterator();
		GenomicAnnotation ga = null;
		while(it.hasNext()) {
			ga = it.next();
			result.append(getSubAlignment(refId, ga.getStart(), ga.getEnd(), ga.inReversedOrientation()));
		}
		return result;
	}

	/**
	 * Compresses alignment by removing columns containing reference gaps
	 */
	public void compress() {
		HashMap<String, StringBuilder> compressSequenceMap = new HashMap<String, StringBuilder>(getAlignedSequenceIds().size());
		List<int []> ungappedRefIslands = getUngappedReferenceIslands();
		Iterator<int []> refIslandIt = ungappedRefIslands.iterator();
		while(refIslandIt.hasNext()) {
			Iterator<AlignedSequence> seqIt = getAlignedSequences().iterator();
			int [] ungappedIsland = refIslandIt.next();
			while(seqIt.hasNext()) {
				AlignedSequence seq = seqIt.next();
				StringBuilder compressedSeqBuilder = compressSequenceMap.get(seq.getId());
				if(compressedSeqBuilder == null) {
					compressedSeqBuilder = new StringBuilder();
					compressSequenceMap.put(seq.getId(), compressedSeqBuilder);
				}
				compressedSeqBuilder.append(seq.getSequenceBases().substring(ungappedIsland[0], ungappedIsland[1]));
			}
		}

		Iterator<String> seqIdIt = compressSequenceMap.keySet().iterator();
		while(seqIdIt.hasNext()) {
			String seqId = seqIdIt.next();
			AlignedSequence seq = getAlignedSequence(seqId);
			seq.setSequenceBases(compressSequenceMap.get(seqId).toString());
		}
		isRefGapped = false;
	}

	public void remove(List<String> removeList) {
		Iterator<String> it = removeList.iterator();

		while(it.hasNext()) {
			String toRemove = it.next();
			System.out.print("To remove " + toRemove + " ... ");
			AlignedSequence removedSeq = alignments.remove(toRemove);
			System.out.println(removedSeq == null ? " not found " : "yes");
		}
	}

	public MultipleAlignment getSubAlignment(int refStart, int refEnd, boolean reverse) {
		return getSubAlignment(getReferenceId(), refStart, refEnd, reverse);
	}

	/**
	 * Computes the number of sequences that have aligning bases as opposed to being pure gap
	 * @return Number of sequences that are not gap only.
	 */
	public int getNumberOfAligningSequences() {
		int i = 0;
		Iterator<AlignedSequence> seqIt = getAlignedSequences().iterator();
		while(seqIt.hasNext()) {
			AlignedSequence seq = seqIt.next();
			if (seq.getUngappedLength() > 0 ) {
				i++;
			}
		}

		return i;
	}

	/**
	 *
	 * @param refId - The sequence id in which coordinates are given
	 * @param refStart - Start of region of interest in the refId sequence. Positions start in 1.
	 * @param refEnd - End of region of interest in the refId sequence. Positions start in 1.
	 * @param reverse - If true the multiple alignment will return all sequences in the opposite strand
	 * @return The multiple alignment of the region of interest or an empty alignment if the coordinates
	 * 			do not overlap the refId aligned region.
	 */
	public MultipleAlignment getSubAlignment(String refId, int refStart, int refEnd, boolean reverse) {
		MultipleAlignment result = new MultipleAlignment();
		BasicGenomicAnnotation target = new BasicGenomicAnnotation("ReferenceTarget");
		target.setStart(refStart);
		target.setEnd(refEnd);

		AlignedSequence ref = alignments.get(refId) ;
		if(ref == null) {
			throw new IllegalArgumentException("Bad refId, multiple alignment does not include " + refId);
		}
		target.takeIntersection(ref);
		//System.out.println("RefStart " + refStart + " ref.getRegionStart() " + ref.getRegionStart() + " refEnd " + refEnd );
		int stringStartPos = ref.getGapAdjustedCoordinate(target.getStart() - ref.getRegionStart());
		int stringEndPos   = ref.getGapAdjustedCoordinate(target.getEnd()   - ref.getRegionStart());
		//System.out.println("gap adjusted start " + stringStartPos + " gap adjusted end " +stringEndPos );

		Iterator<AlignedSequence> it = getAlignedSequences().iterator();
		SequenceRegion tmpRegion = null;
		while(it.hasNext()) {
			AlignedSequence seq = it.next();
			//System.out.print(seq.getId());
			tmpRegion = seq.getRegion(stringStartPos, stringEndPos);
			//System.out.println(" " + tmpRegion.getSequenceBases());
			if(reverse){
				tmpRegion.reverse();
			}
			AlignedSequence newSeq = new AlignedSequence(tmpRegion.getContainingSequenceId());
			newSeq.setRegionStart(seq.getRegionStart() + seq.getSequencePosition(stringStartPos));
			newSeq.setRegionEnd(seq.getRegionStart() + seq.getSequencePosition(stringEndPos));
			newSeq.setChromosome(seq.getChromosome());
			newSeq.setSequence(tmpRegion);
			newSeq.setName(tmpRegion.getContainingSequenceId());
			result.addSequence(newSeq);
		}

		result.setIOHelper(this.ioHelper);
		return result;

	}
	public void append(MultipleAlignment other) {
		if(alignments.isEmpty()) {
			Iterator<AlignedSequence> otherSeqIt =  other.getAlignedSequences().iterator();
			//System.out.println("ailignments are empty, setting to others");
			while(otherSeqIt.hasNext()) {
				AlignedSequence seq = otherSeqIt.next();
				AlignedSequence seqClone = new AlignedSequence(seq);
				alignments.put(seq.getId(), seqClone);
			}
		} else {
			//Iterator<AlignedSequence> it = other.getAlignedSequences().iterator();
			Iterator<AlignedSequence> it = getAlignedSequences().iterator();
			AlignedSequence otherSeq = null;
			AlignedSequence thisSeq = null;
			while(it.hasNext()) {
				thisSeq = it.next();
				//System.out.println("This seq " + thisSeq.getId() + " bases " + thisSeq.getSequenceBases());
				otherSeq  = other.getAlignedSequence(thisSeq.getName());
				//System.out.println("Other seq " + otherSeq.getId() + " bases " + otherSeq.getSequenceBases());
				if (otherSeq == null) {
					throw new IllegalStateException("Sequence " + thisSeq.getName() +" has no entry in other alignment, can't append alignment");
				}
				thisSeq.appendToSequence(otherSeq.getSequenceBases());
				thisSeq.setEnd(thisSeq.getEnd() + otherSeq.getUngappedLength());
			}
		}
	}

	public static class AlignedSequence extends SequenceRegion {
		private static final Pattern GAP_PATTERN = Pattern.compile("-+");
		private static final Pattern UNGAP_PATTERN = Pattern.compile("[^-]+");
		private int totalLength;
		Stack<int[]> ungappedRegions;

		public AlignedSequence(Sequence sequence) {
			super(sequence.getId());
			setId(sequence.getId());
			setName(sequence.getId());
			setSequence(sequence);
			setRegionStart(0);
			setRegionEnd(getLength());
			setSequenceBases(getSequenceBases());

		}

		public void setSequenceBases(String bases ) {
			super.setSequenceBases(bases);
		}

		public AlignedSequence(String containingSequenceId) {
			super(containingSequenceId);
			setId(containingSequenceId);
		}

		public int getLength() {
			if(getEncodedSequence() != null && getEncodedSequence().length > 0) {
				return getEncodedSequence().length;
			} else if (getVectorEncodedSequence() != null) {
				return getVectorEncodedSequence().getColumnDimension();
			} else if (getSequenceBases() != null && getSequenceBases().length() > 0) {
				return getSequenceBases().length();
			}else {
				return super.getLength();
			}
		}

		public int getUngappedLength() {
			int ungappedLength = 0;
			char [] basesArr = getSequenceBases().toCharArray();
			int gaps = 0;
			for(int i = 0; i < basesArr.length; i++) {
				if(basesArr[i] == '-') {
					gaps++;
				}

			}
			ungappedLength = getLength() - gaps;
			return ungappedLength;
		}


		public void setStrand(String strand) { super.setForwardStrand("-".equals(strand)); }

		public float getPercentGaps() {
			return 1f - getUngappedLength()/(float)getLength();
		}

		public int getTotalGaps() { return getLength() - getUngappedLength();}

		public List<Integer> getGapSizes() {
			Matcher m = GAP_PATTERN.matcher(getSequenceBases());
			ArrayList<Integer> gapSizes = new ArrayList<Integer>();
			while(m.find()) {
				m.start();
				m.end();
				gapSizes.add(m.group().length());
			}
			return gapSizes;
		}

		/**
		 * Finds the list of ungapped regions in the sequence.
		 * @return A list of two sized integer arrays with the start and end of the ungapped regions. The convention here is semi closed
		 * 		   intervals, each list items is of the form [start, end).
		 */
		public List<int []> findUngappedChunks() {
			Matcher m = UNGAP_PATTERN.matcher(getSequenceBases());
			ArrayList<int []> gapUngapp = new ArrayList<int []>();
			//System.out.println("Sequence for " + getId() + ": " + getSequenceBases());
			while(m.find()) {
				int [] ungappedReg = {m.start(), m.end()};
				//System.out.println("\treg["+m.start()+"-"+m.end()+"]: " + getSequenceBases().substring(m.start(), m.end()));
				gapUngapp.add(ungappedReg);
			}
			//System.out.println("");

			return gapUngapp;
		}

		/**
		 * Maps coordinates relative to start of alignment (0 based) to the gapped adjusted string position in the aligned sequence string.
		 * @param position to map. It is the inverse of getSequencePostion.
		 * @return
		 * @see #getSequencePosition(int)
		 */
		public int getGapAdjustedCoordinate(int position) {
			List<int []> ungappedChunks = findUngappedChunks();

			Iterator<int []> chunkIt = ungappedChunks.iterator();
			int gapNum = 0;
			int [] lastChunk = {0,0};
			while(chunkIt.hasNext()) {
				int [] chunk = chunkIt.next();
				int gapsBetweenChunks = chunk[0] - lastChunk[1];
				if(chunk[0] - gapNum - gapsBetweenChunks <=  position) {
					gapNum += gapsBetweenChunks;
				} else {
					break;
				}
				lastChunk = chunk;
			}

			return position + gapNum;
		}

		/**
		 * Maps a string position to a sequence position relative to the start of the alignment. It returns the position once gaps are
		 * ignored. It is the inverse of getGapAdjustedCoordinate
		 * @param coordinate
		 * @return position in sequence.
		 * @see #getGapAdjustedCoordinate(int);
		 */
		public int getSequencePosition(int coordinate) {
			List<int []> ungappedChunks = findUngappedChunks();

			Iterator<int []> chunkIt = ungappedChunks.iterator();
			int gapNum = 0;
			int [] lastChunk = {0,0};
			while(lastChunk[1] < coordinate && chunkIt.hasNext()) {
				int [] chunk = chunkIt.next();
				int gapsBetweenChunks = chunk[0] - Math.min(lastChunk[1], coordinate);
				if(chunk[0]  <=  coordinate) {
					gapNum += gapsBetweenChunks;
				}
				lastChunk = chunk;
			}

			return coordinate - gapNum;
		}

		public int getTotalLength() {
			return totalLength;
		}

		public void setTotalLength(int totalLength) {
			this.totalLength = totalLength;
		}


	}

	/**
	 * Finds the list of ungapped regions in the reference sequence.
	 * @return A list of two sized integer arrays with the start and end of the ungapped regions. The convention here is semi closed
	 * 		   intervals, each list items is of the form [start, end).
	 */
	public List<int[]> getUngappedReferenceIslands() {
		//System.out.println("Reference: " + getReferenceId() + " reference sequence " + getReference());
		return getReference().findUngappedChunks();
	}

	/**
	 * Finds the list of ungapped sequence (as opposed to alignment gaps) regions in the reference sequence.
	 * @return A list of two sized integer arrays with the start and end of the ungapped regions. The convention here is semi closed
	 * 		   intervals, each list items is of the form [start, end).
	 */
	public List<int[]> getUngappedSequenceReferenceIslands() {
		//System.out.println("Reference: " + getReferenceId() + " reference sequence " + getReference());
		return getReference().findUngappedSequenceChunks();
	}

	public void introduceGaps(boolean canGapReference, int maxNumberOfGaps) throws IllegalAccessException {
		Random r = new Random();

		int alnLength = length();
		List<String> seqIds = new ArrayList<String>(getAlignedSequenceIds());
		if(!canGapReference) {
			seqIds.remove(getReferenceId());
		}

		for (int i = 0; i < alnLength; i++) {
			List<String> availableSeqIds = new ArrayList<String>(seqIds);

			int numOfGaps = r.nextInt(maxNumberOfGaps);

			for(int j = 0; j < numOfGaps; j++) {
				int toGapIdx = r.nextInt(availableSeqIds.size());
				String seqToGap = availableSeqIds.remove(toGapIdx);
				Sequence seq = alignments.get(seqToGap);
				seq.setCharAt(i, '-');
			}
		}
		
	}





}
