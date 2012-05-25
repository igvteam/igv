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

import Jama.Matrix;
import edu.mit.broad.prodinfo.chromosome.BasicGenomicAnnotation;
import edu.mit.broad.prodinfo.genomicplot.GenomicAnnotation;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Random;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Sequence {
    private static final Pattern SOFT_MASKED_PAT = Pattern.compile("[acgt]+");
    private static final Pattern SEQUENCE_GAP_PAT = Pattern.compile("[Nn]+");
    private static final Pattern UNGAPPED_SEQUENCE_PAT = Pattern.compile("[^Nn]+");

    private String id;
    private StringBuilder sequenceBases;
    private short[] encodedSequence;
    private Matrix vectorEncodedSequence;
    private boolean forwardStrand = true;
    private boolean encodeIgnoreCase = false;

    public static final char[] SHORT_READ_FLOW = {'T', 'A', 'C', 'G'};
    public static final char[] LONG_READ_FLOW = {'G', 'A', 'T', 'C'};

    public static final short SHORT_ENCODED_A = 0;
    public static final short SHORT_ENCODED_C = 1;
    public static final short SHORT_ENCODED_G = 2;
    public static final short SHORT_ENCODED_T = 3;
    public static final short SHORT_ENCODED_GAP = 4;
    public static final short SHORT_ENCODED_a = 5;
    public static final short SHORT_ENCODED_c = 6;
    public static final short SHORT_ENCODED_g = 7;
    public static final short SHORT_ENCODED_t = 8;
    public static final short SHORT_ENCODED_N = 9;

    public Sequence(String id) {
        super();
        this.id = id;
        sequenceBases = new StringBuilder();
    }

    public Sequence(String id, boolean isLarge) {
        this(id);
        if (isLarge) {
            sequenceBases = new StringBuilder(850000000);
        }
    }

    public String getId() {
        return id;
    }

    protected void setForwardStrand(boolean isForwardStrand) {
        this.forwardStrand = isForwardStrand;
    }

    public void setId(String id) {
        this.id = id;
    }

    public void unloadSequence() {
        sequenceBases.delete(0, sequenceBases.length());
        sequenceBases.trimToSize();
    }

    public int getLength() {
        return sequenceBases.length();
    }

    public void setCharAt(int position, char newCharacter) throws IllegalAccessException {
        if (sequenceBases == null) {
            throw new IllegalAccessException("This methods is only implemented for un encoded sequences. sequenceBases is null, it must be non null.");
        }

        sequenceBases.setCharAt(position, newCharacter);
    }

    public String getSequenceBases() {
        String bases = "";
        if (sequenceBases.length() > 0) {
            bases = sequenceBases.toString();
        } else if (encodedSequence != null && encodedSequence.length > 0) {
            StringBuffer buf = new StringBuffer(encodedSequence.length);
            if (!encodeIgnoreCase) {
                for (int i = 0; i < encodedSequence.length; i++) {
                    switch (encodedSequence[i]) {
                        case SHORT_ENCODED_c:
                            buf.append("c");
                            break;
                        case SHORT_ENCODED_C:
                            buf.append("C");
                            break;
                        case SHORT_ENCODED_G:
                            buf.append("G");
                            break;
                        case SHORT_ENCODED_g:
                            buf.append("g");
                            break;
                        case SHORT_ENCODED_a:
                            buf.append("a");
                            break;
                        case SHORT_ENCODED_A:
                            buf.append("A");
                            break;
                        case SHORT_ENCODED_T:
                            buf.append("T");
                            break;
                        case SHORT_ENCODED_t:
                            buf.append("t");
                            break;
                        case SHORT_ENCODED_GAP:
                            buf.append("-");
                            break;
                        default:
                            buf.append("N");
                            break;
                    }
                }

            } else {
                for (int i = 0; i < encodedSequence.length; i++) {
                    switch (encodedSequence[i]) {
                        case SHORT_ENCODED_A:
                            buf.append("A");
                            break;
                        case SHORT_ENCODED_C:
                            buf.append("C");
                            break;
                        case SHORT_ENCODED_G:
                            buf.append("G");
                            break;
                        case SHORT_ENCODED_T:
                            buf.append("T");
                            break;
                        case SHORT_ENCODED_GAP:
                            buf.append("-");
                            break;
                        default:
                            buf.append("N");
                            break;
                    }
                }
            }
            bases = buf.toString();
        } else if (vectorEncodedSequence != null) {
            Random random = new Random();
            int seqLength = vectorEncodedSequence.getColumnDimension();
            StringBuffer buf = new StringBuffer(seqLength);
            for (int j = 0; j < seqLength; j++) {
                double draw = random.nextDouble();
                int drawedBase = -1;
                double cummulativeProbability = 0;
                for (int i = 0; i < 4; i++) {
                    cummulativeProbability += vectorEncodedSequence.get(i, j);
                    if (draw <= cummulativeProbability) {
                        drawedBase = i;
                        break;
                    }
                }

                switch (drawedBase) {
                    case SHORT_ENCODED_A:
                        buf.append("A");
                        break;
                    case SHORT_ENCODED_C:
                        buf.append("C");
                        break;
                    case SHORT_ENCODED_G:
                        buf.append("G");
                        break;
                    case SHORT_ENCODED_T:
                        buf.append("T");
                        break;
                    default:
                        buf.append("-");
                        break;
                }
            }
            bases = buf.toString();
        }

        return bases;
    }

    public void setSequenceBases(String sequence) {
        this.sequenceBases = new StringBuilder(sequence);
    }

    public void appendToSequence(String partialSequence) {
        sequenceBases.append(partialSequence);
    }

    public void appendToSequence(char c) {
        sequenceBases.append(c);
    }

    public boolean isGap(int position) {
        boolean isGap = false;
        if (encodedSequence != null && encodedSequence.length >= position) {
            isGap = encodeIgnoreCase && (encodedSequence[position] == 4);
        } else if (vectorEncodedSequence != null && vectorEncodedSequence.getColumnDimension() >= position) {
            double maxProb = 0;

            for (int i = 0; i < 4; i++) {
                double baseProbability = vectorEncodedSequence.get(i, position);
                if (maxProb < baseProbability) {
                    maxProb = baseProbability;
                }
            }

            isGap = maxProb == 0;
        } else if (sequenceBases != null && sequenceBases.length() >= position) {
            isGap = '-' == sequenceBases.charAt(position);
        }

        return isGap;
    }

    public List<GenomicAnnotation> getSoftmaskedRegions() {
        List<GenomicAnnotation> softMaskedRegions = new ArrayList<GenomicAnnotation>();
        if (sequenceBases == null) {
            return softMaskedRegions;
        }

        Matcher m = SOFT_MASKED_PAT.matcher(sequenceBases);
        int i = 1;
        while (m.find()) {
            BasicGenomicAnnotation softMaskedReg = new BasicGenomicAnnotation(getId() + "_SoftmaskedReg_" + i++);
            softMaskedReg.setStart(m.start() + 1);
            softMaskedReg.setEnd(m.end() + 1);
            softMaskedReg.setChromosome(getId());
            softMaskedRegions.add(softMaskedReg);
        }
        return softMaskedRegions;
    }

    public short[] getEncodedSequence() {
        return encodedSequence;
    }

    public Matrix getVectorEncodedSequence() {
        return vectorEncodedSequence;
    }

    public Matrix encodeSequenceAsVector() {

        vectorEncodedSequence = new Matrix(4, sequenceBases.length());

        for (int j = 0; j < sequenceBases.length(); j++) {
            char c = sequenceBases.charAt(j);
            if ('a' == c || 'A' == c) {
                vectorEncodedSequence.set(SHORT_ENCODED_A, j, 1);
            } else if ('C' == c || 'c' == c) {
                vectorEncodedSequence.set(SHORT_ENCODED_C, j, 1);
            } else if ('G' == c || 'g' == c) {
                vectorEncodedSequence.set(SHORT_ENCODED_G, j, 1);
            } else if ('T' == c || 't' == c) {
                vectorEncodedSequence.set(SHORT_ENCODED_T, j, 1);
            }
        }

        return vectorEncodedSequence;
    }

    public short[] encodeSequence() {
        if (encodedSequence != null) {
            return encodedSequence;
        }

        short[] encodedSeq = new short[sequenceBases.length()];

        for (int i = 0; i < sequenceBases.length(); i++) {
            char c = sequenceBases.charAt(i);
            if ('c' == c) {
                encodedSeq[i] = SHORT_ENCODED_c;
            } else if ('C' == c) {
                encodedSeq[i] = SHORT_ENCODED_C;
            } else if ('G' == c) {
                encodedSeq[i] = SHORT_ENCODED_G;
            } else if ('g' == c) {
                encodedSeq[i] = SHORT_ENCODED_g;
            } else if ('a' == c) {
                encodedSeq[i] = SHORT_ENCODED_a;
            } else if ('A' == c) {
                encodedSeq[i] = SHORT_ENCODED_A;
            } else if ('t' == c) {
                encodedSeq[i] = SHORT_ENCODED_t;
            } else if ('T' == c) {
                encodedSeq[i] = SHORT_ENCODED_T;
            } else {
                encodedSeq[i] = SHORT_ENCODED_N;
            }
        }

        this.encodedSequence = encodedSeq;
        encodeIgnoreCase = false;
        return encodedSeq;
    }

    public short[] encodeSequenceIgnoreCase() {
        if (encodedSequence != null) {
            return encodedSequence;
        }

        short[] encodedSeq = new short[sequenceBases.length()];

        for (int i = 0; i < sequenceBases.length(); i++) {
            char c = sequenceBases.charAt(i);
            if ('a' == c || 'A' == c) {
                encodedSeq[i] = SHORT_ENCODED_A;
            } else if ('C' == c || 'c' == c) {
                encodedSeq[i] = SHORT_ENCODED_C;
            } else if ('G' == c || 'g' == c) {
                encodedSeq[i] = SHORT_ENCODED_G;
            } else if ('T' == c || 't' == c) {
                encodedSeq[i] = SHORT_ENCODED_T;
            } else if ('-' == c) {
                encodedSeq[i] = SHORT_ENCODED_GAP;
            } else {
                encodedSeq[i] = SHORT_ENCODED_GAP;  //All other letter codes are cosidered gaps.
            }
        }

        this.encodedSequence = encodedSeq;
        encodeIgnoreCase = true;
        return encodedSeq;
    }


    public List<Short> compute454Flow(char[] flowOrder) {
        return compute454Flow(flowOrder, 1, sequenceBases.length());
    }

    public List<Short> compute454Flow(char[] flowOrder, int start, int end) {
        char[] subSeqChrs = new char[end - start];
        sequenceBases.getChars(start - 1, end - 1, subSeqChrs, 0);
        ArrayList<Short> flow = new ArrayList<Short>();
        int seqIdx = 0;
        int flowIdx = 0;
        char flowLetter = 0;
        short flowRead = 0;
        while (seqIdx < subSeqChrs.length) {
            flowLetter = flowOrder[flowIdx % 4];
            char seqChar = Character.toUpperCase(subSeqChrs[seqIdx]);
            if (seqChar != 'A' && seqChar != 'T' && seqChar != 'G' && seqChar != 'C') {
                System.err.println("Sequence Character " + seqChar + " is not a base, ignoring it");
                seqIdx++;
            }
            if (flowLetter == seqChar) {
                flowRead++;
                seqIdx++;
            } else {
                flow.add(flowRead);
                flowIdx++;
                flowRead = 0;
            }

        }
        // TODO Auto-generated method stub
        return flow;
    }

    public float gcContent() {
        return computeGCContent(sequenceBases.toString());
    }

    public boolean contains(String motif) {
        Sequence revMotifSeq = new Sequence("rev_motif");
        revMotifSeq.setSequenceBases(motif);
        revMotifSeq.reverse();
        String revMotif = revMotifSeq.getSequenceBases();

        return getSequenceBases().contains(motif) || getSequenceBases().contains(revMotif);
    }

    public static float computeGCContent(String dnaString) {
        int gcs = 0;
        for (int i = 0; i < dnaString.length(); i++) {
            char c = Character.toUpperCase(dnaString.charAt(i));
            if ('C' == c || 'G' == c) {
                gcs++;
            }
        }

        return ((float) gcs) / ((float) dnaString.length());
    }

    public void reverse() {
        String seqBases = getSequenceBases();
        StringBuilder reversedSeq = new StringBuilder(seqBases.length());
        for (int j = seqBases.length() - 1; j >= 0; j--) {
            char c = seqBases.charAt(j);
            if ('c' == c) {
                reversedSeq.append('g');
            } else if ('C' == c) {
                reversedSeq.append('G');
            } else if ('G' == c) {
                reversedSeq.append('C');
            } else if ('g' == c) {
                reversedSeq.append('c');
            } else if ('a' == c) {
                reversedSeq.append('t');
            } else if ('A' == c) {
                reversedSeq.append('T');
            } else if ('t' == c) {
                reversedSeq.append('a');
            } else if ('T' == c) {
                reversedSeq.append('A');
            } else {
                reversedSeq.append(c);
            }
        }

        sequenceBases = reversedSeq;
        if (encodedSequence != null) {
            encodedSequence = null;
            encodeSequenceIgnoreCase();
        } else if (vectorEncodedSequence != null) {
            vectorEncodedSequence = null;
            encodeSequenceAsVector();
        }
        forwardStrand = false;
    }

    public static String reverseSequence(String seq) {
        Sequence tmpSeq = new Sequence("tmp");
        tmpSeq.setSequenceBases(seq);
        tmpSeq.reverse();
        return tmpSeq.getSequenceBases();
    }

    public void uppercase() {
        StringBuilder uppercasedSeq = new StringBuilder(sequenceBases.toString().toUpperCase());
        sequenceBases = uppercasedSeq;
    }

    /**
     * @param start - start of the region to extract the starting base will be included starting at 0
     * @param end   - the end of the region to extract, the base at position end will not be included.
     * @return
     */
    public SequenceRegion getRegion(int start, int end) {
        SequenceRegion region = new SequenceRegion(getId());
        region.setRegionStart(start);
        region.setRegionEnd(end);
        String bases = getSequenceBases();
        if (bases != null && bases.length() > 0) {
            region.setSequenceBases(bases.substring(start, end));
        }

        return region;
    }

    public void getRegion(SequenceRegion region) {
        region.setSequenceBases(sequenceBases.substring(region.getStart(), region.getEnd()));
    }

    public void getRegions(List<? extends SequenceRegion> regions) {
        Iterator<? extends SequenceRegion> it = regions.iterator();
        while (it.hasNext()) {
            SequenceRegion reg = it.next();
            getRegion(reg);
        }

    }

    public SequenceRegion extractRegionBasedOnGC(float targetGC, int size, int buffer) {
        SequenceRegion theRegion = null;
        float closestGC = 0;

        System.out.println("extracting region based on GC for sequence " + getId() + " of size " + getSequenceBases().length() + " buffer " + buffer);

        if (getSequenceBases() == null || getSequenceBases().length() < buffer) {
            return theRegion;
        }

        for (int i = buffer; i < (getSequenceBases().length() - buffer - size); i += buffer) {
            String currentSequence = getSequenceBases().substring(i - buffer, i + size + buffer - 1);
            float currentGC = computeGCContent(currentSequence);
            System.out.println("position " + i + " currentGC " + currentGC + ", targetGC " + targetGC + ",  closestGC distance to target " + Math.abs(closestGC - targetGC) + " currentGC distance to target " + Math.abs(targetGC - currentGC));
            if (Math.abs(closestGC - targetGC) > Math.abs(targetGC - currentGC)) {
                closestGC = currentGC;
                theRegion = new SequenceRegion(getId());
                theRegion.setSequenceBases(currentSequence);
                theRegion.setRegionStart(i - buffer);
                theRegion.setEnd(i + size + buffer - 1);
            }
        }
        System.out.println("Returning closest GC region: " + theRegion + " GC% " + closestGC);
        return theRegion;
    }

    public boolean isForwardStrand() {
        return forwardStrand;
    }

    public StringBuilder getSequenceBuilder() {
        return sequenceBases == null || sequenceBases.length() == 0
                ? new StringBuilder(getSequenceBases())
                : sequenceBases;
    }

    public int getGapsSize() {
        int totalGaps = 0;
        if (sequenceBases != null && sequenceBases.length() > 0) {
            char[] baseArray = sequenceBases.toString().toUpperCase().toCharArray();
            for (int i = 0; i < baseArray.length; i++) {
                if (baseArray[i] == 'N') {
                    totalGaps++;
                }
            }
        }
        return totalGaps;
    }

    public WindowSlider getSlider(int windowSize, int overlap) {
        return WindowSlider.getSlider(this, windowSize, overlap);
    }

    public List<SequenceRegion> chunk(int chunkSize, int chunkOverlap) {
        int numOfChunks = (int) Math.floor(getLength() / ((float) (chunkSize - chunkOverlap)));
        int chunkStart = 0;
        int chunkEnd = 0;
        System.out.println("\tChunking " + getId() + " number of chunks " + numOfChunks);
        List<SequenceRegion> chunks = new ArrayList<SequenceRegion>(numOfChunks);

        for (int i = 0; i < numOfChunks; i++) {
            chunkStart = i * (chunkSize - chunkOverlap);
            chunkEnd = chunkStart + chunkSize;
            SequenceRegion chunk = getRegion(chunkStart, chunkEnd);
            chunks.add(chunk);
        }
        if (chunkEnd < getLength()) {
            chunkStart = chunkEnd - chunkOverlap;
            SequenceRegion lastChunk = getRegion(chunkStart, getLength());
            chunks.add(lastChunk);
        }

        System.out.println("\tObtained " + chunks.size() + " chunks");
        return chunks;
    }

    public int getEnd() {
        return getLength();
    }

    /**
     * Finds the list of ungapped regions in the sequence.
     *
     * @return A list of two sized integer arrays with the start and end of the ungapped regions. The convention here is semi closed
     *         intervals, each list items is of the form [start, end).
     */
    public List<int[]> findUngappedSequenceChunks() {
        Matcher m = UNGAPPED_SEQUENCE_PAT.matcher(getSequenceBases());
        ArrayList<int[]> gapUngapp = new ArrayList<int[]>();
        //System.out.println("Sequence for " + getId() + ": " + getSequenceBases());
        while (m.find()) {
            int[] ungappedReg = {m.start(), m.end()};
            //System.out.println("\treg["+m.start()+"-"+m.end()+"]: " + getSequenceBases().substring(m.start(), m.end()));
            gapUngapp.add(ungappedReg);
        }
        //System.out.println("");

        return gapUngapp;
    }


}


