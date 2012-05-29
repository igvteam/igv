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

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Stack;

//import Jama.Matrix;
import edu.mit.broad.prodinfo.chromosome.BasicGenomicAnnotation;
import edu.mit.broad.prodinfo.datastrutures.IntervalTree;
import edu.mit.broad.prodinfo.datastrutures.IntervalTree.Node;
import edu.mit.broad.prodinfo.genomicplot.GenomicAnnotation;
import edu.mit.broad.prodinfo.genomicplot.ParseException;
import edu.mit.broad.prodinfo.sequence.SequenceRegion;
import org.broad.igv.util.ParsingUtils;
import org.broad.tribble.util.SeekableStream;
import org.broad.tribble.util.SeekableStreamFactory;

public class MAFAlignment extends MultipleAlignment {
    private IntervalTree<MAFMultipleAlignmentBlock> alignmentBlockTree;
    private List<String> sequenceIds;


    MAFHeader header;
    //Map<Integer, Long> index;
    IntervalTree<Long> index;
    private String referenceChromosome;

    public MAFAlignment() {
        super();
        header = new MAFHeader();
        alignmentBlockTree = new IntervalTree<MAFMultipleAlignmentBlock>();
        sequenceIds = new ArrayList<String>();
    }

    public MAFAlignment(IntervalTree<Long> index) throws IOException, ParseException {
        this();
        this.index = index;
        //System.out.println("Loaded index file " + index.keySet().size());
    }

    public MAFHeader getHeader() {
        return header;
    }

    // MultipleAlignment overloading.
    public void encode() {
        //Iterator<MAFMultipleAlignmentBlock> alnBlockIt = alignmentBlocks.iterator();
        Iterator<Node<MAFMultipleAlignmentBlock>> alnBlockIt = alignmentBlockTree.iterator();
        while (alnBlockIt.hasNext()) {
            //MultipleAlignment block = alnBlockIt.next();
            Node<MAFMultipleAlignmentBlock> blockNode = alnBlockIt.next();
            blockNode.getValue().encode();
        }

    }

    public void encodeAsMatrix() {
        //Iterator<MAFMultipleAlignmentBlock> alnBlockIt = alignmentBlocks.iterator();
        Iterator<Node<MAFMultipleAlignmentBlock>> alnBlockIt = alignmentBlockTree.iterator();
        while (alnBlockIt.hasNext()) {
            //MultipleAlignment block = alnBlockIt.next();
            Node<MAFMultipleAlignmentBlock> blockNode = alnBlockIt.next();
            blockNode.getValue().encodeAsMatrix();
        }
    }

    public void reverse() {
        //Iterator<MAFMultipleAlignmentBlock> alnBlockIt = alignmentBlocks.iterator();
        Iterator<Node<MAFMultipleAlignmentBlock>> alnBlockIt = alignmentBlockTree.iterator();
        while (alnBlockIt.hasNext()) {
            //MultipleAlignment block = alnBlockIt.next();
            Node<MAFMultipleAlignmentBlock> blockNode = alnBlockIt.next();
            blockNode.getValue().reverse();
        }
        /*
          Stack<MAFMultipleAlignmentBlock> reversedBlocks = new Stack<MAFMultipleAlignmentBlock>();
          while(alignmentBlocks.size() > 0) {
              MAFMultipleAlignmentBlock block = alignmentBlocks.pop();
              block.reverse();
              reversedBlocks.push(block);
          }
          alignmentBlocks = reversedBlocks;
          */
    }

    public int length() {
        int length = 0;
        //if(!alignmentBlocks.isEmpty()) {
        if (alignmentBlockTree.size() > 0) {
            MAFMultipleAlignmentBlock first = alignmentBlockTree.min().getValue();//alignmentBlocks.get(0);
            MAFMultipleAlignmentBlock last = alignmentBlockTree.max().getValue();//alignmentBlocks.get(alignmentBlocks.size() - 1);

            length = last.getReferenceEnd() - first.getReferenceStart();
        }

        return length;
    }

    /**
     * Assumes that encode has been called and each aligned sequence is encoded
     */
    public Map<String, Short> getColumn(int i) {
        Map<String, Short> col = getGapColumn();
        Iterator<Node<MAFMultipleAlignmentBlock>> overlappingNodeIt = alignmentBlockTree.overlappers(i, i + 1);

        if (overlappingNodeIt.hasNext()) {
            MAFMultipleAlignmentBlock containingBlock = overlappingNodeIt.next().getValue();
            Map<String, Short> closestCol = containingBlock.getColumn(i);
            Iterator<String> closestAlignSeqIdIt = closestCol.keySet().iterator();
            while (closestAlignSeqIdIt.hasNext()) {
                String seqId = closestAlignSeqIdIt.next();
                col.put(seqId, closestCol.get(seqId));
            }
        }

        return col;
    }

    public void addShortEncodedColumn(Map<String, Short> col) {
        throw new RuntimeException("Not yet implemented, to edit an alignment use a base MultipleAlignment MAF are ReadOnly");
    }

    public void addShortEncodedRegion(Map<String, short[]> region) {
        throw new RuntimeException("Not yet implemented, to edit an alignment use a base MultipleAlignment MAF are ReadOnly");
    }

    public void addSequence(AlignedSequence seq) {
        throw new RuntimeException("Not yet implemented, to edit an alignment use a base MultipleAlignment MAF are ReadOnly");
    }

    public void addSequences(List<AlignedSequence> sequences) {
        throw new RuntimeException("Not yet implemented, to edit an alignment use a base MultipleAlignment MAF are ReadOnly");
    }

    public boolean isEmpty() {
        return alignmentBlockTree.isEmpty();
    }

    public void write(BufferedWriter bw) throws IOException {
        if (getIOHelper() != null) {
            super.write(bw);
        } else {
            header.addVariableValuePair("ref", getReferenceId());
            header.write(bw);
            Iterator<Node<MAFMultipleAlignmentBlock>> blockIt = alignmentBlockTree.iterator();
            while (blockIt.hasNext()) {
                blockIt.next().getValue().write(bw);
            }
        }

    }

    public MAFAlignment getSubAlignment(int refStart, int refEnd, boolean reverse) {
        MAFAlignment subAln = new MAFAlignment();
        subAln.setReferenceId(getReferenceId());
        subAln.sequenceIds = this.sequenceIds;
        Iterator<MAFMultipleAlignmentBlock> overlapperIt =
                new IntervalTree.ValuesIterator<MAFMultipleAlignmentBlock>(alignmentBlockTree.overlappers(refStart, refEnd));
        GenomicAnnotation target = new BasicGenomicAnnotation(getReference());
        target.setStart(refStart);
        target.setEnd(refEnd);
        while (overlapperIt.hasNext()) {
            MAFMultipleAlignmentBlock overlappingBlock = overlapperIt.next();
            subAln.addBlock(overlappingBlock.trim(target));
        }
        return subAln;
    }

    public AlignedSequence getAlignedSequence(String sequenceId, boolean fillInGapsBetweenBlocks) {
        //System.out.println("Is alignmentBlockTree empty? " + alignmentBlockTree.isEmpty());
        if (!sequenceIds.contains(sequenceId) || alignmentBlockTree.isEmpty()) {
            //System.out.println("seqids " + sequenceIds + " does not contain <" + sequenceId +"> is contained in sequenceIds " + sequenceIds.contains(sequenceId));
            return null;
        }

        AlignedSequence seq = new AlignedSequence(sequenceId);
        seq.setId(sequenceId);
        if (getReferenceId().equals(sequenceId)) {
            //System.out.println("Dealing with ref seq");
            seq.setStart(alignmentBlockTree.min().getValue().getReferenceStart());
            seq.setChromosome(alignmentBlockTree.min().getValue().getAlignedSequence(sequenceId).getChromosome());
            seq.setEnd(alignmentBlockTree.max().getValue().getReferenceEnd());
        }
        Iterator<MAFMultipleAlignmentBlock> it = new IntervalTree.ValuesIterator<MAFMultipleAlignmentBlock>(alignmentBlockTree.iterator());
        MAFMultipleAlignmentBlock lastBlock = null;
        while (it.hasNext()) {
            MAFMultipleAlignmentBlock block = it.next();
            AlignedSequence segment = block.getAlignedSequence(sequenceId);
            if (lastBlock != null && lastBlock.getReferenceEnd() < block.getReferenceStart() && fillInGapsBetweenBlocks) {
                //System.out.println("last block " + (lastBlock != null ? lastBlock.getReferenceStart()+"-"+lastBlock.getReferenceEnd():" not yet ") + " new block " + block.getReferenceStart()+"-"+block.getReferenceEnd() + ", appending " + (block.getReferenceStart() - lastBlock.getReferenceEnd()) + " Ns");
                for (int i = 0; i < block.getReferenceStart() - lastBlock.getReferenceEnd(); i++) {
                    //System.out.println("\t\tappending gaps N");
                    seq.appendToSequence('N');
                }
            }
            seq.appendToSequence(segment.getSequenceBases());
            /*
               if(segment != null) {
                   seq.appendToSequence(segment.getSequenceBases());
               } else {
                   for(int i = 0; i < blockNode.getValue().length(); i++) {
                       seq.appendToSequence('-');
                   }
               }
               */
            lastBlock = block;
        }
        return seq;
    }

    public AlignedSequence getAlignedSequence(String sequenceId) {
        return getAlignedSequence(sequenceId, true);
    }

    public List<AlignedSequence> getAlignedSequences(boolean fillGapsBetweenBlocks) {
        ArrayList<AlignedSequence> sequences = new ArrayList<AlignedSequence>(sequenceIds.size());
        Iterator<String> idIt = sequenceIds.iterator();
        while (idIt.hasNext()) {
            String id = idIt.next();
            sequences.add(getAlignedSequence(id, fillGapsBetweenBlocks));
        }
        return sequences;
    }

    public List<AlignedSequence> getAlignedSequences() {
        return getAlignedSequences(true);
    }

    public boolean overlaps(GenomicAnnotation annotation) {
        return index.overlappers(annotation.getStart(), annotation.getEnd()).hasNext();
    }

    public boolean contains(GenomicAnnotation annotation) {
        return index.min().getStart() < annotation.getStart() && index.max().getEnd() > annotation.getEnd();
    }

    public List<String> getAlignedSequenceIds() {
        return sequenceIds;
    }

    public int getReferenceStart() {
        Node<MAFMultipleAlignmentBlock> first = alignmentBlockTree.min();
        return first == null ? -1 : first.getValue().getReferenceStart();
    }

    public int getReferenceEnd() {
        Node<MAFMultipleAlignmentBlock> last = alignmentBlockTree.max();
        return last == null ? -1 : last.getValue().getReferenceEnd();
    }

    public void compress() {
        Iterator<Node<MAFMultipleAlignmentBlock>> nodeIt = alignmentBlockTree.iterator();
        while (nodeIt.hasNext()) {
            nodeIt.next().getValue().compress();
        }
    }

    public MultipleAlignment toMultipleAlignment(boolean fillGapsBetweenBlocks) {
        MultipleAlignment ma = new MultipleAlignment();
        ma.setReferenceId(getReferenceId());
        ma.setRefGapped(true);
        ma.addSequences(getAlignedSequences(fillGapsBetweenBlocks));
        ma.getReference().setStart(getReferenceStart());
        ma.getReference().setEnd(getReferenceEnd());
        //System.out.println("to multiplealn #aln seqs: " + ma.length() + " ref? " + ma.getReferenceId());
        return ma;
    }

    public MultipleAlignment toMultipleAlignment() {
        return toMultipleAlignment(true);
    }

    public MultipleAlignment sampleColumns(int colNum, int numConsecutiveCols) {
        Random r = new Random();
        int treeSize = alignmentBlockTree.size();
        MultipleAlignment sampledAlignment = new MultipleAlignment();
        sampledAlignment.setReferenceId(getReferenceId());
        for (int i = 0; i < colNum - numConsecutiveCols + 1; i = numConsecutiveCols + i) {
            int blockId = r.nextInt(treeSize);
            System.out.println("randomly selected block " + blockId + " out of " + treeSize + " blocks");
            MAFMultipleAlignmentBlock block = alignmentBlockTree.findByIndex(blockId).getValue();
            MultipleAlignment sampledCols = block.sampleColumns(1, numConsecutiveCols);
            Iterator<String> seqIdIt = getAlignedSequenceIds().iterator();
            while (seqIdIt.hasNext()) {
                String seqId = seqIdIt.next();
                if (sampledCols.getAlignedSequence(seqId) == null) {
                    AlignedSequence missingSeq = new AlignedSequence(seqId);
                    for (int j = 0; j < numConsecutiveCols; j++) {
                        missingSeq.appendToSequence('-');
                    }
                    sampledCols.addSequence(missingSeq);
                }
            }
            sampledAlignment.append(sampledCols); //changed without test 5/23/08
        }

        return sampledAlignment;
    }

    public void load(SeekableStream handle, int referenceStart, int referenceEnd, List<String> sequencesToLoad) throws IOException, ParseException {


        long offset = getClosestOffset(referenceStart);
        handle.seek(offset);
        boolean okToAdd = true;
        BasicGenomicAnnotation reference = new BasicGenomicAnnotation("reference");
        reference.setStart(referenceStart);
        reference.setEnd(referenceEnd);
        sequenceIds = new ArrayList<String>();

        if (sequencesToLoad != null) {
            Iterator<String> seqIt = sequencesToLoad.iterator();
            while (seqIt.hasNext()) {
                sequenceIds.add(seqIt.next());
            }
        }

        String line = null;
        String[] lineInfo = null;
        Stack<MAFMultipleAlignmentBlock> alignmentBlockStack = new Stack<MAFMultipleAlignmentBlock>();

        BufferedReader reader = new BufferedReader(new InputStreamReader(handle));

        while ((line = reader.readLine()) != null) {
            //Ignore all comment lines
            if (line.startsWith("#") || line.trim().length() == 0) {
                continue;
            }

            //System.out.println(line);
            if (line.startsWith("a ")) {
                //System.err.println("New alignment: " + line);
                //First check last alignment to see if it should be kept.
                MAFMultipleAlignmentBlock lastMA = alignmentBlockStack.isEmpty() ? null : alignmentBlockStack.pop();
                if (lastMA != null) {
                    alignmentBlockTree.put(lastMA.getReferenceStart(), lastMA.getReferenceEnd(), lastMA.trim(reference));
                }
                /*
                    if(lastMA != null && referenceEnd > 0 && !lastMA.overlaps(reference)) {
                        alignments.pop(); // multiple alignment just built does overlap desired region.
                        AlignedSequence refAln = lastMA.getAlignedSequence(reference.getContainingSequenceId());
                        System.out.println("Reference seq<"+reference.getContainingSequenceId() +"> lastMA aligned seqs : " + lastMA.getAlignedSequenceIds());
                        if (refAln.getStart() > referenceEnd) {
                            break;
                        }
                    }
                    */
                MAFMultipleAlignmentBlock ma = new MAFMultipleAlignmentBlock();
                alignmentBlockStack.push(ma);
                lineInfo = line.substring(2).split("\\s");
                ma.setAlignmentInfoFromRawData(lineInfo);
                okToAdd = true;
            } else if (line.startsWith("s ")) {
                if (alignmentBlockStack.isEmpty() || !okToAdd) {
                    continue;
                }
                //System.err.println("\tAlignment aligned seq " + line);
                MAFMultipleAlignmentBlock ma = alignmentBlockStack.peek();
                lineInfo = line.substring(2).split("\\s+");
                AlignedSequence seq = ma.createSequenceFromRawData(lineInfo);

                if (getReferenceId() == null) {
                    setReferenceId(seq.getId());
                    setReferenceChromosome(seq.getChromosome());
                }

                if (ma.getReferenceId() == null) {
                    ma.setReferenceId(seq.getId());
                }

                if (sequencesToLoad == null || sequencesToLoad.size() == 0 || sequencesToLoad.contains(seq.getId())) {
                    //System.err.println("YES: sequence " + seq.getId() + " chr " + seq.getChromosome() + " sequenceToLoad is null? " + (sequencesToLoad == null) + " is empty? " + (sequencesToLoad.size() == 0 ) + " or contains sequence " + sequencesToLoad.contains(seq.getId()));
                    ma.addSequence(seq);
                } else {
                    //System.err.println("NO: sequence " + seq.getId() + " chr " + seq.getChromosome() + " sequenceToLoad is null? " + (sequencesToLoad == null) + " is empty? " + (sequencesToLoad.size() == 0 ) + " or contains sequence " + sequencesToLoad.contains(seq.getId()));
                    continue;
                }

                //A bit wasteful, and should fix, but finally check if reference location is OK to add
                if (getReferenceId().equals(seq.getId()) && seq.getEnd() <= referenceStart) {
                    alignmentBlockStack.pop();
                    okToAdd = false;
                } else if (getReferenceId().equals(seq.getId()) && seq.getStart() >= referenceEnd) {
                    alignmentBlockStack.pop();
                    break;
                } else {
                    //alignmentBlockStack.push(ma);
                    if (!getAlignedSequenceIds().contains(seq.getId())) {
                        sequenceIds.add(seq.getId());
                    }
                }

            } else if (line.startsWith("i ")) {
                //We do not handle information lines yet.
                continue;
            } else if (line.startsWith("q ")) {
                //We do not handle quality lines yet.
                continue;
            } else if (line.startsWith("e ")) {
                //We do not support e lines yet.
                continue;
            } else {
                throw new ParseException("Invalid alignment line <" + line + ">");
            }

        }
        //Handle last alignment block
        MAFMultipleAlignmentBlock lastMA = alignmentBlockStack.isEmpty() ? null : alignmentBlockStack.pop();
        if (lastMA != null && lastMA.getReferenceStart() < referenceEnd && lastMA.getReferenceEnd() > referenceStart) {
            alignmentBlockTree.put(lastMA.getReferenceStart(), lastMA.getReferenceEnd(), lastMA.trim(reference));
            //System.err.println("put last alignment: ");
        }
    }

    public IntervalTree<Long> getIndex() {
        return index;
    }

    public void setBlocks(List<MAFMultipleAlignmentBlock> blocks) {
        alignmentBlockTree = new IntervalTree<MAFMultipleAlignmentBlock>();
        Iterator<MAFMultipleAlignmentBlock> blockIt = blocks.iterator();
        while (blockIt.hasNext()) {
            addBlock(blockIt.next());
        }

    }

    public void addBlock(MAFMultipleAlignmentBlock block) throws IllegalArgumentException {
        Iterator<Node<MAFMultipleAlignmentBlock>> overlappers =
                alignmentBlockTree.overlappers(block.getReferenceStart(), block.getReferenceEnd());
        if (overlappers != null && overlappers.hasNext()) {
            throw new IllegalArgumentException("A block in the alignment to append already overlaps exisiting alignment block. Can't append");
        }
        Iterator<String> alnSeqIdIt = block.getAlignedSequenceIds().iterator();
        while (alnSeqIdIt.hasNext()) {
            String seqId = alnSeqIdIt.next();
            if (!sequenceIds.contains(seqId)) {
                sequenceIds.add(seqId);
            }
        }

        alignmentBlockTree.put(block.getReferenceStart(), block.getReferenceEnd(), block);
    }

    public void addBlocks(List<MAFMultipleAlignmentBlock> blocks) throws IllegalArgumentException {
        addBlocks(blocks.iterator());
    }

    public void addBlocks(Iterator<MAFMultipleAlignmentBlock> blockIt) throws IllegalArgumentException {
        while (blockIt.hasNext()) {
            addBlock(blockIt.next());
        }
    }

    public void append(MAFAlignment mafAln) throws IllegalArgumentException {
        addBlocks(mafAln.getBlockIterator());
    }

    public void clear() {
        alignmentBlockTree = new IntervalTree<MAFMultipleAlignmentBlock>();
    }

    long getClosestOffset(int position) {
        Node<Long> node = index.max(position, position + 1);
        return node == null ? 0 : node.getValue();
    }

    public Iterator<MAFMultipleAlignmentBlock> getBlockIterator() {
        final Iterator<Node<MAFMultipleAlignmentBlock>> nodeIt = alignmentBlockTree.iterator();
        return new Iterator<MAFMultipleAlignmentBlock>() {

            public boolean hasNext() {
                return nodeIt.hasNext();
            }

            public MAFMultipleAlignmentBlock next() {
                return nodeIt.next().getValue();
            }

            public void remove() {
                nodeIt.remove();
            }

        };


    }

    private void setReferenceChromosome(String chromosome) {
        this.referenceChromosome = chromosome;

    }

    private String getReferenceChromosome() {
        return this.referenceChromosome;
    }

    private void padN(int start, int end) {
        if (getReference() == null) {
            return;
        }

        int startPadLength = getReferenceStart() - start;
        int endPadLength = end - getReferenceEnd();
        if (startPadLength > 0) {
            MAFMultipleAlignmentBlock pad = new MAFMultipleAlignmentBlock();
            pad.setReferenceId(getReferenceId());
            AlignedSequence ref = new AlignedSequence(getReference().getContainingSequenceId());
            ref.setStart(start);
            ref.setEnd(getReferenceStart());
            ref.setId(getReferenceId());
            pad.addAlignedSequence(ref);


            ref.setSequenceBases(buildPad(startPadLength));
            addBlock(pad);
        }

        if (endPadLength > 0) {
            MAFMultipleAlignmentBlock pad = new MAFMultipleAlignmentBlock();
            pad.setReferenceId(getReferenceId());
            AlignedSequence ref = new AlignedSequence(getReference().getContainingSequenceId());
            ref.setStart(getReferenceEnd() + 1);
            ref.setEnd(end);
            ref.setId(getReferenceId());
            pad.addAlignedSequence(ref);


            ref.setSequenceBases(buildPad(endPadLength));
            addBlock(pad);
        }

    }

    private String buildPad(int size) {
        StringBuilder padSeq = new StringBuilder(size);
        for (int i = 0; i < size; i++) {
            padSeq.append("N");
        }

        return padSeq.toString();
    }

    public static class MAFMultipleAlignmentBlock extends MultipleAlignment {
        int pass;

        public MAFMultipleAlignmentBlock() {
            super();
            setRefGapped(true);
        }

        public MAFMultipleAlignmentBlock(MultipleAlignment ma) {
            this();
            setReferenceId(ma.getReferenceId());
            Iterator<AlignedSequence> it = ma.getAlignedSequences().iterator();
            while (it.hasNext()) {
                addAlignedSequence(it.next());
            }

        }

        public int getPass() {
            return pass;
        }


        public AlignedSequence createSequenceFromRawData(String[] data) {
            String[] seqNameInfo = data[0].split("\\.");
            AlignedSequence aln = new AlignedSequence(seqNameInfo[0].intern());
            if (seqNameInfo.length > 1) {
                aln.setChromosome(seqNameInfo[1]);
            }
            aln.setId(seqNameInfo[0]);
            aln.setName(seqNameInfo[0]);
            aln.setRegionStart(Integer.parseInt(data[1]));
            aln.setRegionEnd(aln.getRegionStart() + Integer.parseInt(data[2]));
            aln.setStrand(data[3]);
            aln.setTotalLength(Integer.parseInt(data[4]));
            aln.setSequenceBases(data[5]);
            return aln;
        }

        public AlignedSequence addSequenceFromRawData(String[] data) {
            AlignedSequence aln = createSequenceFromRawData(data);
            addSequence(aln);
            return aln;
        }

        public void addAlignedSequence(AlignedSequence aln) {
            addSequence(aln);
        }

        public AlignedSequence getAlignedSequence(String sequenceId) {
            AlignedSequence seq = super.getAlignedSequence(sequenceId);
            if (seq == null) {
                seq = new AlignedSequence(sequenceId);
                int length = length();
                for (int i = 0; i < length; i++) {
                    seq.appendToSequence('-');
                }
            }

            return seq;
        }

        public void setAlignmentInfoFromRawData(String[] data) {
            for (int i = 0; i < data.length; i++) {
                String[] nameValPair = data[i].split("=");
                if (nameValPair[0].equalsIgnoreCase("pass")) {
                    setPass(Integer.parseInt(nameValPair[1]));
                } else if (nameValPair[0].equalsIgnoreCase("score")) {
                    setScore(Float.parseFloat(nameValPair[1]));
                } else {
                    System.err.println("Unsuported alignment attribute: " + nameValPair[0] + " .... ignoring");
                }
            }
        }

        public void setPass(int pass) {
            this.pass = pass;
        }

        public MAFMultipleAlignmentBlock trim(GenomicAnnotation target) {
            //System.out.println("target: " + target.getLocationString() + " reference: " + getReference().getLocationString());
            if (target.getStart() <= getReferenceStart() && target.getEnd() >= getReferenceEnd()) {
                //System.out.println("No trim is nencesary");
                return this;
            }

            return new MAFMultipleAlignmentBlock(getSubAlignment(target.getStart(), target.getEnd(), false));

        }


    }

    public static class MAFHeader {
        MAFAlignment alignment;
        String version;
        String scoring;
        String program;

        String runParameters;

        Hashtable<String, String> otherVariables = new Hashtable<String, String>();

        private static final long serialVersionUID = 239451013421586L;

        protected void setVariablesFromRawData(String[] data) {
            for (int i = 0; i < data.length; i++) {
                String[] variableValuePair = data[i].split("=");
                if (variableValuePair[0].equalsIgnoreCase("version")) {
                    version = variableValuePair[1];
                } else if (variableValuePair[0].equalsIgnoreCase("scoring")) {
                    scoring = variableValuePair[1];
                } else if (variableValuePair[0].equalsIgnoreCase("program")) {
                    program = variableValuePair[1];
                } else {
                    otherVariables.put(variableValuePair[0].toLowerCase(), variableValuePair[1]);
                }

            }
        }

        public String getRunParameters() {
            return runParameters;
        }

        public void addVariableValuePair(String variable, String value) {
            otherVariables.put(variable, value);
        }

        public void write(BufferedWriter bw) throws IOException {
            bw.write("##maf");
            if (version != null && version.length() > 0) {
                bw.write(" version=");
                bw.write(version);
            }
            if (scoring != null && scoring.length() > 0) {
                bw.write(" scoring=");
                bw.write(scoring);
            }
            if (program != null && program.length() > 0) {
                bw.write(" program=");
                bw.write(program);
            }

            Iterator<String> varIt = otherVariables.keySet().iterator();
            while (varIt.hasNext()) {
                String var = varIt.next();
                bw.write(" ");
                bw.write(var);
                bw.write("=");
                bw.write(otherVariables.get(var));
            }
            bw.newLine();

            if (runParameters != null) {
                bw.write("# ");
                bw.write(runParameters);
            }

            bw.newLine();
            bw.newLine();
        }

        public void setRunParameters(String runParameters) {
            this.runParameters = runParameters;
        }


        public void setProgram(String program) {
            this.program = program;
        }

        public void setScoring(String scoring) {
            this.scoring = scoring;
        }

        public void setVersion(String version) {
            this.version = version;
        }

        public String getProgram() {
            return program;
        }

        public String getScoring() {
            return scoring;
        }

        public String getVersion() {
            return version;
        }

    }

    protected void setAlignedSequences(List<String> seqIds) {
        sequenceIds = seqIds;
    }


    private Map<String, Short> getGapColumn() {
        HashMap<String, Short> gapCol = new HashMap<String, Short>(sequenceIds.size());
        Iterator<String> sequenceIdIt = sequenceIds.iterator();
        while (sequenceIdIt.hasNext()) {
            gapCol.put(sequenceIdIt.next(), GAP_CODE);
        }

        return gapCol;
    }

}
