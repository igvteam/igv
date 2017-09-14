/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broad.igv.feature.genome.fasta;

import org.broad.igv.feature.genome.ChromosomeNameComparator;
import org.broad.igv.feature.genome.Sequence;

import java.io.IOException;
import java.util.*;

/**
 * Implementation of Sequence backed by an indexed fasta file
 *
 * @author Jim Robinson
 * @date 3/24/12
 */

public class FastaDirectorySequence implements Sequence {

    Map<String, FastaIndexedSequence> sequenceMap;
    List<String> chromosomeNames;
    Map<String, Integer> chrLengths;

    public FastaDirectorySequence(String directoryPath, String[] fastaFiles) throws IOException {
        readIndexes(directoryPath, fastaFiles);
    }

    private void readIndexes(String directoryPath, String[] fastaFiles) throws IOException {
        sequenceMap = new LinkedHashMap<>();
        for (String file : fastaFiles) {
            String fastaPath = directoryPath + "/" + file;
            FastaIndexedSequence fastaSequence = new FastaIndexedSequence(fastaPath);
            for (String chr : fastaSequence.getChromosomeNames()) {
                sequenceMap.put(chr, fastaSequence);
            }
        }

        chromosomeNames = new ArrayList<>();
        for (FastaIndexedSequence fastaSequence : getFastaSequences()) {
            chromosomeNames.addAll(fastaSequence.getChromosomeNames());
        }
        Collections.sort(chromosomeNames, ChromosomeNameComparator.get());

        chrLengths = new HashMap<>(chromosomeNames.size());
        for (FastaIndexedSequence fastaSequence : getFastaSequences()) {
            for (String chr : fastaSequence.getChromosomeNames()) {
                int length = fastaSequence.getChromosomeLength(chr);
                chrLengths.put(chr, length);
            }
        }
    }

    public Collection<FastaIndexedSequence> getFastaSequences() {
        return sequenceMap.values();
    }

    public byte[] getSequence(String chr, int start, int end, boolean useCache) {


        if (!sequenceMap.containsKey(chr)) {
            return null;
        }
        return sequenceMap.get(chr).getSequence(chr, start, end, useCache);
    }


    @Override
    public List<String> getChromosomeNames() {
        return chromosomeNames;
    }

    @Override
    public byte getBase(String chr, int position) {
        throw new RuntimeException("getBase() is not implemented for class " + FastaIndexedSequence.class.getName());
    }

    @Override
    public int getChromosomeLength(String chrname) {
        return chrLengths.get(chrname);
    }

    @Override
    public boolean isRemote() {

        for(FastaIndexedSequence seq : sequenceMap.values()) {
            if(seq.isRemote()) return true;
        }
        return false;
    }


}
