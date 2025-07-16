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

import htsjdk.samtools.seekablestream.SeekableStream;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.genome.Sequence;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.collections.ByteArrayList;
import org.broad.igv.util.stream.IGVSeekableStreamFactory;

import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Implementation of Sequence backed by an indexed fasta file
 *
 * @author jrobinso
 * @date 8/7/11
 */
public class FastaNonIndexedSequence implements Sequence {

    static Logger log = LogManager.getLogger(FastaNonIndexedSequence.class);


    private ArrayList<String> chromoNamesList;
    private Map<String, byte[]> sequenceMap;
    private Map<String, Integer> chromSizesMap;
    private List<Chromosome> chromosomes;

    public FastaNonIndexedSequence(String path) throws IOException {
        read(path);
    }

    private void read(String path) throws IOException {
        try (BufferedReader reader = ParsingUtils.openBufferedReader(path)) {

            String line;
            String lastName = null;
            sequenceMap = new HashMap<>();
            chromoNamesList = new ArrayList<>();
            chromSizesMap = new HashMap<>();
            chromosomes = new ArrayList<>();
            int idx = 0;
            ByteArrayList byteArrayList = new ByteArrayList(100000);
            while ((line = reader.readLine()) != null) {
                if (line.startsWith(">")) {
                    String name = line.substring(1).trim();
                    chromoNamesList.add(name);
                    if (lastName != null) {
                        final int size = byteArrayList.size();
                        sequenceMap.put(lastName, byteArrayList.toArray());
                        chromSizesMap.put(lastName, size);
                        chromosomes.add(new Chromosome(idx++, lastName, size));
                    }
                    lastName = name;
                    byteArrayList = new ByteArrayList(100000);
                } else {
                    byteArrayList.addAll(line.getBytes());
                }
            }
        }
    }


    /**
     * Return the sequence for the query interval as a byte array.  Coordinates are "ucsc" style (0 based)
     * @param chr
     * @param qstart
     * @param qend
     * @return
     */

    public byte[] getSequence(String chr, int qstart, int qend) {

        byte [] seq = sequenceMap.get(chr);
        if (seq == null) {
            log.warn("No fasta sequence entry for: " + chr);
            return null;
        }
        if(qstart < 0 || qstart >= seq.length) {
            throw new IllegalArgumentException("Start position out of bounds: " + qstart);
        }
        if(qend <= qstart || qend >= seq.length) {
            throw new IllegalArgumentException("End position out of bounds: " + qend);
        }
        return java.util.Arrays.copyOfRange(seq, qstart, qend);
    }


    @Override
    public byte getBase(String chr, int position) {
        throw new RuntimeException("getBase() is not implemented for class " + FastaNonIndexedSequence.class.getName());
    }


    @Override
    public List<String> getChromosomeNames() {
        return chromoNamesList;
    }

    @Override
    public int getChromosomeLength(String chrname) {
        return chromSizesMap.get(chrname);
    }

    @Override
    public List<Chromosome> getChromosomes() {
        return chromosomes;
    }
}
