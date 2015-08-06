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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.feature.genome;

import org.apache.log4j.Logger;
import org.broad.igv.util.ObjectCache;

import java.util.Hashtable;
import java.util.List;

/**
 * A wrapper class that provides caching for on-disk, queried, and web-service Sequence implementations.
 *
 * @author jrobinso
 */
public class SequenceWrapper implements Sequence {

    private static Logger log = Logger.getLogger(SequenceWrapper.class);
    private static boolean cacheSequences = true;
    private static int tileSize = 1000000;

    private Sequence sequence;
    private ObjectCache<String, SequenceTile> sequenceCache = new ObjectCache<String, SequenceTile>(50);


    public SequenceWrapper(Sequence sequence) {
        this.sequence = sequence;
    }

    public byte getBase(String chr, int position) {
        if (cacheSequences) {
            int tileNo = position / tileSize;

            // Get first chunk
            SequenceTile tile = getSequenceTile(chr, tileNo);
            int offset = position - tile.getStart();
            byte[] bytes = tile.bytes;
            if (offset > 0 && offset < bytes.length) {
                return bytes[offset];
            } else {
                return 0;
            }

        } else {
            // TODO -- implement or disable
            return sequence.getBase(chr, position);
        }
    }

    @Override
    public List<String> getChromosomeNames() {
        return sequence.getChromosomeNames();
    }

    @Override
    public int getChromosomeLength(String chrname) {
        return sequence.getChromosomeLength(chrname);
    }

    /**
     * Return the reference dna sequence for the exact interval specified.
     *
     * @param chr
     * @param start
     * @param end
     * @return
     */
    public byte[] getSequence(String chr, int start, int end) {
        if (cacheSequences) {
            byte[] seqbytes = new byte[end - start];
            int startTile = start / tileSize;
            int endTile = end / tileSize;

            // Get first chunk
            SequenceTile tile = getSequenceTile(chr, startTile);
            if (tile == null) {
                return null;
            }

            byte[] tileBytes = tile.getBytes();
            if (tileBytes == null) {
                return null;
            }

            int fromOffset = start - tile.getStart();
            int toOffset = 0;

            // A negative offset means the requested start is < the the first tile start.  This situation can arise at the
            // left end of chromosomes.  In this case we want to copy the first tile to some offset location in the
            // destination sequence array.
            if (fromOffset < 0) {
                toOffset = -fromOffset;
                fromOffset = 0;
            }

            // # of bytes to copy.  Note that only one of fromOffset or toOffset is non-zero.
            int nBytes = Math.min(tileBytes.length - Math.abs(fromOffset), seqbytes.length - Math.abs(toOffset));

            // Copy first chunk
            System.arraycopy(tileBytes, fromOffset, seqbytes, toOffset, nBytes);

            // If multiple chunks ...
            for (int t = startTile + 1; t <= endTile; t++) {
                tile = getSequenceTile(chr, t);
                if (tile == null) {
                    break;
                }

                int nNext = Math.min(seqbytes.length - nBytes, tile.getSize());

                System.arraycopy(tile.getBytes(), 0, seqbytes, nBytes, nNext);
                nBytes += nNext;
            }

            return seqbytes;
        } else {
            return sequence.getSequence(chr, start, end);
        }
    }


    private SequenceTile getSequenceTile(String chr, int tileNo) {
        String key = getKey(chr, tileNo);
        SequenceTile tile = sequenceCache.get(key);

        if (tile == null) {
            int start = tileNo * tileSize;
            int end = start + tileSize; // <=  UCSC coordinate conventions (end base not inclusive)

            if (end <= start) {
                return null;
            }

            byte[] seq = sequence.getSequence(chr, start, end);
            tile = new SequenceTile(start, seq);
            sequenceCache.put(key, tile);
        }

        return tile;
    }

    /**
     * Generate unique key to be used to store/retrieve
     * tiles. We combined the chr and tileNo, with a
     * delimiter in between to ensure that
     * chr1 12 doesn't clash with chr11 2
     * @param chr
     * @param tileNo
     * @return
     */
    static String getKey(String chr, int tileNo) {
        return chr + "/" + tileNo;
    }

    /**
     * This accessor provided to support unit tests.
     *
     * @param aChunkSize
     */
    static void setTileSize(int aChunkSize) {
        tileSize = aChunkSize;
    }

    /**
     * Accessor to support unit tests.
     *
     * @param aCacheSequences
     */
    static void setCacheSequences(boolean aCacheSequences) {
        cacheSequences = aCacheSequences;
    }

    public void clearCache() {
        sequenceCache.clear();
    }

    static class SequenceTile {

        private int start;
        private byte[] bytes;

        SequenceTile(int start, byte[] bytes) {
            this.start = start;
            this.bytes = bytes;
        }

        public int getStart() {
            return start;
        }

        public int getSize() {
            return bytes == null ? 0 : bytes.length;
        }

        public byte[] getBytes() {
            return bytes;
        }
    }

    /**
     * Translates sequence URLs that might be cached on client machines.  This method should be retired eventually,
     * as caches expire.
     * <p/>
     * Also modifies URLs to Broad hosted sequences that will use byte range requests if byte-range requests are
     * disabled.  This hack is neccessary for the Partners network, which does not forward the byte-range header.
     * <p/>
     * Older "sequence servlet" request URLs
     * http://www.broad.mit.edu/igv/SequenceServlet/
     * http://www.broadinstitute.org/igv/sequence
     * <p/>
     * Direct URLS  (uses byte range requests)
     * http://www.broadinstitute.org/igvdata/annotations/seq/
     * http://igvdata.broadinstitute.org/genomes/seq
     *
     * @param url
     * @return
     */

    private static Hashtable<String, String> sequenceUrlCache = new Hashtable();


    /**
     * Some rather ugly code to maintain backward compatibility.  Does 2 things
     * (1) domain swap  (mit -> broadinstitute)
     * (2) removes references to SequenceServlet, there are 2 forms
     * <p/>
     * This method can be removed when its verified that references to the MIT domain and sequence servlet have
     * been removed from all genomes.
     *
     * @param url
     * @return
     */
    public static String checkSequenceURL(String url) {

        String key = url;
        String convertedURL = sequenceUrlCache.get(key);
        if (convertedURL == null) {
            convertedURL = url;

            // Legacy URLs -- this code can be removed when all .genome files are updated.
            convertedURL = convertedURL.replace("broad.mit.edu", "broadinstitute.org");

            // Replace all references to the old SequenceServlet with direct references to a sequence directory.
            convertedURL = convertedURL.replace(
                    "http://www.broadinstitute.org/igv/SequenceServlet",
                    "http://igvdata.broadinstitute.org/genomes/seq");

            convertedURL = convertedURL.replace(
                    "http://www.broadinstitute.org/igv/sequence",
                    "http://igvdata.broadinstitute.org/genomes/seq");


            if (!url.equals(convertedURL)) {
                log.info("Converting sequence URL: " + url + " -> " + convertedURL);
            }

            sequenceUrlCache.put(key, convertedURL);
        }

        return convertedURL;

    }


}
