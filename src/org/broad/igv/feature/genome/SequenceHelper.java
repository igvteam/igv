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
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.feature.genome;

import org.apache.log4j.Logger;
import org.broad.igv.util.*;

import java.util.*;

/**
 * @author jrobinso
 */
public class SequenceHelper {

    private static Logger log = Logger.getLogger(SequenceHelper.class);
    private static boolean cacheSequences = true;
    private static int chunkSize = 30000;


    //Genome genome;
    Sequence sequence;
    private ObjectCache<String, SequenceChunk> sequenceCache = new ObjectCache(50);


    public SequenceHelper(String seqpath) {

        if (seqpath == null) {

        } else {
            seqpath = convertSequenceURL(seqpath);
            if (seqpath.startsWith("http://www.broadinstitute.org/igv/SequenceServlet") ||
                    seqpath.startsWith("http://www.broadinstitute.org/igv/sequence")) {
                sequence = new ServletSequence(seqpath);
            } else {
                sequence = new IGVSequence(seqpath);
            }
        }
    }

    public SequenceHelper(Sequence sequence) {
        this.sequence = sequence;
    }

    /**
     * Return the sequence in Color Space (SOLID alignment encoding)
     *
     * @param genome
     * @param chr
     * @param start
     * @param end
     * @return
     */
    public byte[] readCSSequence(String genome, String chr, int start, int end) {
        // We need to know the base just to the left of the start
        int csStart = (start == 0 ? 0 : start - 1);
        byte[] baseSequence = sequence.readSequence(chr, csStart, end);
        if (baseSequence == null || baseSequence.length == 0) {
            return baseSequence;
        }

        byte[] csSequence = new byte[end - start];
        int i = 0;
        int c1 = start == 0 ? 0 : baseToCS(baseSequence[i++]);
        for (; i < baseSequence.length; i++) {
            int c2 = baseToCS(baseSequence[i]);
            csSequence[i] = (byte) (c1 ^ c2);
        }
        return csSequence;

    }

    private static int baseToCS(byte base) {
        switch (base) {
            case 'A':
            case 'a':
                return 0;
            case 'C':
            case 'c':
                return 1;
            case 'T':
            case 't':
                return 2;
            case 'G':
            case 'g':
                return 3;
        }
        return -1;
    }


    /**
     * Return the reference dna sequence for the exact interval specified.
     *
     * @param genome
     * @param chr
     * @param start
     * @param end
     * @return
     */
    public byte[] readSequence(String genome, String chr, int start, int end) {

        if (sequence != null) {
            return getSequence(chr, start, end);

        }
        return new byte[0];
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
            int startTile = start / chunkSize;
            int endTile = end / chunkSize;

            // Get first chunk
            SequenceChunk chunk = getSequenceChunk(chr, startTile);
            int offset = start - chunk.getStart();
            byte[] seqBytes = chunk.getBytes();
            if (seqBytes == null) {
                return null;
            }

            // # of bytes to return, minimum of requested sequence lenth or bytes available
            int nBytes = Math.min(seqBytes.length - offset, seqbytes.length);

            // Copy first chunk
            System.arraycopy(chunk.getBytes(), offset, seqbytes, 0, nBytes);

            // If multiple chunks ...
            for (int tile = startTile + 1; tile <= endTile; tile++) {
                chunk = getSequenceChunk(chr, tile);

                int nNext = Math.min(seqbytes.length - nBytes, chunk.getSize());

                System.arraycopy(chunk.getBytes(), 0, seqbytes, nBytes, nNext);
                nBytes += nNext;
            }

            return seqbytes;
        } else {
            return sequence.readSequence(chr, start, end);
        }
    }


    private SequenceChunk getSequenceChunk(String chr, int tileNo) {
        String key = getKey(chr, tileNo);
        SequenceChunk chunk = sequenceCache.get(key);

        if (chunk == null) {
            int start = tileNo * chunkSize;
            int end = start + chunkSize; // <=  UCSC coordinate conventions (end base not inclusive)

            byte[] seq = sequence.readSequence(chr, start, end);
            chunk = new SequenceChunk(start, seq);
            sequenceCache.put(key, chunk);
        }

        return chunk;
    }


    static String getKey(String chr, int tileNo) {
        return chr + tileNo;
    }

    /**
     * This accessor provided to support unit tests.
     *
     * @param aChunkSize
     */
    static void setChunkSize(int aChunkSize) {
        chunkSize = aChunkSize;
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

    static class SequenceChunk {

        private int start;
        private byte[] bytes;

        SequenceChunk(int start, byte[] bytes) {
            this.start = start;
            this.bytes = bytes;
        }

        public int getStart() {
            return start;
        }

        public int getSize() {
            return bytes.length;
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


    public static String convertSequenceURL(String url) {

        final boolean useByteRange = IGVHttpClientUtils.useByteRange();
        String key = url + useByteRange;
        String convertedURL = sequenceUrlCache.get(key);
        if (convertedURL == null) {
            convertedURL = url;
            // Legacy URLs -- this code can be removed when all .genome files are updated.
            convertedURL = convertedURL.replace(
                    "broad.mit.edu",
                    "broadinstitute.org");


            convertedURL = convertedURL.replace(
                    "http://www.broadinstitute.org/igv/SequenceServlet",
                    "http://igvdata.broadinstitute.org/genomes/seq");


            if (!useByteRange) {
                // Translate our URLS to use the servlet
                convertedURL = convertedURL.replace(
                        "http://www.broadinstitute.org/igvdata/annotations/seq",
                        "http://www.broadinstitute.org/igv/sequence");

                convertedURL = convertedURL.replace(
                        "http://igvdata.broadinstitute.org/genomes/seq",
                        "http://www.broadinstitute.org/igv/sequence");

                convertedURL = convertedURL.replace(
                        "http://igv.broadinstitute.org/genomes/seq",
                        "http://www.broadinstitute.org/igv/sequence");
            }


            if (!url.equals(convertedURL)) {
                log.info("Converting sequence URL: " + url + " -> " + convertedURL);
            }

            sequenceUrlCache.put(key, convertedURL);
        }

        return convertedURL;

    }


}
