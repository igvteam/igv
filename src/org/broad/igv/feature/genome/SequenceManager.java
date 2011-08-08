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
import org.broad.igv.ui.IGV;
import org.broad.igv.util.*;
import org.broad.igv.util.stream.IGVSeekableStreamFactory;
import org.broad.tribble.util.SeekableStream;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * @author jrobinso
 */
public class SequenceManager {

    private static Logger log = Logger.getLogger(SequenceManager.class);
    private static boolean cacheSequences = true;
    private static int chunkSize = 30000;
    private static ObjectCache<String, SequenceChunk> sequenceCache = new ObjectCache(50);

    /**
     * Return the sequence in Color Space (SOLID alignment encoding)
     *
     * @param genome
     * @param chr
     * @param start
     * @param end
     * @return
     */
    public static byte[] readCSSequence(String genome, String chr, int start, int end) {
        // We need to know the base just to the left of the start
        int csStart = (start == 0 ? 0 : start - 1);
        byte[] baseSequence = readSequence(genome, chr, csStart, end);
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
    public static byte[] readSequence(String genome, String chr, int start, int end) {

        File rootDirectory = null;

        try {
            GenomeDescriptor descriptor = IGV.getInstance().getGenomeManager().getGenomeDescriptor(genome);
            if (descriptor != null) {
                String location = descriptor.getSequenceLocation();
                if (location != null) {
                    location = convertSequenceURL(location, genome);
                    return getSequence(descriptor, chr, start, end, location);
                }
            }
        } catch (Exception e) {
            log.error(e.getMessage(), e);

        }
        return new byte[0];
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
    private static byte[] getSequence(GenomeDescriptor genome, String chr, int start, int end, String location) {
        if (cacheSequences) {
            byte[] sequence = new byte[end - start];
            int startTile = start / chunkSize;
            int endTile = end / chunkSize;

            // Get first chunk
            SequenceChunk chunk = getSequenceChunk(genome, chr, startTile, location);
            int offset = start - chunk.getStart();
            byte[] seqBytes = chunk.getBytes();
            if (seqBytes == null) {
                return null;
            }

            // # of bytes to return, minimum of requested sequence lenth or bytes available
            int nBytes = Math.min(seqBytes.length - offset, sequence.length);

            // Copy first chunk
            System.arraycopy(chunk.getBytes(), offset, sequence, 0, nBytes);

            // If multiple chunks ...
            for (int tile = startTile + 1; tile <= endTile; tile++) {
                chunk = getSequenceChunk(genome, chr, tile, location);

                int nNext = Math.min(sequence.length - nBytes, chunk.getSize());

                System.arraycopy(chunk.getBytes(), 0, sequence, nBytes, nNext);
                nBytes += nNext;
            }

            return sequence;
        } else {
            return readSequence(genome, chr, start, end, location);
        }
    }


    private static SequenceChunk getSequenceChunk(GenomeDescriptor genome, String chr, int tileNo, String location) {
        String key = getKey(genome.getId(), chr, tileNo);
        SequenceChunk chunk = sequenceCache.get(key);

        if (chunk == null) {
            int start = tileNo * chunkSize;
            int end = start + chunkSize; // <=  UCSC coordinate conventions (end base not inclusive)

            byte[] seq = readSequence(genome, chr, start, end, location);
            chunk = new SequenceChunk(start, seq);
            sequenceCache.put(key, chunk);
        }

        return chunk;
    }

    /**
     * Read and return the genomic sequence for the specified interval as a byte array.
     *
     * @param genome
     * @param chr
     * @param start
     * @param end
     * @param path
     * @return
     */
    private static byte[] readSequence(GenomeDescriptor genome, String chr, int start, int end, String path) {

        if (path.startsWith("http://www.broadinstitute.org/igv/SequenceServlet") ||
                path.startsWith("http://www.broadinstitute.org/igv/sequence")) {
            return SequenceServletWrapper.readBytes(path, chr, start, end);
        } else {

            if (!path.endsWith("/")) {
                path = path + "/";
            }

            String fn = chr + ".txt";
            if (genome.isChrNamesAltered()) {
                fn = getChrFileName(fn);
            }

            String seqFile = path + fn;

            return readSequence(chr, start, end, seqFile);


        }
    }


    /**
     * Note: This method is "package" scope to permit unit testing.
     *
     * @param chr
     * @param start
     * @param end
     * @param seqFile
     * @return
     */
    static byte[] readSequence(String chr, int start, int end, String seqFile) {

        SeekableStream is = null;
        try {


            is = IGVSeekableStreamFactory.getStreamFor(seqFile);

            byte[] bytes = new byte[end - start];
            is.seek(start);
            is.read(bytes);
            return bytes;


        } catch (Exception ex) {
            log.error("Error reading genome sequence from: " + seqFile, ex);
            return null;
        } finally {
            if (is != null) {
                try {
                    is.close();
                } catch (IOException ex) {
                    log.error("Error closing sequence file.", ex);
                }
            }
        }
    }


    private static Map<String, String> chrFileNameCache = new HashMap();

    private static String getChrFileName(String fn) {
        String chrFN = chrFileNameCache.get(fn);
        if (chrFN == null) {
            chrFN = FileUtils.legalFileName(fn);
            chrFileNameCache.put(fn, chrFN);
        }
        return chrFN;
    }

    static String getKey(String genome, String chr, int tileNo) {
        return genome + chr + tileNo;
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

    public static void clearCache() {
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


    public static String convertSequenceURL(String url, String genome) {

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
