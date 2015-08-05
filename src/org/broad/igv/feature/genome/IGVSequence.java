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

package org.broad.igv.feature.genome;

import htsjdk.samtools.seekablestream.SeekableStream;
import org.apache.log4j.Logger;
import org.broad.igv.feature.Cytoband;
import org.broad.igv.util.stream.IGVSeekableStreamFactory;

import java.io.IOException;
import java.util.*;

/**
 * Represents a sequence database composed of plain text files with no white space, one per chromosome, in a directory.
 * This is the original IGV "sequence" format, replaced in favor of indexed fasta files.
 *
 * @author jrobinso
 * @Date 8/8/11
 */

public class IGVSequence implements Sequence {

    private static Logger log = Logger.getLogger(IGVSequence.class);

    private String dirPath;
    private Map<String, String> chrFileNameCache = new HashMap();
    private HashMap<String, Integer> chromosomeLengths;
    private List<String> chromosomeNames;

    public IGVSequence(String dirPath) {
        if (!dirPath.endsWith("/")) {
            dirPath = dirPath + "/";
        }
        this.dirPath = dirPath;
    }

    public byte[] getSequence(String chr, int start, int end) {

        String fn = getChrFileName(chr);
        String seqFile = dirPath + fn;

        SeekableStream is = null;
        try {


            is = IGVSeekableStreamFactory.getInstance().getStreamFor(seqFile);

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


    @Override
    public byte getBase(String chr, int position) {
        throw new RuntimeException("getBase() is not implemented for class " + FastaIndexedSequence.class.getName());
    }

    @Override
    public List<String> getChromosomeNames() {
        return chromosomeNames;
    }

    @Override
    public int getChromosomeLength(String chrname) {
        return chromosomeLengths.get(chrname);
    }

    /**
     * Get a "legal" chromosome file name from the chr name.  This method supports "old" style .genome
     * files.
     *
     * @param chr
     * @return
     */
    private String getChrFileName(String chr) {
        String chrFN = chrFileNameCache.get(chr);
        if (chrFN == null) {
            chrFN = chr;
            for (Map.Entry<String, String> entry : illegalChar.entrySet()) {
                chrFN = chrFN.replaceAll(entry.getValue(), entry.getKey());
            }
            chrFN += ".txt";
            chrFileNameCache.put(chr, chrFN);
        }
        return chrFN;
    }

    /**
     * Generate chromosomes from the list of cytobands.  This method is provided for backward compatibility for
     * pre V2.1 genomes.
     *
     * @param chrCytoMap
     * @param chromosomesAreOrdered
     */
    public void generateChromosomes(LinkedHashMap<String, List<Cytoband>> chrCytoMap, boolean chromosomesAreOrdered) {

        chromosomeLengths = new HashMap<String, Integer>();
        chromosomeNames = new ArrayList<String>(chrCytoMap.size());
        for (Map.Entry<String, List<Cytoband>> entry : chrCytoMap.entrySet()) {
            String chr = entry.getKey();
            chromosomeNames.add(chr);

            List<Cytoband> cytobands = entry.getValue();
            int length = cytobands.get(cytobands.size() - 1).getEnd();
            chromosomeLengths.put(chr, length);
        }

        if (!chromosomesAreOrdered) {
            Collections.sort(chromosomeNames, ChromosomeNameComparator.get());
        }
    }


    static Map<String, String> illegalChar = new HashMap();

    static {
        illegalChar.put("_qm_", "\\?");
        illegalChar.put("_fbr_", "\\[");
        illegalChar.put("_rbr_", "]");
        illegalChar.put("_fsl_", "/");
        illegalChar.put("_bsl_", "\\\\");
        illegalChar.put("_eq_", "=");
        illegalChar.put("_pl_", "\\+");
        illegalChar.put("_lt_", "<");
        illegalChar.put("_gt_", ">");
        illegalChar.put("_co_", ":");
        illegalChar.put("_sc_", ";");
        illegalChar.put("_dq_", "\"");
        illegalChar.put("_sq_", "'");
        illegalChar.put("_st_", "\\*");
        illegalChar.put("_pp_", "\\|");
    }

}
