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

package org.broad.igv.feature.genome;

import org.apache.log4j.Logger;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.ChromosomeImpl;
import org.broad.igv.feature.Cytoband;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.stream.IGVSeekableStreamFactory;
import org.broad.tribble.util.SeekableStream;

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
        for (Map.Entry<String, List<Cytoband>> entry : chrCytoMap.entrySet()) {
            String chr = entry.getKey();
            List<Cytoband> cytobands = entry.getValue();

            int length = cytobands.get(cytobands.size() - 1).getEnd();
            final ChromosomeImpl chromosome = new ChromosomeImpl(chr, length);
            chromosomeLengths.put(chr, length);
        }

        this.chromosomeNames = new LinkedList<String>(chromosomeLengths.keySet());
        if (!chromosomesAreOrdered) {
            Collections.sort(chromosomeNames, new ChromosomeComparator());
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
