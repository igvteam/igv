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

package org.broad.igv.data.expression;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.FeatureDB;
import org.broad.igv.feature.Locus;
import org.broad.igv.session.Session;
import org.broad.igv.ui.IGV;
import org.broad.igv.util.ParsingUtils;
import htsjdk.tribble.Feature;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.*;

/**
 * @author jrobinso
 * @date Oct 9, 2010
 */
public class GeneToLocusHelper {

    private static Logger log = Logger.getLogger(GeneToLocusHelper.class);

    static String LOCUS_START_DELIMITER = "|@";
    static String LOCUS_END_DELIMITER = "|";

    Map<String, List<Locus>> probeLocusMap;

    /**
     * Create an instance and loaded the supplied probe mapping file.
     *
     * @param probeResource - file path or URL to a bed file containing the probe mappings.  CAN BE NULL.
     * @throws IOException
     */
    public GeneToLocusHelper(String probeResource) throws IOException {

        if (probeResource != null && probeResource.trim().length() > 0) {
            loadProbeMap(probeResource);
        }

        // If a probe file is supplied,  see if there is a user default.  Only do this if the custom file option
        // is set and a probe mapping file has been supplied.
        boolean use_probe_mf;
        String userMappingFile;
        if (!Globals.isHeadless()){
            Session session = IGV.getInstance().getSession();
            use_probe_mf = session.getPreferenceAsBoolean(PreferenceManager.USE_PROBE_MAPPING_FILE);
            userMappingFile = session.getPreference(PreferenceManager.PROBE_MAPPING_FILE);
        }else{
            use_probe_mf = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.USE_PROBE_MAPPING_FILE);
            userMappingFile = PreferenceManager.getInstance().get(PreferenceManager.PROBE_MAPPING_FILE);
        }
        
        if (use_probe_mf) {
            if (userMappingFile != null && userMappingFile.trim().length() > 0) {
                loadProbeMap(userMappingFile);
            }
        }

    }

    /**
     * Return a list of loci mapping to the given probe.
     *
     * @param probeId     - the probe (or gene) ID
     * @param description - optional description field.  Some formats allow encoding of lcous in description fiel
     * @param genomeId    - the genome
     * @return
     */
    public List<Locus> getLoci(String probeId, String description, String genomeId) {


        if (probeLocusMap != null) {
            return probeLocusMap.get(probeId);
        }

        // Search for locus in description string.  This relies on the special
        // IGV convention for specifying loci  (e.g  |@chrX:1000-2000|
        if ((description != null) && (description.length() > 3)) {
            String[] locusStrings = getExplicitLocusStrings(description);
            if (locusStrings != null) {

                List<Locus> loci = new ArrayList(locusStrings.length);
                for (String ls : locusStrings) {
                    ls = ls.trim();
                    Locus locus = getLocus(ls);
                    if ((locus != null) && locus.isValid()) {
                        loci.add(locus);
                    }
                }
                return loci;
            }
        }

        // Search for locus from the probe name itself.
        Locus locus = getLocus(probeId);
        if ((locus != null) && locus.isValid()) {
            return Arrays.asList(locus);
        }


        // See if the probes can be mapped to genes
        String[] genes = ProbeToLocusMap.getInstance().getLociForProbe(probeId, genomeId);
        if (genes != null) {
            List<Locus> loci = new ArrayList(genes.length);
            for (String g : genes) {
                locus = getLocus(g);
                if (locus != null) {
                    loci.add(locus);
                }
            }
            return loci;
        }
        return null;
    }


    /**
     * Search for a locus explicitly specified in the description field.
     * A locus can be specified either directrly, as a UCSC style locus string
     * (e.g. chrX:1000-2000), or indirectly as a HUGO gene symbol (e.g. egfr).
     * The locus string is distinguished by the  delimiters |@ and |.
     */
    private String[] getExplicitLocusStrings(String description) {

        // Search for locus in description string
        int startIndex = description.indexOf(LOCUS_START_DELIMITER);
        if (startIndex < 0) {
            return null;
        } else {
            startIndex += 2;
        }

        int endIndex = description.indexOf(LOCUS_END_DELIMITER, startIndex + 1);
        if (endIndex < 0) {

            // Assume the locus extends to the end of the string
            endIndex = description.length();
        }

        if (endIndex > startIndex + 3) {
            String locusString = description.substring(startIndex, endIndex);
            if (locusString.contains(",")) {
                return locusString.split(",");
            } else {
                return new String[]{locusString};
            }
        }
        return null;

    }

    /**
     * Return a locus from the gene (e.g. EGFR) or locus (e.g. chr1:1-100) string
     */
    Locus getLocus(String geneOrLocusString) {
        Locus locus = Locus.fromString(geneOrLocusString);
        if (locus != null && locus.isValid()) {
            return locus;
        } else {
            // Maybe its a gene or feature
            Feature gene = FeatureDB.getFeature(geneOrLocusString);
            if (gene != null) {
                return new Locus(gene.getChr(), gene.getStart(), gene.getEnd());
            }
        }
        return null;
    }


    private void loadProbeMap(String probeResource) {

        BufferedReader reader = null;
        try {
            probeLocusMap = new HashMap<String, List<Locus>>(50000);
            int maxErrors = 10;
            int errorCount = 0;
            reader = ParsingUtils.openBufferedReader(probeResource);
            String nextLine;
            while ((nextLine = reader.readLine()) != null) {
                if (nextLine.startsWith("#") || nextLine.startsWith("browser") || nextLine.startsWith("track")) {
                    continue;
                }
                String[] tokens = Globals.tabPattern.split(nextLine, -1);
                if (tokens.length > 3) {
                    try {
                        String chr = tokens[0];
                        int start = (int) Double.parseDouble(tokens[1]);
                        int end = (int) Double.parseDouble(tokens[2]);
                        String probe = tokens[3];
                        Locus locus = new Locus(chr, start, end);
                        probeLocusMap.put(probe, Arrays.asList(locus));
                    } catch (NumberFormatException e) {
                        log.info("Skipping line: " + nextLine);
                        errorCount++;
                        if (errorCount > maxErrors) {
                            probeLocusMap = null;
                            log.info("Too many errors.  Aborting.");
                        }
                    }
                }
            }
        } catch (IOException e) {

        } finally {
            if (reader != null) {
                try {
                    reader.close();
                } catch (IOException e) {

                }
            }
        }
    }

}
