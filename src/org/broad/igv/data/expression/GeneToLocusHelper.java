/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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

package org.broad.igv.data.expression;

import org.broad.igv.feature.*;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.tribble.Feature;
import org.broad.tribble.readers.AsciiLineReader;

import java.io.IOException;
import java.util.*;

/**
 * @author jrobinso
 * @date Oct 9, 2010
 */
public class GeneToLocusHelper {

    static String LOCUS_START_DELIMITER = "|@";
    static String LOCUS_END_DELIMITER = "|";

    Map<String, IGVFeature> probeLocusMap;

    public GeneToLocusHelper(String probeResource, Genome genome) throws IOException {

        if (probeResource != null && probeResource.trim().length() > 0) {
            BEDFileParser parser = new BEDFileParser(genome);
            ResourceLocator rl = new ResourceLocator(probeResource);
            AsciiLineReader reader = ParsingUtils.openAsciiReader(rl);
            List<Feature> features = parser.loadFeatures(reader);
            reader.close();
            probeLocusMap = new HashMap(features.size() * 2);
            for (Feature f : features) {
                probeLocusMap.put(((IGVFeature) f).getName(), (IGVFeature) f);
            }
        }

    }

    public List<Locus> getLoci(String probeId, String description) {

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

        String[] genes = ProbeToLocusMap.getInstance().getLociForProbe(probeId);
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
        Locus locus = new Locus(geneOrLocusString);
        if (locus.isValid()) {
            return locus;
        } else {

            if (probeLocusMap != null) {
                Feature f = probeLocusMap.get(geneOrLocusString);
                if (f == null) {
                    return null;
                } else {
                    return new Locus(f.getChr(), f.getStart(), f.getEnd());
                }
            }

            // Maybe its a gene or feature
            Feature gene = FeatureDB.getFeature(geneOrLocusString);
            if (gene != null) {
                return new Locus(gene.getChr(), gene.getStart(), gene.getEnd());
            }
        }
        return null;
    }
}
