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
package org.broad.igv.data.expression;


import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.exceptions.LoadResourceFromServerException;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.IGVHttpClientUtils;
import org.broad.tribble.readers.AsciiLineReader;
import org.broad.igv.util.ParsingUtils;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Map;
import java.util.zip.GZIPInputStream;

/**
 * @author: Marc-Danie Nazaire
 * private static String affyMappingURL =
 * SERVER_URL + "/igv/resources/probes/affy_probe_gene_mapping.txt.gz";
 * private static String agilentMappingURL =
 * SERVER_URL + "/igv/resources/probes/agilent_probe_gene_mapping.txt.gz";
 * private static String illuminaMappingURL =
 * SERVER_URL + "/igv/resources/probes/illumina_probe_gene_mapping.txt.gz";
 */
public class ProbeToLocusMap {
    public static final String SERVER_URL = "http://www.broadinstitute.org";

    private static Logger log = Logger.getLogger(ProbeToLocusMap.class);
    private static String affyGenesMappingURL =
            SERVER_URL + "/igvdata/probes/affy/affy_probe_gene_mapping.txt.gz";
    private static String affyHumanMappingURL =
            SERVER_URL + "/igvdata/probes/affy/affy_human_mappings.txt.gz";
    private static String affyMouseMappingURL =
            SERVER_URL + "/igvdata/probes/affy/affy_mouse_mappings.txt.gz";
    private static String affyOtherMappingURL =
            SERVER_URL + "/igvdata/probes/affy/affy_other_mappings.txt.gz";
    private static String agilentGenesMappingURL =
            SERVER_URL + "/igvdata/probes/agilent/agilent_probe_gene_mapping.txt.gz";
    private static String agilentHumanMappingURL =
            SERVER_URL + "/igvdata/probes/agilent/agilent_human_mappings.txt.gz";
    private static String agilentMouseMappingURL =
            SERVER_URL + "/igvdata/probes/agilent/agilent_mouse_mappings.txt.gz";
    private static String agilentOtherMappingURL =
            SERVER_URL + "/igvdata/probes/agilent/agilent_other_mappings.txt.gz";
    private static String illuminaMappingURL =
            SERVER_URL + "/igvdata/probes/illumina/illumina_allMappings.txt.gz";
    private static String illuminaGenesMappingURL =
            SERVER_URL + "/igvdata/probes/illumina/illumina_probe_gene_mapping.txt.gz";
    private static String methylationGeneMappingURL =
            SERVER_URL + "/igvdata/probes/meth/methylation_pobeToGene.tab.gz";
    private static String methylationLociMappingURL =
            SERVER_URL + "/igvdata/probes/meth/methylation_probeToLoci.mappings.txt.gz";
    private static Map<String, String[]> rnaiMap;
    private static ProbeToLocusMap instance;
    Map<String, Map<String, String[]>> probeMaps = new HashMap();
    MappingUrlCache mappingUrlCache = new MappingUrlCache();

    enum Platform {

        Affymetrix, Agilient, Illumina, Methylation, Mirna, unknown
    }

    ;

    /**
     * Method description
     *
     * @return
     */
    public static synchronized ProbeToLocusMap getInstance() {
        if (instance == null) {
            instance = new ProbeToLocusMap();
        }

        return instance;
    }

    public void clearProbeMappings() {
        if (probeMaps != null) {
            probeMaps.clear();
        }
        if (mappingUrlCache != null) {
            mappingUrlCache.clear();
        }
    }

    /**
     * Method description
     *
     * @param urlString
     * @param map
     * @throws Exception
     */
    public void loadMapping(String urlString, Map<String, String[]> map) {
        AsciiLineReader bufReader = null;
        InputStream is = null;
        try {
            if (IGVHttpClientUtils.isURL(urlString)) {
                URL url = new URL(urlString);
                is = IGVHttpClientUtils.openConnectionStream(url);
            } else {
                is = new FileInputStream(urlString);
            }
            if (urlString.endsWith("gz")) {
                is = new GZIPInputStream(is);
            }
            bufReader = new AsciiLineReader(is);
            loadMapping(bufReader, map);

        }
        catch (Exception e) {
            log.error("Error loading probe mapping", e);

            throw new LoadResourceFromServerException(e.getMessage(), urlString, e.getClass().getSimpleName());
        } finally {
            if (is != null) {
                try {
                    is.close();
                } catch (IOException e) {
                    log.error("Error closing probe mapping stream", e);
                }
            }
        }
    }

    public void loadMapping(AsciiLineReader bufReader, Map<String, String[]> map) throws IOException {
        String line;
        String[] result = new String[100];
        while ((line = bufReader.readLine()) != null) {
            int nTokens = ParsingUtils.split(line, result, '\t');
            if (nTokens != 2) {
                continue;
            }
            // check if more than one gene symbol found
            String[] genes = result[1].split("///");

            map.put(result[0].trim(), genes);
        }
    }

    /**
     * Method description
     *
     * @return
     */
    public Map<String, String[]> getRNAiProbeMap() {
        return rnaiMap;
    }

    public String getMappingURL(String genomeId, Platform platform) {

        if (platform == Platform.unknown) {
            return null;
        }

        String mappingUrl = mappingUrlCache.getMappingUrl(genomeId, platform);
        if (mappingUrl == null) {

            boolean mapToGenes = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.PROBE_MAPPING_KEY);
            if (!mapToGenes) {
                boolean hasLociMapping = checkForLociMapping(platform, genomeId);
                if (!hasLociMapping) {
                    MessageUtils.showMessage(
                            "<html>" + platform.toString() + " probe locations are not available for the selected genome " +
                                    " (" + genomeId + "). <br>Expression data will be mapped to gene locations.");
                    mapToGenes = true;
                }
            }

            if (platform == Platform.Affymetrix) {
                if (mapToGenes) {
                    mappingUrl = affyGenesMappingURL;
                } else if (genomeId.startsWith("hg")) {
                    mappingUrl = affyHumanMappingURL;

                } else if (genomeId.startsWith("mm")) {
                    mappingUrl = affyMouseMappingURL;

                } else {
                    mappingUrl = affyOtherMappingURL;
                }
            } else if (platform == Platform.Agilient) {
                if (mapToGenes) {
                    mappingUrl = agilentGenesMappingURL;
                } else if (genomeId.startsWith("hg")) {
                    mappingUrl = agilentHumanMappingURL;

                } else if (genomeId.startsWith("mm")) {
                    mappingUrl = agilentMouseMappingURL;
                } else {
                    mappingUrl = agilentOtherMappingURL;
                }
            } else if (platform == Platform.Illumina) {
                if (mapToGenes) {
                    mappingUrl = illuminaGenesMappingURL;
                } else {
                    mappingUrl = illuminaMappingURL;
                }
            } else if (platform == Platform.Methylation) {
                if (mapToGenes) {
                    mappingUrl = methylationGeneMappingURL;
                } else {
                    mappingUrl = methylationLociMappingURL;
                }
            } else {
                return null;
            }
            mappingUrlCache.put(genomeId, platform, mappingUrl);
        }
        return mappingUrl;

    }

    public static Platform getPlatform(String probeId) {
        if (probeId.endsWith("_at") || probeId.endsWith("_st")) {
            return Platform.Affymetrix;
        } else if (probeId.startsWith("A_")) {
            return Platform.Agilient;
        } else if (probeId.startsWith("ILMN_") || probeId.startsWith("GI_") || probeId.startsWith(
                "NM_") || probeId.startsWith("XM_")) {
            return Platform.Illumina;
        } else if (probeId.startsWith("cg")) {
            return Platform.Methylation;
        }
        return Platform.unknown;
    }

    //mm9 (Affymetrix), mm5 (Agilent) or mm8 (Illumina)

    private static boolean checkForLociMapping(Platform platform, String genomeId) {

        boolean hasLociMapping = true;
        if (genomeId.startsWith("hg") && !genomeId.equals("hg18")) {
            hasLociMapping = false;
        } else if (genomeId.startsWith("mm")) {
            switch (platform) {
                case Affymetrix:
                    if (!genomeId.equals("mm9")) {
                        hasLociMapping = false;
                    }
                    break;

                case Agilient:
                    if (!genomeId.equals("mm5")) {
                        hasLociMapping = false;
                    }
                    break;

                case Illumina:
                    if (!genomeId.equals("mm8")) {
                        hasLociMapping = false;
                    }
                    break;
            }
        }
        return hasLociMapping;
    }

    public String[] getLociForProbe(String probeId) {

        String genomeId = IGV.getInstance().getGenomeManager().getGenomeId();
        Platform platform = getPlatform(probeId);

        if (platform == Platform.unknown) {
            return null;
        }

        String mappingURL = getMappingURL(genomeId, platform);

        if (mappingURL == null) {
            return null;
        }

        Map<String, String[]> pMap = probeMaps.get(mappingURL);
        if (pMap == null) {
            pMap = new HashMap(500000);
            loadMapping(mappingURL, pMap);
            probeMaps.put(mappingURL, pMap);
        }

        if (log.isDebugEnabled() && pMap == null) {
            log.debug(("Null probeMap for:" + mappingURL));
        }

        return pMap == null ? null : (String[]) pMap.get(probeId);
    }


    // TODO -- this could be generalized to a Generic 2 key cache

    static class MappingUrlCache {

        Map<Platform, Map<String, String>> cache = new Hashtable();

        public void clear() {
            cache.clear();
        }

        public String getMappingUrl(String genomeId, Platform platform) {
            if (cache.containsKey(platform)) {
                Map<String, String> urlMap = cache.get(platform);
                return urlMap.get(genomeId);
            } else {
                return null;
            }
        }

        public void put(String genomeId, Platform platform, String mappingUrl) {
            if (!cache.containsKey(platform)) {
                cache.put(platform, new Hashtable());
            }
            cache.get(platform).put(genomeId, mappingUrl);
        }
    }


}
