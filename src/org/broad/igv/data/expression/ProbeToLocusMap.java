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
import org.broad.igv.exceptions.LoadResourceFromServerException;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.HttpUtils;
import htsjdk.tribble.readers.AsciiLineReader;

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

    enum Platform {Affymetrix, Agilient, Illumina, Methylation, Mirna, unknown}

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

    private static ProbeToLocusMap instance;

    private Map<String, Map<String, String[]>> probeMaps = new HashMap();

    private MappingUrlCache mappingUrlCache = new MappingUrlCache();


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


    public void loadMapping(String urlString, Map<String, String[]> map) {
        AsciiLineReader bufReader = null;
        InputStream is = null;
        try {
            if (HttpUtils.isRemoteURL(urlString)) {
                URL url = new URL(urlString);
                is = HttpUtils.getInstance().openConnectionStream(url);
            } else {
                is = new FileInputStream(urlString);
            }
            if (urlString.endsWith("gz")) {
                is = new GZIPInputStream(is);
            }
            bufReader = new AsciiLineReader(is);
            loadMapping(bufReader, map);

        } catch (Exception e) {
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
        while ((line = bufReader.readLine()) != null) {
            String[] result = Globals.tabPattern.split(line, -1);
            int nTokens = result.length;
            if (nTokens != 2) {
                continue;
            }
            // check if more than one gene symbol found
            String[] genes = result[1].split("///");

            map.put(result[0].trim(), genes);
        }
    }

    /**
     * Return the URL to the mapping file for the give genomeId and platform.
     *
     * @param genomeId
     * @param platform
     * @return
     */
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

    public String[] getLociForProbe(String probeId, String genomeId) {

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
