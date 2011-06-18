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

package org.broad.igv.track;

import org.apache.log4j.Logger;
import org.broad.igv.feature.*;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.ui.IGV;
import org.broad.tribble.Feature;
import org.broad.igv.sam.reader.FeatureIndex;
import org.broad.igv.sam.reader.SamUtils;
import org.broad.igv.ui.WaitCursorManager;

import java.io.*;
import java.util.*;

/**
 * User: jrobinso
 * Date: Apr 21, 2010
 */
public class IndexedBEDFeatureSource implements FeatureSource {

    static Logger log = Logger.getLogger(IndexedBEDFeatureSource.class);

    File bedFile;
    FeatureIndex featureIndex;
    BEDFileParser parser;
    HashMap<String, String> chrMappings = new HashMap();

    public IndexedBEDFeatureSource(File samFile, File indexFile, Genome genome) {
        this.bedFile = samFile;
        if (indexFile.exists()) {
            featureIndex = new FeatureIndex(indexFile);
        }
        parser = new BEDFileParser(genome);
        initChrMap();
    }

    private void initChrMap() {
        if (featureIndex != null) {
            Genome genome = IGV.getInstance().getGenomeManager().getCurrentGenome();
            if (genome != null) {
                Set<String> seqNames = featureIndex.getIndexedChromosomes();
                if (seqNames != null) {
                    for (String chr : seqNames) {
                        String alias = genome.getChromosomeAlias(chr);
                        chrMappings.put(alias, chr);
                    }

                }
            }
        }
    }

    public Class getFeatureClass() {
        return IGVFeature.class;
    }

    public Iterator<Feature> getFeatures(String chr, int start, int end) {
        return loadFeatures(chr, start, end).iterator();
    }

    public List<LocusScore> getCoverageScores(String chr, int i, int i1, int zoom) {
        return null;
    }

    public int getFeatureWindowSize() {
        return 0;  
    }

    public void setFeatureWindowSize(int size) {
        // ignored
    }

    private FeatureIndex getIndex() {
        if (featureIndex == null) {
            featureIndex = SamUtils.getIndexFor(bedFile.getAbsolutePath());
        }
        return featureIndex;
    }

    public Set<String> getSequenceNames() {
        FeatureIndex idx = getIndex();
        if (idx == null) {
            return null;
        } else {
            return idx.getIndexedChromosomes();
        }

    }


    private List<Feature> loadFeatures(String chr, int start, int end) {

        final List<Feature> features = new ArrayList();
        if (featureIndex == null) {
            featureIndex = SamUtils.getIndexFor(bedFile.getAbsolutePath());
        }

        if (featureIndex == null) {
            throw new UnsupportedOperationException("SAM files must be indexed to support query methods");
        }
        
        String chrAlias = chrMappings.containsKey(chr) ? chrMappings.get(chr) : chr;
        
        if (!featureIndex.containsChromosome(chrAlias)) {
            return features;
        }
        WaitCursorManager.CursorToken ct = WaitCursorManager.showWaitCursor();
        try {


            // If contained == false (include overlaps) we need to adjust the start to
            // ensure we get features that extend into this segment.
            int startAdjustment = featureIndex.getLongestFeature(chrAlias);
            int startTileNumber = Math.max(0, (start - startAdjustment)) / featureIndex.getTileWidth();

            FeatureIndex.TileDef seekPos = featureIndex.getTileDef(chrAlias, startTileNumber);

            if (seekPos != null) {
                FileInputStream is = null;
                try {
                    // Skip to the start of the chromosome (approximate)
                    is = new FileInputStream(bedFile);
                    is.getChannel().position(seekPos.getStartPosition());

                    BufferedReader br = new BufferedReader(new InputStreamReader(is));
                    String nextLine = "";
                    while ((nextLine = br.readLine()) != null) {
                        String[] tokens = nextLine.split("\t");
                        Feature f = parser.parseLine(tokens, tokens.length);
                        if (f.getStart() > end || !f.getChr().equals(chrAlias)) {
                            break;
                        } else if (f.getEnd() < start) {
                            continue;
                        }
                        features.add(f);
                    }


                } catch (IOException ex) {
                    log.error("Error opening sam file", ex);
                }
                finally {
                    if (is != null) {
                        try {
                            is.close();
                        } catch (IOException e) {
                        }
                    }
                }
            }
        } finally {
            WaitCursorManager.removeWaitCursor(ct);
        }
        return features;
    }

    public boolean hasIndex() {
        if (featureIndex == null) {
            featureIndex = SamUtils.getIndexFor(bedFile.getAbsolutePath());
        }
        return featureIndex != null;
    }


}
