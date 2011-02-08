/**
 * Copyright (c) 2010-2011 by Fred Hutchinson Cancer Research Center.  All Rights Reserved.

 * This software is licensed under the terms of the GNU Lesser General
 * Public License (LGPL), Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.

 * THE SOFTWARE IS PROVIDED "AS IS." FRED HUTCHINSON CANCER RESEARCH CENTER MAKES NO
 * REPRESENTATIONS OR WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED,
 * INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
 * PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS,
 * WHETHER OR NOT DISCOVERABLE.  IN NO EVENT SHALL FRED HUTCHINSON CANCER RESEARCH
 * CENTER OR ITS TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR
 * ANY DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR
 * CONSEQUENTIAL DAMAGES, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS,
 * REGARDLESS OF  WHETHER FRED HUTCHINSON CANCER RESEARCH CENTER SHALL BE ADVISED,
 * SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE
 * FOREGOING.
 */

/*
 *
 */
package org.broad.igv.track;

import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.*;
import org.broad.igv.renderer.SpliceJunctionRenderer;
import org.broad.igv.sam.*;
import org.broad.igv.sam.reader.CachingQueryReader;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.renderer.DataRenderer;
import org.broad.igv.tdf.TDFDataSource;

import java.awt.*;
import java.util.*;
import java.util.List;

/**
 * @author dhmay
 */
public class SpliceJunctionFinderTrack extends AbstractTrack {

    private static Logger log = Logger.getLogger(SpliceJunctionFinderTrack.class);

    char[] nucleotides = {'a', 'c', 'g', 't', 'n'};
    public static Color lightBlue = new Color(0, 0, 150);
    private float[] bgColorComps = new float[3];

    AlignmentDataManager dataManager;
    SpliceJunctionRenderer spliceJunctionRenderer;
    PreferenceManager prefs;


    public SpliceJunctionFinderTrack(String id, String name) {
        super(id, name);
        super.setDataRange(new DataRange(0, 0, 60));
        spliceJunctionRenderer = new SpliceJunctionRenderer();
        prefs = PreferenceManager.getInstance();
    }

    public void setDataManager(AlignmentDataManager dataManager) {
        this.dataManager = dataManager;
    }


    @Override
    public void setDataRange(DataRange axisDefinition) {
        // Explicitly setting a data range turns off auto-scale
        super.setDataRange(axisDefinition);
    }


    public void render(RenderContext context, Rectangle rect) {

        //todo: add preference specifically for splice junctions
        float maxRange = PreferenceManager.getInstance().getAsFloat(PreferenceManager.SAM_MAX_VISIBLE_RANGE);
        float minVisibleScale = (maxRange * 1000) / 700;

        if (context.getScale() < minVisibleScale) {

            AlignmentInterval interval = null;
            if (dataManager != null) {
                interval = dataManager.getLoadedInterval(context);
            }


            if (interval != null && interval.contains(context.getGenomeId(), context.getChr(), (int) context.getOrigin(),
                    (int) context.getEndLocation())) {



                Map<Integer, Map<Integer, SpliceJunction>> posStartEndJunctionsMap =
                        new HashMap<Integer, Map<Integer, SpliceJunction>>();
                Map<Integer, Map<Integer, SpliceJunction>> negStartEndJunctionsMap =
                        new HashMap<Integer, Map<Integer, SpliceJunction>>();



                List<AlignmentInterval.Row> alignmentRows = interval.getAlignmentRows();

                for (AlignmentInterval.Row row : alignmentRows)
                {
                    //todo: is this dangerous?  Might something else depend on this not being reset?
                    row.resetIdx();

                    while (row.hasNext())
                    {
                        Alignment alignment = row.nextAlignment();



                        Map<Integer, Map<Integer, SpliceJunction>> startEndJunctionsMapThisStrand =
                                alignment.isNegativeStrand() ? negStartEndJunctionsMap : posStartEndJunctionsMap;
                        AlignmentBlock[] blocks = alignment.getAlignmentBlocks();
                        if (blocks.length < 2)
                            continue;
                        int prevBlockStart = -1;
                        int prevBlockEnd = -1;
                        for (AlignmentBlock block : blocks)
                        {
                            int blockEnd = block.getEnd();
                            int blockStart = block.getStart();
                            if (prevBlockEnd != -1)
                            {
                                Map<Integer, SpliceJunction> endJunctionsMap =
                                        startEndJunctionsMapThisStrand.get(prevBlockEnd);
                                if (endJunctionsMap == null)
                                {
                                    endJunctionsMap = new HashMap<Integer, SpliceJunction>();
                                    startEndJunctionsMapThisStrand.put(prevBlockEnd, endJunctionsMap);
                                }
                                SpliceJunction junction = endJunctionsMap.get(blockStart);
                                if (junction == null)
                                {
                                    junction = new SpliceJunction(prevBlockEnd, blockStart);
                                    endJunctionsMap.put(blockStart, junction);
                                }
                                junction.addRead(prevBlockStart, blockEnd);
                            }
                            prevBlockStart = blockStart;
                            prevBlockEnd = blockEnd;
                        }
                    }
                }

                List<IGVFeature> spliceJunctionFeatures = new ArrayList<IGVFeature>();
                //Build splice junction features for neg, then pos strand
                for (boolean isNeg : new boolean[] {true, false})
                {
                    Map<Integer, Map<Integer, SpliceJunction>> startEndJunctionsMapThisStrand =
                            isNeg ? negStartEndJunctionsMap : posStartEndJunctionsMap;

                    //note: this is designed to put the features in order.  But technically it won't -- they'll
                    //be in order by junction start, which isn't the start of the feature because it doesn't
                    //take into account flanking reads.  No problem for this renderer
                    for (int junctionStart : startEndJunctionsMapThisStrand.keySet())
                    {
                        Map<Integer, SpliceJunction> endJunctionMap = startEndJunctionsMapThisStrand.get(junctionStart);
                        for (int junctionEnd : endJunctionMap.keySet())
                        {
                            SpliceJunction junction = endJunctionMap.get(junctionEnd);
                            BasicFeature feature = new BasicFeature(context.getChr(),
                                    junction.flankingStart, junction.flankingEnd);
                            feature.setStrand(isNeg ? Strand.NEGATIVE : Strand.POSITIVE);
                            feature.setScore(junction.depth);
                            feature.addExon(new Exon(context.getChr(),
                                    junction.flankingStart, junction.start,
                                    isNeg ? Strand.NEGATIVE : Strand.POSITIVE
                            ));
                            feature.addExon(new Exon(context.getChr(),
                                    junction.end, junction.flankingEnd,
                                    isNeg ? Strand.NEGATIVE : Strand.POSITIVE
                            ));
                            spliceJunctionFeatures.add(feature);
                        }
                    }
                }

                Collections.sort(spliceJunctionFeatures, new Comparator<IGVFeature>()
                {
                    public int compare(IGVFeature o1, IGVFeature o2) {
                        return o1.getStart() - o2.getStart();
                    }
                });
                spliceJunctionRenderer.render(spliceJunctionFeatures, context, rect, this);
            }
        }
    }


    protected class SpliceJunction {
        int depth = 0;
        int start, end;
        //start of beginning flanking region, end of end flanking region
        int flankingStart, flankingEnd;

        public SpliceJunction(int start, int end) {
            this.start = start;
            this.end = end;
            this.flankingStart = start;            
            this.flankingEnd = end;
        }

        public void addRead(int readStart, int readEnd) {
            depth++;
            flankingStart = Math.min(readStart, flankingStart);
            flankingEnd = Math.max(readEnd, flankingEnd);
        }

    }


    public boolean isLogNormalized() {
        return false;
    }

    public float getRegionScore(String chr, int start, int end, int zoom, RegionScoreType type, ReferenceFrame frame) {
        return 0;
    }


}
