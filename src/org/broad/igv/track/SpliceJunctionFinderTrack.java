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
import org.broad.igv.ui.IGVMainFrame;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.renderer.DataRange;

import java.awt.*;
import java.io.IOException;
import java.util.*;
import java.util.List;

/**
 * @author dhmay
 */
public class SpliceJunctionFinderTrack extends FeatureTrack {

    private static Logger log = Logger.getLogger(SpliceJunctionFinderTrack.class);

    AlignmentDataManager dataManager;
    SpliceJunctionRenderer spliceJunctionRenderer;
    PreferenceManager prefs;
    RenderContext context;

    public SpliceJunctionFinderTrack(String id, String name, AlignmentDataManager dataManager) {
        super(id, name);
        init(new SpliceJunctionFinderFeatureSource());
        super.setDataRange(new DataRange(0, 0, 60));
        setRendererClass(SpliceJunctionRenderer.class);
        this.dataManager = dataManager;
        prefs = PreferenceManager.getInstance();
    }

    public void setDataManager(AlignmentDataManager dataManager) {
        this.dataManager = dataManager;
    }

    protected class SpliceJunctionFinderFeatureSource implements FeatureSource {

        public SpliceJunctionFinderFeatureSource()
        {
        }

        public Iterator getFeatures(String chr, int start, int end) throws IOException {
            List<IGVFeature> spliceJunctionFeatures = new ArrayList<IGVFeature>();

            if (!shouldShowFeatures())
                return null;

            AlignmentInterval interval = null;
            if (dataManager != null) {
                //This method is called in a Runnable in its own thread, so we can enter a long while loop here
                //and make sure the features get loaded
//                while (interval == null)
//                {
                    interval = dataManager.getLoadedInterval(context);
                    while (dataManager.isLoading())
                        try {

                            Thread.sleep(50);
                        }
                        catch (InterruptedException e) {}
                    interval = dataManager.getLoadedInterval(context);
//                }
                //sometimes, interval does not contain the requested interval.  In that case, we display no features.
                //I don't understand this.  Maybe dataManager is getting multiple requests with different extents?

//if (!interval.contains(context.getGenomeId(), context.getChr(), (int) context.getOrigin(),
//                    (int) context.getEndLocation()))
//{
//    RenderContext newContext = context;
//    new StringBuilder(oldContext + "" + newContext);
//}
            }
            if (interval == null)
                new StringBuilder();


            if (interval != null && interval.contains(context.getGenomeId(), context.getChr(), (int) context.getOrigin(),
                    (int) context.getEndLocation())) {

                Map<Integer, Map<Integer, SpliceJunctionFeature>> posStartEndJunctionsMap =
                        new HashMap<Integer, Map<Integer, SpliceJunctionFeature>>();
                Map<Integer, Map<Integer, SpliceJunctionFeature>> negStartEndJunctionsMap =
                        new HashMap<Integer, Map<Integer, SpliceJunctionFeature>>();



                List<AlignmentInterval.Row> alignmentRows = interval.getAlignmentRows();

                for (AlignmentInterval.Row row : alignmentRows)
                {
                    //todo: is this dangerous?  Might something else depend on this not being reset?
                    row.resetIdx();

                    while (row.hasNext())
                    {
                        Alignment alignment = row.nextAlignment();

                        boolean isNegativeStrand = alignment.isNegativeStrand();
                        Object strandAttr = alignment.getAttribute("XS");
                        if (strandAttr != null)
                            isNegativeStrand = strandAttr.toString().charAt(0) == '-';

                        Map<Integer, Map<Integer, SpliceJunctionFeature>> startEndJunctionsMapThisStrand =
                                isNegativeStrand ? negStartEndJunctionsMap : posStartEndJunctionsMap;
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
                                Map<Integer, SpliceJunctionFeature> endJunctionsMap =
                                        startEndJunctionsMapThisStrand.get(prevBlockEnd);
                                if (endJunctionsMap == null)
                                {
                                    endJunctionsMap = new HashMap<Integer, SpliceJunctionFeature>();
                                    startEndJunctionsMapThisStrand.put(prevBlockEnd, endJunctionsMap);
                                }
                                SpliceJunctionFeature junction = endJunctionsMap.get(blockStart);
                                if (junction == null)
                                {
                                    junction = new SpliceJunctionFeature(chr, prevBlockEnd, blockStart,
                                            isNegativeStrand ? Strand.NEGATIVE : Strand.POSITIVE);
                                    endJunctionsMap.put(blockStart, junction);
                                    spliceJunctionFeatures.add(junction);
                                }
                                junction.addRead(prevBlockStart, blockEnd);
                            }
                            prevBlockStart = blockStart;
                            prevBlockEnd = blockEnd;
                        }
                    }
                }

                Collections.sort(spliceJunctionFeatures, new Comparator<IGVFeature>()
                {
                    public int compare(IGVFeature o1, IGVFeature o2) {
                        return o1.getStart() - o2.getStart();
                    }
                });
//  if (spliceJunctionFeatures.isEmpty())
//      new StringBuilder("");
            }
            else
            {
//                new StringBuilder("");
            }
//         if (spliceJunctionFeatures.isEmpty())
//      new StringBuilder("").toString();
            return spliceJunctionFeatures.iterator();

        }

        public List<LocusScore> getCoverageScores(String chr, int start, int end, int zoom) {
            return null;  //To change body of implemented methods use File | Settings | File Templates.
        }

        public int getFeatureWindowSize() {
            return 0;  //To change body of implemented methods use File | Settings | File Templates.
        }

        public void setFeatureWindowSize(int size) {
            //To change body of implemented methods use File | Settings | File Templates.
        }

        public Class getFeatureClass() {
            return null;  //To change body of implemented methods use File | Settings | File Templates.
        }
    }

    /**
     * Relies on context
     * @return
     */
    protected boolean shouldShowFeatures()  {
        //todo: add preference specifically for splice junctions
        float maxRange = PreferenceManager.getInstance().getAsFloat(PreferenceManager.SAM_MAX_VISIBLE_RANGE);
        float minVisibleScale = (maxRange * 1000) / 700;

        if (context == null || (context.getScale() > minVisibleScale))
            return false;
        return true;
    }


    

        // Render features in the given input rectangle.

    protected void renderFeatures(RenderContext context, Rectangle inputRect) {

        if (featuresLoading) {
            return;
        }

        if (log.isDebugEnabled()) {
            log.debug("renderFeatures: " + getName());
        }

        String chr = context.getChr();
        int start = (int) context.getOrigin();
        int end = (int) context.getEndLocation() + 1;

        PackedFeatures packedFeatures = packedFeaturesMap.get(context.getReferenceFrame().getName());

        if (packedFeatures == null || !packedFeatures.containsInterval(chr, start, end)) {
            this.context = context;
            if (shouldShowFeatures()) {
                featuresLoading = true;
                loadFeatures(chr, start, end, context);

            }
            else if (packedFeatures != null)
                packedFeaturesMap.put(context.getReferenceFrame().getName(), null);
            if (!IGVMainFrame.getInstance().isExportingSnapshot()) {
                return;
            }
        }

       renderFeatureImpl(context, inputRect, packedFeatures);
    }

    protected void renderFeatureImpl(RenderContext context, Rectangle inputRect, PackedFeatures packedFeatures) {
        //todo: add preference specifically for splice junctions
        float maxRange = PreferenceManager.getInstance().getAsFloat(PreferenceManager.SAM_MAX_VISIBLE_RANGE);
        float minVisibleScale = (maxRange * 1000) / 700;

        if (context.getScale() > minVisibleScale)
            return;
        if (getDisplayMode() == DisplayMode.EXPANDED) {
            List<PackedFeatures.FeatureRow> rows = packedFeatures.getRows();
            if (rows != null && rows.size() > 0) {

                int nLevels = rows.size();
                synchronized (levelRects) {

                    levelRects.clear();

                    // Divide rectangle into equal height levels
                    double h = inputRect.getHeight() / nLevels;
                    Rectangle rect = new Rectangle(inputRect.x, inputRect.y, inputRect.width, (int) h);
                    int i = 0;
                    for (PackedFeatures.FeatureRow row : rows) {
                        levelRects.add(new Rectangle(rect));
                        getRenderer().render(row.features, context, rect, this);
                        if (selectedFeatureRowIndex == i) {
                            Graphics2D fontGraphics =
                                    (Graphics2D) context.getGraphic2DForColor(SELECTED_FEATURE_ROW_COLOR).create();
                            fontGraphics.fillRect(rect.x, rect.y, rect.width, rect.height);
                        }
                        rect.y += h;
                        i++;
                    }
                }
            }
        } else {
            List<IGVFeature> features = packedFeatures.getFeatures();
            if (features != null) {
                getRenderer().render(features, context, inputRect, this);
            }
        }
    }


    @Override
    public void setDataRange(DataRange axisDefinition) {
        // Explicitly setting a data range turns off auto-scale
        super.setDataRange(axisDefinition);
    }




    public boolean isLogNormalized() {
        return false;
    }

    public float getRegionScore(String chr, int start, int end, int zoom, RegionScoreType type, ReferenceFrame frame) {
        return 0;
    }


}
