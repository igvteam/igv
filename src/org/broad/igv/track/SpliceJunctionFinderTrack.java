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

import com.jidesoft.swing.JidePopupMenu;
import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.*;
import org.broad.igv.renderer.SpliceJunctionRenderer;
import org.broad.igv.sam.*;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.renderer.DataRange;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;
import java.util.*;
import java.util.List;

/**
 * @author dhmay
 * Finds splice junctions in real time and renders them as Features
 */
public class SpliceJunctionFinderTrack extends FeatureTrack {

    private static Logger log = Logger.getLogger(SpliceJunctionFinderTrack.class);

    AlignmentDataManager dataManager;
    PreferenceManager prefs;
    RenderContext context;

    protected int minReadFlankingWidth = 0;
    protected int minJunctionCoverage = 1;


    public SpliceJunctionFinderTrack(String id, String name, AlignmentDataManager dataManager) {
        super(id, name);
        init(new SpliceJunctionFinderFeatureSource());
        super.setDataRange(new DataRange(0, 0, 60));
        setRendererClass(SpliceJunctionRenderer.class);
        this.dataManager = dataManager;
        prefs = PreferenceManager.getInstance();
    }

    /**
     * Uses the datamanager to build SpliceJunctionFeatures for the alignments in a given range
     */
    protected class SpliceJunctionFinderFeatureSource implements FeatureSource {

        public SpliceJunctionFinderFeatureSource()
        {
        }

        public Iterator getFeatures(String chr, int start, int end) throws IOException {
            List<SpliceJunctionFeature> spliceJunctionFeatures = new ArrayList<SpliceJunctionFeature>();

            if (!shouldShowFeatures())
                return null;

            AlignmentInterval interval = null;
            if (dataManager != null) {
                //This method is called in a Runnable in its own thread, so we can enter a long while loop here
                //and make sure the features get loaded, without hanging the interface
                interval = dataManager.getLoadedInterval(context);
                if (dataManager.isLoading())
                {
                    while (dataManager.isLoading())
                        try {
                            Thread.sleep(50);
                        }
                        catch (InterruptedException e) {}
                    interval = dataManager.getLoadedInterval(context);
                }
            }
            //interval really shouldn't be null at this point

            if (interval != null && interval.contains(context.getGenomeId(), context.getChr(),
                    (int) context.getOrigin(), (int) context.getEndLocation())) {
                //we need to keep the positive and negative strand junctions separate, since
                //they don't represent the same thing and are rendered separately
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
                        //Any alignment with 2 or more blocks is considered to be a splice junction
                        Alignment alignment = row.nextAlignment();
                        AlignmentBlock[] blocks = alignment.getAlignmentBlocks();
                        if (blocks.length < 2)
                            continue;

                        //there may be other ways in which this is indicated. May have to code for them later
                        boolean isNegativeStrand = false;
                        Object strandAttr = alignment.getAttribute("XS");
                        if (strandAttr != null)
                            isNegativeStrand = strandAttr.toString().charAt(0) == '-';

                        Map<Integer, Map<Integer, SpliceJunctionFeature>> startEndJunctionsMapThisStrand =
                                isNegativeStrand ? negStartEndJunctionsMap : posStartEndJunctionsMap;

                        int flankingStart = -1;
                        int junctionStart = -1;
                        //for each pair of blocks, create or add evidence to a splice junction
                        for (AlignmentBlock block : blocks)
                        {
                            int flankingEnd = block.getEnd();
                            int junctionEnd = block.getStart();
                            if (junctionStart != -1)
                            {
                                //only proceed if the flanking regions are both bigger than the minimum
                                if (minReadFlankingWidth == 0 ||
                                        ((junctionStart - flankingStart >= minReadFlankingWidth) &&
                                        (flankingEnd - junctionEnd >= minReadFlankingWidth))) {
                                    Map<Integer, SpliceJunctionFeature> endJunctionsMap =
                                            startEndJunctionsMapThisStrand.get(junctionStart);
                                    if (endJunctionsMap == null)
                                    {
                                        endJunctionsMap = new HashMap<Integer, SpliceJunctionFeature>();
                                        startEndJunctionsMapThisStrand.put(junctionStart, endJunctionsMap);
                                    }
                                    SpliceJunctionFeature junction = endJunctionsMap.get(junctionEnd);
                                    if (junction == null)
                                    {
                                        junction = new SpliceJunctionFeature(chr, junctionStart, junctionEnd,
                                                isNegativeStrand ? Strand.NEGATIVE : Strand.POSITIVE);
                                        endJunctionsMap.put(junctionEnd, junction);
                                        spliceJunctionFeatures.add(junction);
                                    }
                                    junction.addRead(flankingStart, flankingEnd);
                                }
                            }
                            flankingStart = junctionEnd;
                            junctionStart = flankingEnd;
                        }
                    }
                }

                //get rid of any features without enough coverage
                if (minJunctionCoverage > 1)
                {
                    List<SpliceJunctionFeature> coveredFeatures =
                            new ArrayList<SpliceJunctionFeature>(spliceJunctionFeatures.size());
                    for (SpliceJunctionFeature feature : spliceJunctionFeatures)
                        if (feature.getJunctionDepth() >= minJunctionCoverage)
                            coveredFeatures.add(feature);
                    spliceJunctionFeatures = coveredFeatures;
                }

                //Sort by increasing beginning of start flanking region, as required by the renderer
                Collections.sort(spliceJunctionFeatures, new Comparator<IGVFeature>()
                {
                    public int compare(IGVFeature o1, IGVFeature o2) {
                        return o1.getStart() - o2.getStart();
                    }
                });
            }
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
     * Determine whether we should show the features. Relies on context
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



    /**
     * Render features in the given input rectangle.
     * @param context
     * @param inputRect
     */
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
            if (!IGV.getInstance().isExportingSnapshot()) {
                return;
            }
        }

       renderFeatureImpl(context, inputRect, packedFeatures);
    }

    /**
     * Render the track decorations and tell the renderer to render the features if appropriate
     * todo: make the horizontal center line appear even if not rendering features
     * @param context
     * @param inputRect
     * @param packedFeatures
     */
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

    /**
     * Add a MenuItem to control the minimum flanking width for reads used in finding junctions.
     * If either the start OR the end flanking region is less than this, the read is not used
     * @param menu
     * @return
     */
    public JMenuItem addFlankingWidthTresholdItem(JPopupMenu menu) {
        JMenuItem flankingWidthItem = new JMenuItem("Set minimum read flanking width...");

        flankingWidthItem.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {

                String value = JOptionPane.showInputDialog("Minimum start and end flanking region width: ",
                        Float.valueOf(minReadFlankingWidth));
                if (value == null) {
                    return;
                }
                try {
                    int tmp = Integer.parseInt(value);
                    minReadFlankingWidth = tmp;
                    IGV.getInstance().repaintDataPanels();
                }
                catch (Exception exc) {
                    //log
                }

            }
        });
        menu.add(flankingWidthItem);

        return flankingWidthItem;
    }


    /**
     * Add a MenuItem to control the minimum flanking width for reads used in finding junctions.
     * If either the start OR the end flanking region is less than this, the read is not used
     * @param menu
     * @return
     */
    public JMenuItem addJunctionCoverageTresholdItem(JPopupMenu menu) {
        JMenuItem junctionDepthItem = new JMenuItem("Set minimum junction coverage...");

        junctionDepthItem.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {

                String value = JOptionPane.showInputDialog("Minimum coverage depth for displayed junctions: ",
                        Float.valueOf(minJunctionCoverage));
                if (value == null) {
                    return;
                }
                try {
                    int tmp = Integer.parseInt(value);
                    minJunctionCoverage = tmp;
                    IGV.getInstance().repaintDataPanels();
                }
                catch (Exception exc) {
                    //log
                }

            }
        });
        menu.add(junctionDepthItem);

        return junctionDepthItem;
    }

    /**
     * Override to return a specialized popup menu
     *
     * @return
     */
    @Override
    public JPopupMenu getPopupMenu(TrackClickEvent te) {

        JPopupMenu popupMenu = new JidePopupMenu();

        JLabel popupTitle = new JLabel("  " + getName(), JLabel.CENTER);

        Font newFont = popupMenu.getFont().deriveFont(Font.BOLD, 12);
        popupTitle.setFont(newFont);
        if (popupTitle != null) {
            popupMenu.add(popupTitle);
        }

        popupMenu.addSeparator();

        ArrayList<Track> tmp = new ArrayList();
        tmp.add(this);
        TrackMenuUtils.addStandardItems(popupMenu, tmp, te);
        popupMenu.addSeparator();
        addFlankingWidthTresholdItem(popupMenu);
        addJunctionCoverageTresholdItem(popupMenu);
        return popupMenu;
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
