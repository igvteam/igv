/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.ui;

import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.IExon;
import org.broad.igv.renderer.SashimiJunctionRenderer;
import org.broad.igv.sam.*;
import org.broad.igv.track.*;
import org.broad.igv.ui.color.ColorPalette;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.ui.event.AlignmentTrackEvent;
import org.broad.igv.ui.panel.*;
import org.broad.igv.ui.util.UIUtilities;

import javax.swing.*;
import javax.swing.event.PopupMenuEvent;
import javax.swing.event.PopupMenuListener;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.geom.Rectangle2D;
import java.io.File;
import java.util.*;
import java.util.List;

/**
 * Window for displaying sashimi style junction plot
 * See http://genes.mit.edu/burgelab/miso/docs/sashimi.html
 * <p/>
 * User: jacob
 * Date: 2013-Jan-11
 */
public class SashimiPlot extends JFrame {

    private List<SpliceJunctionFinderTrack> spliceJunctionTracks;
    private ReferenceFrame frame;
    private AlignmentDataManager dataManager;  // <= retained to release self from event bus

    /**
     * The minimum allowed origin of the frame. We set scrolling
     * limits based on initialization
     */
    private final double minOrigin;

    /**
     * The maximum allow end of the frame. We set scrolling
     * limits based on initialization
     */
    private final double maxEnd;



    private static final List<Color> plotColors;

    static {
        ColorPalette palette = ColorUtilities.getDefaultPalette();
        plotColors = Arrays.asList(palette.getColors());
    }


    public SashimiPlot(ReferenceFrame iframe, Collection<? extends AlignmentTrack> alignmentTracks, FeatureTrack geneTrack) {
        getGlassPane().setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
        int minJunctionCoverage = PreferenceManager.getInstance().getAsInt(PreferenceManager.SAM_JUNCTION_MIN_COVERAGE);

        this.frame = new ReferenceFrame(iframe);

        minOrigin = this.frame.getOrigin();
        maxEnd = this.frame.getEnd();
        int minZoom = this.frame.getZoom();

        initSize(frame.getWidthInPixels());

        BoxLayout boxLayout = new BoxLayout(getContentPane(), BoxLayout.Y_AXIS);
        getContentPane().setLayout(boxLayout);

        //Add control elements to the top
        ZoomSliderPanel controlPanel = new ZoomSliderPanel(this.frame);
        controlPanel.setMinZoomLevel(minZoom);

        Dimension dimSize = new Dimension(200, 30);
        controlPanel.setPreferredSize(dimSize);
        controlPanel.setMinimumSize(dimSize);
        controlPanel.setMaximumSize(dimSize);

        getContentPane().add(controlPanel);
        controlPanel.addMouseListener(new MouseAdapter() {
            @Override
            public void mouseClicked(MouseEvent e) {
                SashimiPlot.this.repaint();
            }
        });

        spliceJunctionTracks = new ArrayList<SpliceJunctionFinderTrack>(alignmentTracks.size());
        int colorInd = 0;

        for(AlignmentTrack alignmentTrack: alignmentTracks){


            AlignmentDataManager oldDataManager = alignmentTrack.getDataManager();
            MemoryAlignmentDataManager dataManager = new MemoryAlignmentDataManager(oldDataManager, oldDataManager.getSpliceJunctionLoadOptions());

            SpliceJunctionFinderTrack spliceJunctionTrack = new SpliceJunctionFinderTrack(alignmentTrack.getResourceLocator(), alignmentTrack.getName(), dataManager, true);

            spliceJunctionTrack.setRendererClass(SashimiJunctionRenderer.class);

            Color color = plotColors.get(colorInd);
            colorInd = (colorInd + 1) % plotColors.size();
            spliceJunctionTrack.setColor(color);

            TrackComponent<SpliceJunctionFinderTrack> trackComponent = new TrackComponent<SpliceJunctionFinderTrack>(frame, spliceJunctionTrack);

            initSpliceJunctionComponent(trackComponent, dataManager, oldDataManager.getCoverageTrack(), minJunctionCoverage);

            getContentPane().add(trackComponent);
            spliceJunctionTracks.add(spliceJunctionTrack);
        }

        Axis axis = createAxis(frame);
        getContentPane().add(axis);

        SelectableFeatureTrack geneTrackClone = new SelectableFeatureTrack(geneTrack);
        TrackComponent<SelectableFeatureTrack> geneComponent = new TrackComponent<SelectableFeatureTrack>(frame, geneTrackClone);


        getContentPane().add(geneComponent);

        initGeneComponent(frame.getWidthInPixels(), geneComponent, geneTrackClone);
        validate();
    }

    private void initSize(int width) {
        setSize(width, 500);
    }

    private Axis createAxis(ReferenceFrame frame) {
        Axis axis = new Axis(frame);

        Dimension maxDim = new Dimension(Integer.MAX_VALUE, 25);
        axis.setMaximumSize(maxDim);
        Dimension prefDim = new Dimension(maxDim);
        prefDim.setSize(frame.getWidthInPixels(), prefDim.height);
        axis.setPreferredSize(prefDim);

        return axis;
    }

    private void initGeneComponent(int prefWidth, TrackComponent<SelectableFeatureTrack> geneComponent, FeatureTrack geneTrack) {

        geneTrack.setDisplayMode(Track.DisplayMode.SQUISHED);

        //Hacky way of clearing packed features
        geneTrack.setVisibilityWindow(geneTrack.getVisibilityWindow());
        RenderContext context = new RenderContextImpl(geneComponent, null, frame, null);
        geneTrack.preload(context);


        Dimension maxGeneDim = new Dimension(Integer.MAX_VALUE, geneTrack.getNumberOfFeatureLevels() * geneTrack.getSquishedRowHeight() + 10);
        geneComponent.setMaximumSize(maxGeneDim);
        Dimension prefGeneDim = new Dimension(maxGeneDim);
        prefGeneDim.setSize(prefWidth, prefGeneDim.height);
        geneComponent.setPreferredSize(prefGeneDim);

        GeneTrackMouseAdapter ad2 = new GeneTrackMouseAdapter(geneComponent);
        geneComponent.addMouseListener(ad2);
        geneComponent.addMouseMotionListener(ad2);
    }

    private void initSpliceJunctionComponent(TrackComponent<SpliceJunctionFinderTrack> trackComponent, IAlignmentDataManager dataManager, CoverageTrack coverageTrack, int minJunctionCoverage) {
        JunctionTrackMouseAdapter ad1 = new JunctionTrackMouseAdapter(trackComponent);
        trackComponent.addMouseListener(ad1);
        trackComponent.addMouseMotionListener(ad1);

        getRenderer(trackComponent.track).setDataManager(dataManager);
        getRenderer(trackComponent.track).setCoverageTrack(coverageTrack);

        dataManager.setMinJunctionCoverage(minJunctionCoverage);

        getRenderer(trackComponent.track).setBackground(getBackground());
    }

    private SashimiJunctionRenderer getRenderer(SpliceJunctionFinderTrack spliceJunctionTrack) {
        return (SashimiJunctionRenderer) spliceJunctionTrack.getRenderer();
    }

    /**
     * Should consider using this elsewhere. Single component
     * which contains a single track
     */
    private static class TrackComponent<T extends Track> extends JComponent {

        private T track;
        private ReferenceFrame frame;

        public TrackComponent(ReferenceFrame frame, T track) {
            this.frame = frame;
            this.track = track;
        }


        @Override
        public void paintComponent(Graphics g) {
            super.paintComponent(g);
            Rectangle visibleRect = getVisibleRect();
            RenderContext context = new RenderContextImpl(this, (Graphics2D) g, frame, visibleRect);
            track.render(context, visibleRect);
        }

    }

    private class JunctionTrackMouseAdapter extends TrackComponentMouseAdapter<SpliceJunctionFinderTrack> {

        JunctionTrackMouseAdapter(TrackComponent<SpliceJunctionFinderTrack> trackComponent) {
            super(trackComponent);
        }

        @Override
        protected void handleDataClick(MouseEvent e) {
            //Show data of some sort?
        }

        @Override
        protected IGVPopupMenu getPopupMenu(MouseEvent e) {
            IGVPopupMenu menu = new IGVPopupMenu();

            CoverageTrack covTrack = getRenderer(this.trackComponent.track).getCoverageTrack();
            JMenuItem setCoverageDataRange = CoverageTrack.addDataRangeItem(SashimiPlot.this, null, Arrays.asList(covTrack));
            setCoverageDataRange.setText("Set Coverage Data Range");
            menu.add(setCoverageDataRange);

            JMenuItem minJunctionCoverage = new JMenuItem("Set Min Junction Coverage");
            minJunctionCoverage.setToolTipText("Junctions below this threshold will be removed from view");
            minJunctionCoverage.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    /** The popup sets it per track and is not persistent
                     *
                     * On top of this, our "Set Max Junction Coverage Range" just changes the view scaling, it doesn't
                     * filter anything, which is different behavior than the minimum. This might be confusing.
                     *
                     * Did some refactoring, going to set it for the track and clear that tracks data, but the setting
                     * will not persist at all. May be weird for users. Still has the problem that max/min do very different
                     * things.
                     */

                    IAlignmentDataManager dataManager = getRenderer(trackComponent.track).getDataManager();
                    SpliceJunctionHelper.LoadOptions loadOptions = dataManager.getSpliceJunctionLoadOptions();

                    String input = JOptionPane.showInputDialog("Set Minimum Junction Coverage", loadOptions.minJunctionCoverage);
                    if (input == null || input.length() == 0) return;
                    try {
                        int newMinJunctionCoverage = Integer.parseInt(input);
                        dataManager.setMinJunctionCoverage(newMinJunctionCoverage);

                        trackComponent.track.onAlignmentTrackEvent(new AlignmentTrackEvent(this, AlignmentTrackEvent.Type.SPLICE_JUNCTION));
                        trackComponent.repaint();
                    } catch (NumberFormatException ex) {
                        JOptionPane.showMessageDialog(SashimiPlot.this, input + " is not an integer");
                    }


                }
            });

            JMenuItem maxJunctionCoverageRange = new JMenuItem("Set Max Junction Coverage Range");
            maxJunctionCoverageRange.setToolTipText("The thickness of each line will be proportional to the coverage, up until this value");

            maxJunctionCoverageRange.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    String input = JOptionPane.showInputDialog("Set Max Junction Coverage", getRenderer(trackComponent.track).getMaxDepth());
                    if (input == null || input.length() == 0) return;
                    try {
                        int newMaxDepth = Integer.parseInt(input);
                        getRenderer(trackComponent.track).setMaxDepth(newMaxDepth);
                        repaint();
                    } catch (NumberFormatException ex) {
                        JOptionPane.showMessageDialog(SashimiPlot.this, input + " is not an integer");
                    }
                }
            });

            JMenuItem colorItem = new JMenuItem("Set Color");
            colorItem.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    Color color = UIUtilities.showColorChooserDialog(
                            "Select Track Color", trackComponent.track.getColor());
                    SashimiPlot.this.toFront();
                    if (color == null) return;
                    trackComponent.track.setColor(color);
                    trackComponent.repaint();
                }
            });

            JMenuItem saveImageItem = new JMenuItem("Save Image...");
            saveImageItem.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    File defaultFile = new File("Sashimi.png");
                    IGV.getInstance().createSnapshot(SashimiPlot.this.getContentPane(), defaultFile);
                }
            });

            menu.add(minJunctionCoverage);
            menu.add(maxJunctionCoverageRange);
            menu.add(colorItem);
            menu.add(saveImageItem);

            return menu;
        }
    }

    private class GeneTrackMouseAdapter extends TrackComponentMouseAdapter<SelectableFeatureTrack> {

        GeneTrackMouseAdapter(TrackComponent<SelectableFeatureTrack> trackComponent) {
            super(trackComponent);
        }

        @Override
        protected void handleDataClick(MouseEvent e) {
            trackComponent.track.handleDataClick(createTrackClickEvent(e));
            Set<IExon> selectedExon = trackComponent.track.getSelectedExons();
            for (SpliceJunctionFinderTrack spliceTrack : spliceJunctionTracks) {
                getRenderer(spliceTrack).setSelectedExons(selectedExon);
            }
            repaint();
        }

        @Override
        protected IGVPopupMenu getPopupMenu(MouseEvent e) {
            IGVPopupMenu menu = new IGVPopupMenu();
            TrackMenuUtils.addDisplayModeItems(Arrays.<Track>asList(trackComponent.track), menu);
            menu.addPopupMenuListener(new RepaintPopupMenuListener(SashimiPlot.this));
            return menu;
        }
    }

    private abstract class TrackComponentMouseAdapter<T extends Track> extends MouseAdapter {

        protected TrackComponent<T> trackComponent;
        protected PanTool currentTool;

        TrackComponentMouseAdapter(TrackComponent<T> trackComponent) {
            this.trackComponent = trackComponent;
            currentTool = new PanTool(null);
            currentTool.setReferenceFrame(this.trackComponent.frame);
        }


        @Override
        public void mouseDragged(MouseEvent e) {
            if (currentTool.getLastMousePoint() == null) {
                //This shouldn't happen, but does occasionally
                return;
            }
            double diff = e.getX() - currentTool.getLastMousePoint().getX();
            // diff > 0 means moving mouse to the right, which drags the frame towards the negative direction
            boolean hitBounds = SashimiPlot.this.frame.getOrigin() <= minOrigin && diff > 0;
            hitBounds |= SashimiPlot.this.frame.getEnd() >= maxEnd && diff < 0;
            if (!hitBounds) {
                currentTool.mouseDragged(e);
                repaint();
            }
        }

        @Override
        public void mouseReleased(MouseEvent e) {
            if (e.isPopupTrigger()) {
                doPopupMenu(e);
            } else {
                currentTool.mouseReleased(e);
            }
        }

        @Override
        public void mousePressed(MouseEvent e) {
            if (e.isPopupTrigger()) {
                doPopupMenu(e);
            } else {
                currentTool.mousePressed(e);
                super.mousePressed(e);
            }

        }

        protected void doPopupMenu(MouseEvent e) {
            IGVPopupMenu menu = getPopupMenu(e);
            if (menu != null) menu.show(trackComponent, e.getX(), e.getY());
        }

        protected TrackClickEvent createTrackClickEvent(MouseEvent e) {
            return new TrackClickEvent(e, trackComponent.frame);
        }

        @Override
        public void mouseClicked(MouseEvent e) {
            if (e.isPopupTrigger()) {
                doPopupMenu(e);
                return;
            }

            currentTool.mouseClicked(e);
            handleDataClick(e);
        }

        /**
         * Essentially left click
         *
         * @param e
         */
        protected abstract void handleDataClick(MouseEvent e);

        /**
         * Essentially right click
         *
         * @param e
         * @return
         */
        protected abstract IGVPopupMenu getPopupMenu(MouseEvent e);
    }

    /**
     * Show SashimiPlot window, or change settings of {@code currentWindow}
     *
     * @param sashimiPlot
     */
    public static void getSashimiPlot(SashimiPlot sashimiPlot) {
        if (sashimiPlot == null) {
            FeatureTrack geneTrack = null;
            if (IGV.getInstance().getFeatureTracks().size() == 1) {
                geneTrack = IGV.getInstance().getFeatureTracks().get(0);
            } else {
                FeatureTrackSelectionDialog dlg = new FeatureTrackSelectionDialog(IGV.getMainFrame());
                dlg.setTitle("Select Gene Track");
                dlg.setVisible(true);
                if (dlg.getIsCancelled()) return;
                geneTrack = dlg.getSelectedTrack();
            }

            Collection<AlignmentTrack> alignmentTracks = new ArrayList<AlignmentTrack>();
            for (Track track : IGV.getInstance().getAllTracks()) {
                if (track instanceof AlignmentTrack) {
                    alignmentTracks.add((AlignmentTrack) track);
                }
            }

            if (alignmentTracks.size() > 1) {
                TrackSelectionDialog<AlignmentTrack> alDlg = new TrackSelectionDialog<AlignmentTrack>(IGV.getMainFrame(), TrackSelectionDialog.SelectionMode.MULTIPLE, alignmentTracks);
                alDlg.setTitle("Select Alignment Tracks");
                alDlg.setVisible(true);
                if (alDlg.getIsCancelled()) return;

                alignmentTracks = alDlg.getSelectedTracks();
            }

            sashimiPlot = new SashimiPlot(FrameManager.getDefaultFrame(), alignmentTracks, geneTrack);
            //sashimiPlot.setShapeType(shapeType);
            sashimiPlot.setVisible(true);
        } else {
            //sashimiPlot.setShapeType(shapeType);
        }


    }

    private static class RepaintPopupMenuListener implements PopupMenuListener {

        Component component;

        RepaintPopupMenuListener(Component component) {
            this.component = component;
        }

        @Override
        public void popupMenuWillBecomeVisible(PopupMenuEvent e) {
            component.repaint();
        }

        @Override
        public void popupMenuWillBecomeInvisible(PopupMenuEvent e) {
            component.repaint();
        }

        @Override
        public void popupMenuCanceled(PopupMenuEvent e) {
            component.repaint();
        }
    }

    private static class Axis extends JComponent {

        private ReferenceFrame frame;

        Axis(ReferenceFrame frame) {
            this.frame = frame;
        }

        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            Rectangle visibleRect = getVisibleRect();
            RenderContext context = new RenderContextImpl(this, (Graphics2D) g, frame, visibleRect);
            drawGenomicAxis(context, visibleRect);
        }

        /**
         * Draw axis displaying genomic coordinates
         *
         * @param context
         * @param trackRectangle
         */
        private void drawGenomicAxis(RenderContext context, Rectangle trackRectangle) {
            int numTicks = 4;
            int ticHeight = 5;

            double pixelPadding = trackRectangle.getWidth() / 20;
            int yLoc = ticHeight + 1;

            double origin = context.getOrigin();
            double locScale = context.getScale();

            //Pixel start/end positions of ruler
            double startPix = trackRectangle.getX() + pixelPadding;
            double endPix = trackRectangle.getMaxX() - pixelPadding;

            double ticIntervalPix = (endPix - startPix) / (numTicks - 1);
            double ticIntervalCoord = locScale * ticIntervalPix;

            int startCoord = (int) (origin + (locScale * startPix));

            Graphics2D g2D = context.getGraphic2DForColor(Color.black);

            g2D.drawLine((int) startPix, yLoc, (int) endPix, yLoc);

            for (int tic = 0; tic < numTicks; tic++) {
                int xLoc = (int) (startPix + tic * ticIntervalPix);
                g2D.drawLine(xLoc, yLoc, xLoc, yLoc - ticHeight);

                int ticCoord = (int) (startCoord + tic * ticIntervalCoord);
                String text = "" + ticCoord;
                Rectangle2D textBounds = g2D.getFontMetrics().getStringBounds(text, g2D);
                g2D.drawString(text, (int) (xLoc - textBounds.getWidth() / 2), (int) (yLoc + textBounds.getHeight()));
            }


        }
    }
}
