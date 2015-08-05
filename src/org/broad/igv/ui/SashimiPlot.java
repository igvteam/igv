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

        initSize(frame.getWidthInPixels());

        BoxLayout boxLayout = new BoxLayout(getContentPane(), BoxLayout.Y_AXIS);
        getContentPane().setLayout(boxLayout);

        //Add control elements to the top
        getContentPane().add(generateControlPanel(this.frame));

        spliceJunctionTracks = new ArrayList<SpliceJunctionFinderTrack>(alignmentTracks.size());
        int colorInd = 0;

        for (AlignmentTrack alignmentTrack : alignmentTracks) {


            AlignmentDataManager oldDataManager = alignmentTrack.getDataManager();
            MemoryAlignmentDataManager dataManager = new MemoryAlignmentDataManager(oldDataManager, oldDataManager.getSpliceJunctionLoadOptions());

            SpliceJunctionFinderTrack spliceJunctionTrack =
                    new SpliceJunctionFinderTrack(alignmentTrack.getResourceLocator(), alignmentTrack.getName(), dataManager, SpliceJunctionFinderTrack.StrandOption.COMBINE);
            // Override expand/collpase setting -- expanded sashimi plots make no sense
            spliceJunctionTrack.setDisplayMode(Track.DisplayMode.COLLAPSED);

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

    private Component generateControlPanel(ReferenceFrame frame) {
        JPanel controlPanel = new JPanel();

        ZoomSliderPanel zoomSliderPanel = new ZoomSliderPanel(frame);
        zoomSliderPanel.setMinZoomLevel(frame.getZoom());

        zoomSliderPanel.addMouseListener(new MouseAdapter() {
            @Override
            public void mouseClicked(MouseEvent e) {
                SashimiPlot.this.repaint();
            }
        });

        Dimension controlSize = new Dimension(200, 30);

        //JSlider scaleSlider = new JSlider(JSlider.HORIZONTAL);
        //setFixedSize(scaleSlider, controlSize);
        //controlPanel.add(scaleSlider);

        controlPanel.add(zoomSliderPanel);
        setFixedSize(zoomSliderPanel, controlSize);

        Dimension panelSize = controlSize;
        setFixedSize(controlPanel, panelSize);


        BoxLayout layout = new BoxLayout(controlPanel, BoxLayout.X_AXIS);
        controlPanel.setLayout(layout);

        return controlPanel;
    }

    private static void setFixedSize(Component component, Dimension dimension) {
        component.setPreferredSize(dimension);
        component.setMinimumSize(dimension);
        component.setMaximumSize(dimension);
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

        geneTrack.clearPackedFeatures();
        RenderContext context = new RenderContextImpl(geneComponent, null, frame, null);
        geneTrack.setForceLoadSync(true);
        geneTrack.load(context);


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
        getRenderer(trackComponent.track).getCoverageTrack().rescale(trackComponent.frame);

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
        private String toolTipText = null;

        public TrackComponent(ReferenceFrame frame, T track) {
            this.frame = frame;
            this.track = track;
        }

        public void updateToolTipText(TrackClickEvent tce) {
            toolTipText = track.getValueStringAt(tce.getFrame().getChrName(), tce.getChromosomePosition(), tce.getMouseEvent().getY(), tce.getFrame());
            toolTipText = "<html>" + toolTipText;
            setToolTipText(toolTipText);
        }

        @Override
        public void paintComponent(Graphics g) {
            super.paintComponent(g);
            Rectangle visibleRect = getVisibleRect();
            RenderContext context = new RenderContextImpl(this, (Graphics2D) g, frame, visibleRect);
            track.render(context, visibleRect);
        }

    }

    /**
     * Set the minimum junction coverage, per trac,k and is not persistent
     * <p/>
     * Our "Set Max Junction Coverage Range" just changes the view scaling, it doesn't
     * filter anything, which is different behavior than the minimum. This might be confusing.
     *
     * @param trackComponent
     * @param newMinJunctionCoverage
     */
    private void setMinJunctionCoverage(TrackComponent<SpliceJunctionFinderTrack> trackComponent, int newMinJunctionCoverage) {
        IAlignmentDataManager dataManager = getRenderer(trackComponent.track).getDataManager();
        dataManager.setMinJunctionCoverage(newMinJunctionCoverage);
        trackComponent.track.clear();
        trackComponent.repaint();
    }


    /**
     * Set the max coverage depth, which is a graphical scaling parameter for determining how
     * thick the junction arcs will be
     *
     * @param trackComponent
     * @param newMaxDepth
     */
    private void setMaxCoverageDepth(TrackComponent<SpliceJunctionFinderTrack> trackComponent, int newMaxDepth) {
        getRenderer(trackComponent.track).setMaxDepth(newMaxDepth);
        repaint();
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

            final JCheckBoxMenuItem showCoverageData = new JCheckBoxMenuItem("Show Exon Coverage Data");
            showCoverageData.setSelected(PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SASHIMI_SHOW_COVERAGE));
            showCoverageData.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    PreferenceManager.getInstance().put(PreferenceManager.SASHIMI_SHOW_COVERAGE, showCoverageData.isSelected());
                    SashimiPlot.this.repaint();
                }
            });

            CoverageTrack covTrack = getRenderer(this.trackComponent.track).getCoverageTrack();
            JMenuItem setCoverageDataRange = CoverageTrack.addDataRangeItem(SashimiPlot.this, null, Arrays.asList(covTrack));
            setCoverageDataRange.setText("Set Exon Coverage Max");

            JMenuItem maxJunctionCoverageRange = new JMenuItem("Set Junction Coverage Max");
            maxJunctionCoverageRange.setToolTipText("The thickness of each line will be proportional to the coverage, up until this value");

            maxJunctionCoverageRange.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    String input = JOptionPane.showInputDialog("Set Max Junction Coverage", getRenderer(trackComponent.track).getMaxDepth());
                    if (input == null || input.length() == 0) return;
                    try {
                        int newMaxDepth = Integer.parseInt(input);
                        setMaxCoverageDepth(JunctionTrackMouseAdapter.this.trackComponent, newMaxDepth);
                    } catch (NumberFormatException ex) {
                        JOptionPane.showMessageDialog(SashimiPlot.this, input + " is not an integer");
                    }
                }
            });
            JMenuItem minJunctionCoverage = new JMenuItem("Set Junction Coverage Min");
            minJunctionCoverage.setToolTipText("Junctions below this threshold will be removed from view");
            minJunctionCoverage.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    IAlignmentDataManager dataManager = getRenderer(trackComponent.track).getDataManager();
                    SpliceJunctionHelper.LoadOptions loadOptions = dataManager.getSpliceJunctionLoadOptions();

                    String input = JOptionPane.showInputDialog("Set Minimum Junction Coverage", loadOptions.minJunctionCoverage);
                    if (input == null || input.length() == 0) return;
                    try {
                        int newMinJunctionCoverage = Integer.parseInt(input);
                        setMinJunctionCoverage(JunctionTrackMouseAdapter.this.trackComponent, newMinJunctionCoverage);

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

            ButtonGroup shapeGroup = new ButtonGroup();

            JRadioButtonMenuItem textShape = getJunctionCoverageRadioButton("Text", SashimiJunctionRenderer.ShapeType.TEXT);
            textShape.setToolTipText("Show junction coverage as text");
            shapeGroup.add(textShape);

            JRadioButtonMenuItem circleShape = getJunctionCoverageRadioButton("Circle", SashimiJunctionRenderer.ShapeType.CIRCLE);
            circleShape.setToolTipText("Show junction coverage as a circle");
            shapeGroup.add(circleShape);

            JRadioButtonMenuItem ellipseShape = getJunctionCoverageRadioButton("Ellipse", SashimiJunctionRenderer.ShapeType.ELLIPSE);
            ellipseShape.setToolTipText("Show junction coverage as an ellipse");
            shapeGroup.add(ellipseShape);

            JRadioButtonMenuItem noShape = getJunctionCoverageRadioButton("None", SashimiJunctionRenderer.ShapeType.NONE);
            ellipseShape.setToolTipText("Show junction coverage as an ellipse");
            shapeGroup.add(ellipseShape);


            ButtonGroup strandGroup = new ButtonGroup();

            JRadioButtonMenuItem combineStrands = getStrandRadioButton("Combine Strands", SpliceJunctionFinderTrack.StrandOption.COMBINE);
            combineStrands.setToolTipText("Combine junctions from both strands -- best for non-strand preserving libraries.");
            strandGroup.add(combineStrands);

            //  JRadioButtonMenuItem bothStrands = getStrandRadioButton("Both Strands", SpliceJunctionFinderTrack.StrandOption.BOTH);
            //  strandGroup.add(bothStrands);
            //  menu.add(bothStrands);

            JRadioButtonMenuItem plusStrand = getStrandRadioButton("Forward Strand", SpliceJunctionFinderTrack.StrandOption.FORWARD);
            plusStrand.setToolTipText("Show only junctions on the forward read strand  (of first-in-pair for paired reads)");
            strandGroup.add(plusStrand);

            JRadioButtonMenuItem minusStrand = getStrandRadioButton("Reverse Strand", SpliceJunctionFinderTrack.StrandOption.REVERSE);
            plusStrand.setToolTipText("Show only junctions on the reverse read strand  (of first-in-pair for paired reads)");
            strandGroup.add(minusStrand);


            JMenuItem saveImageItem = new JMenuItem("Save Image...");
            saveImageItem.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    File defaultFile = new File("Sashimi.png");
                    IGV.getInstance().createSnapshot(SashimiPlot.this.getContentPane(), defaultFile);
                }
            });

            // Coverage ranges -- these apply to current plot only
            menu.add(new JLabel("Junction Coverage Display"));
            menu.add(setCoverageDataRange);
            menu.add(minJunctionCoverage);
            menu.add(maxJunctionCoverageRange);
            menu.add(colorItem);

            // Coverage data  -- applies to all plots
            menu.add(showCoverageData);

            // Shape options -- all plots
            menu.addSeparator();
            menu.add(textShape);
            menu.add(circleShape);
           // menu.add(ellipseShape);
            menu.add(noShape);

            // Strand options -- applies to all plots
            menu.addSeparator();
            menu.add(combineStrands);
            menu.add(plusStrand);
            menu.add(minusStrand);
            menu.addSeparator();

            menu.add(saveImageItem);

            return menu;
        }

        private JRadioButtonMenuItem getJunctionCoverageRadioButton(String label, final SashimiJunctionRenderer.ShapeType shapeType) {


            JRadioButtonMenuItem item = new JRadioButtonMenuItem(label);
            item.setSelected(SashimiJunctionRenderer.getShapeType() == shapeType);
            item.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    SashimiJunctionRenderer.setShapeType(shapeType);
                    SashimiPlot.this.repaint();
                }
            });
            return item;
        }

        private JRadioButtonMenuItem getStrandRadioButton(String label, final SpliceJunctionFinderTrack.StrandOption option) {

            JRadioButtonMenuItem item = new JRadioButtonMenuItem(label);
            item.setSelected(SpliceJunctionFinderTrack.getStrandOption() == option);
            item.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    SpliceJunctionFinderTrack.setStrandOption(option);
                    for (SpliceJunctionFinderTrack t : spliceJunctionTracks) {
                        t.clear();
                    }
                    repaint();
                }
            });
            return item;
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
        public void mouseMoved(MouseEvent e) {
            trackComponent.updateToolTipText(createTrackClickEvent(e));
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

            List<FeatureTrack> nonJunctionTracks = new ArrayList(IGV.getInstance().getFeatureTracks());
            Iterator iter = nonJunctionTracks.iterator();
            while(iter.hasNext()) {
                if(iter.next() instanceof SpliceJunctionFinderTrack) iter.remove();
            }

            if (nonJunctionTracks.size() == 1) {
                geneTrack = nonJunctionTracks.get(0);
            } else {
                FeatureTrackSelectionDialog dlg = new FeatureTrackSelectionDialog(IGV.getMainFrame(), nonJunctionTracks);
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
                TrackSelectionDialog<AlignmentTrack> alDlg =
                        new TrackSelectionDialog<AlignmentTrack>(IGV.getMainFrame(), TrackSelectionDialog.SelectionMode.MULTIPLE, alignmentTracks);
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
