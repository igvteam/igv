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

package org.broad.igv.sashimi;

import org.broad.igv.event.IGVEvent;
import org.broad.igv.event.IGVEventBus;
import org.broad.igv.event.IGVEventObserver;
import org.broad.igv.event.ViewChange;
import org.broad.igv.feature.IExon;
import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.sam.AlignmentDataManager;
import org.broad.igv.sam.AlignmentTrack;
import org.broad.igv.sam.CoverageTrack;
import org.broad.igv.sam.SpliceJunctionTrack;
import org.broad.igv.track.*;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.color.ColorPalette;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.ui.panel.*;
import org.broad.igv.ui.util.IGVMouseInputAdapter;
import org.broad.igv.ui.util.UIUtilities;

import javax.swing.*;
import javax.swing.event.PopupMenuEvent;
import javax.swing.event.PopupMenuListener;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.geom.Rectangle2D;
import java.io.File;
import java.util.*;
import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.stream.Collectors;

/**
 * Window for displaying sashimi style junction plot
 * See http://genes.mit.edu/burgelab/miso/docs/sashimi.html
 * <p/>
 * User: jacob
 * Date: 2013-Jan-11
 */
public class SashimiPlot extends JFrame implements IGVEventObserver {

    private final SashimiContentPane sashimiContentPane;

    private final SelectableFeatureTrack featureTrack;

    private List<SpliceJunctionTrack> spliceJunctionTracks;

    private ReferenceFrame referenceFrame;

    private IGVEventBus eventBus;

    private static final List<Color> plotColors;

    private Map<Object, SashimiJunctionRenderer> junctionRendererMap;

    static {
        ColorPalette palette = ColorUtilities.getDefaultPalette();
        plotColors = Arrays.asList(palette.getColors());
    }

    public SashimiPlot(ReferenceFrame iframe, Collection<AlignmentTrack> alignmentTracks, FeatureTrack geneTrack) {

        // setContentPane(new SashimiContentPane());
        getGlassPane().setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
        getContentPane().setLayout(new BorderLayout());
        setBackground(Color.white);

        this.eventBus = new IGVEventBus();

        this.referenceFrame = new ReferenceFrame(iframe, eventBus);

        int height = IGV.hasInstance() ? IGV.getInstance().getMainFrame().getHeight() : 800;
        setSize(referenceFrame.getWidthInPixels(), height);

        //Add control elements to the top
        final JPanel controlPanel = generateControlPanel(this.referenceFrame);
        controlPanel.setAlignmentX(Component.CENTER_ALIGNMENT);
        getContentPane().add(controlPanel, BorderLayout.NORTH);


        JPanel sashimiPanel = new JPanel();
        sashimiPanel.setBackground(Color.WHITE);
        BoxLayout boxLayout = new BoxLayout(sashimiPanel, BoxLayout.Y_AXIS);
        sashimiPanel.setLayout(boxLayout);

        spliceJunctionTracks = new ArrayList<>(alignmentTracks.size());

        int colorInd = 0;

        junctionRendererMap = new HashMap<>();

        eventBus.subscribe(ViewChange.class, this);

        for (AlignmentTrack alignmentTrack : alignmentTracks) {

            AlignmentDataManager dataManager = alignmentTrack.getDataManager();

            SpliceJunctionTrack spliceJunctionTrack =
                    new SpliceJunctionTrack(alignmentTrack.getResourceLocator(), alignmentTrack.getName(), dataManager, alignmentTrack, SpliceJunctionTrack.StrandOption.COMBINE);

            // Override expand/collpase setting -- expanded sashimi plots make no sense
            spliceJunctionTrack.setDisplayMode(Track.DisplayMode.COLLAPSED);

            SashimiJunctionRenderer renderer = new SashimiJunctionRenderer();
            spliceJunctionTrack.setRenderer(renderer);
            junctionRendererMap.put(spliceJunctionTrack, renderer);

            Color color = plotColors.get(colorInd);
            colorInd = (colorInd + 1) % plotColors.size();
            spliceJunctionTrack.setColor(color);

            TrackComponent<SpliceJunctionTrack> trackComponent = new TrackComponent<>(referenceFrame, spliceJunctionTrack);
            trackComponent.originalFrame = iframe;

            initSpliceJunctionComponent(trackComponent, dataManager, dataManager.getCoverageTrack());

            sashimiPanel.add(trackComponent);
            spliceJunctionTracks.add(spliceJunctionTrack);

            spliceJunctionTrack.load(iframe);  // <= Must "load" tracks with frame of alignment track (actually just fetches from cache)
        }

        Axis axis = createAxis(referenceFrame);
        sashimiPanel.add(axis);

        featureTrack = new SelectableFeatureTrack(geneTrack);
        TrackComponent<SelectableFeatureTrack> geneComponent = new TrackComponent<>(referenceFrame, featureTrack);
        initGeneComponent(referenceFrame.getWidthInPixels(), geneComponent, featureTrack);

        JScrollPane scrollableGenePane = new JScrollPane(geneComponent);
        scrollableGenePane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        scrollableGenePane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED);

        sashimiContentPane = new SashimiContentPane(sashimiPanel, scrollableGenePane);
        sashimiContentPane.setDividerLocation(2 * height / 3);
        getContentPane().add(sashimiContentPane);

        validate();
    }

    private JPanel generateControlPanel(ReferenceFrame frame) {
        JPanel controlPanel = new JPanel();

        ZoomSliderPanel zoomSliderPanel = new ZoomSliderPanel(frame);
        zoomSliderPanel.setMinZoomLevel(frame.getZoom());

        zoomSliderPanel.addMouseListener(new IGVMouseInputAdapter() {
            @Override
            public void igvMouseClicked(MouseEvent e) {
                SashimiPlot.this.repaint();
            }
        });

        Dimension controlSize = new Dimension(200, 30);

        controlPanel.add(zoomSliderPanel);
        setFixedSize(zoomSliderPanel, controlSize);

        Dimension panelSize = controlSize;
        setFixedSize(controlPanel, panelSize);

        return controlPanel;
    }

    private static void setFixedSize(Component component, Dimension dimension) {
        component.setPreferredSize(dimension);
        component.setMinimumSize(dimension);
        component.setMaximumSize(dimension);
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
        RenderContext context = new RenderContext(geneComponent, null, referenceFrame, null);
        geneTrack.load(context.getReferenceFrame());

        GeneTrackMouseAdapter ad2 = new GeneTrackMouseAdapter(geneComponent);
        geneComponent.addMouseListener(ad2);
        geneComponent.addMouseMotionListener(ad2);
    }

    private void initSpliceJunctionComponent(TrackComponent<SpliceJunctionTrack> trackComponent, AlignmentDataManager dataManager, CoverageTrack coverageTrack) {
        JunctionTrackMouseAdapter ad1 = new JunctionTrackMouseAdapter(trackComponent);
        trackComponent.addMouseListener(ad1);
        trackComponent.addMouseMotionListener(ad1);

        getRenderer(trackComponent.track).setDataManager(dataManager);
        getRenderer(trackComponent.track).setCoverageTrack(coverageTrack);
        getRenderer(trackComponent.track).getCoverageTrack().rescale(trackComponent.originalFrame);
        getRenderer(trackComponent.track).setBackground(getBackground());  // <= this is neccessary
    }

    private SashimiJunctionRenderer getRenderer(SpliceJunctionTrack spliceJunctionTrack) {
        return junctionRendererMap.get(spliceJunctionTrack);
    }

    @Override
    public void receiveEvent(IGVEvent event) {
        List<CompletableFuture<Void>> futures = new ArrayList<>();

        if(!featureTrack.isReadyToPaint(referenceFrame)) {
            futures.add(CompletableFuture.runAsync(() -> featureTrack.load(referenceFrame), IGV.threadExecutor));
        }
        for (SpliceJunctionTrack t : spliceJunctionTracks) {
            if (!t.isReadyToPaint(referenceFrame)) {
                futures.add(CompletableFuture.runAsync(() -> t.load(referenceFrame), IGV.threadExecutor));
            }
        }

        if (futures.isEmpty()) {
            repaint();
        } else {
            CompletableFuture.allOf(futures.toArray(new CompletableFuture[0]))
                    .whenCompleteAsync((value, throwable) -> repaint(), UIUtilities::invokeOnEventThread);
        }
    }

    /**
     * Should consider using this elsewhere. Single component
     * which contains a single track
     */
    private static class TrackComponent<T extends Track> extends JComponent {

        private T track;
        private ReferenceFrame frame;
        private String toolTipText = null;
        public ReferenceFrame originalFrame;

        public TrackComponent(ReferenceFrame frame, T track) {
            this.frame = frame;
            this.track = track;
        }

        public void updateToolTipText(TrackClickEvent tce) {
            toolTipText = track.getValueStringAt(tce.getFrame().getChrName(), tce.getChromosomePosition(), tce.getMouseEvent().getX(), tce.getMouseEvent().getY(), tce.getFrame());
            toolTipText = toolTipText == null ? "" : "<html>" + toolTipText;
            setToolTipText(toolTipText);
        }

        @Override
        public void paintComponent(Graphics g) {
            super.paintComponent(g);
            Rectangle bounds = new Rectangle(getBounds());
            bounds.y = 0;
            RenderContext context = new RenderContext(this, (Graphics2D) g, frame, bounds);
            track.render(context, bounds);
        }

        @Override
        public Dimension getPreferredSize() {
            return new Dimension(100, track.getHeight());
        }
    }

    /**
     * Set the max coverage depth, which is a graphical scaling parameter for determining how
     * thick the junction arcs will be
     *
     * @param trackComponent
     * @param newMaxDepth
     */
    private void setMaxCoverageDepth(TrackComponent<SpliceJunctionTrack> trackComponent, int newMaxDepth) {
        getRenderer(trackComponent.track).setMaxDepth(newMaxDepth);
        repaint();
    }

    private class JunctionTrackMouseAdapter extends TrackComponentMouseAdapter<SpliceJunctionTrack> {

        JunctionTrackMouseAdapter(TrackComponent<SpliceJunctionTrack> trackComponent) {
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
            showCoverageData.setSelected(PreferencesManager.getPreferences().getAsBoolean(Constants.SASHIMI_SHOW_COVERAGE));
            showCoverageData.addActionListener(e16 -> {
                PreferencesManager.getPreferences().put(Constants.SASHIMI_SHOW_COVERAGE, showCoverageData.isSelected());
                SashimiPlot.this.repaint();
            });

            CoverageTrack covTrack = getRenderer(this.trackComponent.track).getCoverageTrack();
            covTrack.setWindowFunction(WindowFunction.max);
            JMenuItem setCoverageDataRange = CoverageTrack.addDataRangeItem(SashimiPlot.this, null, Arrays.asList(covTrack));
            setCoverageDataRange.setText("Set Exon Coverage Max");

            JMenuItem maxJunctionCoverageRange = new JMenuItem("Set Junction Coverage Max");
            maxJunctionCoverageRange.setToolTipText("The thickness of each line will be proportional to the coverage, up until this value");

            maxJunctionCoverageRange.addActionListener(e12 -> {
                String input = JOptionPane.showInputDialog(SashimiPlot.this, "Set Max Junction Coverage", getRenderer(trackComponent.track).getMaxDepth());
                if (input == null || input.length() == 0) return;
                try {
                    int newMaxDepth = Integer.parseInt(input);
                    setMaxCoverageDepth(JunctionTrackMouseAdapter.this.trackComponent, newMaxDepth);
                } catch (NumberFormatException ex) {
                    JOptionPane.showMessageDialog(SashimiPlot.this, input + " is not an integer");
                }
            });

            JMenuItem minJunctionCoverageItem = new JMenuItem("Set Junction Coverage Min");
            minJunctionCoverageItem.setToolTipText("Junctions below this threshold will be removed from view");
            minJunctionCoverageItem.addActionListener(e1 -> {
                int minCov = this.trackComponent.track.getMinJunctionCoverage();
                String input = JOptionPane.showInputDialog(SashimiPlot.this, "Set Minimum Junction Coverage", minCov);
                if (input == null || input.length() == 0) return;
                try {
                    int newMinJunctionCoverage = Integer.parseInt(input);
                    trackComponent.track.setMinJunctionCoverage(newMinJunctionCoverage);
                    trackComponent.track.clear();
                    trackComponent.repaint();
                } catch (NumberFormatException ex) {
                    JOptionPane.showMessageDialog(SashimiPlot.this, input + " is not an integer");
                }
            });


            JMenuItem colorItem = new JMenuItem("Set Color");
            colorItem.addActionListener(e15 -> {
                Color color = UIUtilities.showColorChooserDialog(
                        "Select Track Color", trackComponent.track.getColor());
                SashimiPlot.this.toFront();
                if (color == null) return;
                trackComponent.track.setColor(color);
                trackComponent.repaint();
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

            JRadioButtonMenuItem combineStrands = getStrandRadioButton("Combine Strands", SpliceJunctionTrack.StrandOption.COMBINE);
            combineStrands.setToolTipText("Combine junctions from both strands -- best for non-strand preserving libraries.");
            strandGroup.add(combineStrands);

            //  JRadioButtonMenuItem bothStrands = getStrandRadioButton("Both Strands", SpliceJunctionTrack.StrandOption.BOTH);
            //  strandGroup.add(bothStrands);
            //  menu.add(bothStrands);

            JRadioButtonMenuItem plusStrand = getStrandRadioButton("Forward Strand", SpliceJunctionTrack.StrandOption.FORWARD);
            plusStrand.setToolTipText("Show only junctions on the forward read strand  (of first-in-pair for paired reads)");
            strandGroup.add(plusStrand);

            JRadioButtonMenuItem minusStrand = getStrandRadioButton("Reverse Strand", SpliceJunctionTrack.StrandOption.REVERSE);
            plusStrand.setToolTipText("Show only junctions on the reverse read strand  (of first-in-pair for paired reads)");
            strandGroup.add(minusStrand);


            JMenuItem savePngImageItem = new JMenuItem("Save PNG Image...");
            savePngImageItem.addActionListener(e13 -> {
                File defaultFile = new File("Sashimi.png");
                IGV.getInstance().createSnapshot(SashimiPlot.this.sashimiContentPane, defaultFile);
            });

            JMenuItem saveSvgImageItem = new JMenuItem("Save SVG Image...");
            saveSvgImageItem.addActionListener(e14 -> {
                File defaultFile = new File("Sashimi.svg");
                IGV.getInstance().createSnapshot(SashimiPlot.this.sashimiContentPane, defaultFile);
            });

            // Coverage ranges -- these apply to current plot only
            menu.add(new JLabel("Junction Coverage Display"));
            menu.add(setCoverageDataRange);
            menu.add(minJunctionCoverageItem);
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

            menu.add(savePngImageItem);
            menu.add(saveSvgImageItem);

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

        private JRadioButtonMenuItem getStrandRadioButton(String label, final SpliceJunctionTrack.StrandOption option) {

            JRadioButtonMenuItem item = new JRadioButtonMenuItem(label);
            item.setSelected(SpliceJunctionTrack.getStrandOption() == option);
            item.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    SpliceJunctionTrack.setStrandOption(option);
                    for (SpliceJunctionTrack t : spliceJunctionTracks) {
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
            for (SpliceJunctionTrack spliceTrack : spliceJunctionTracks) {
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
            TrackMenuUtils.addDisplayModeItems(Arrays.asList(trackComponent.track), menu);
            menu.addPopupMenuListener(new RepaintPopupMenuListener(trackComponent));
            return menu;
        }
    }

    private abstract class TrackComponentMouseAdapter<T extends Track> extends IGVMouseInputAdapter {

        protected TrackComponent<T> trackComponent;
        protected PanTool currentTool;

        TrackComponentMouseAdapter(TrackComponent<T> trackComponent) {
            this.trackComponent = trackComponent;
            currentTool = new PanTool(null);
            currentTool.setReferenceFrame(this.trackComponent.frame);
        }


        @Override
        public void mouseDragged(MouseEvent e) {
            currentTool.mouseDragged(e);
            repaint();
        }

        @Override
        public void mouseReleased(MouseEvent e) {
            super.mouseReleased(e);
            if (e.isPopupTrigger()) {
                doPopupMenu(e);
            } else {
                currentTool.mouseReleased(e);
            }
        }

        @Override
        public void mousePressed(MouseEvent e) {
            super.mousePressed(e);
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
        public void igvMouseClicked(MouseEvent e) {
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
     */
    public static void openSashimiPlot() {

        FeatureTrack geneTrack = null;

        List<FeatureTrack> nonJunctionTracks = new ArrayList(IGV.getInstance().getFeatureTracks());
        Iterator iter = nonJunctionTracks.iterator();
        while (iter.hasNext()) {
            if (iter.next() instanceof SpliceJunctionTrack) iter.remove();
        }

        if (nonJunctionTracks.size() == 1) {
            geneTrack = nonJunctionTracks.get(0);
        } else {
            FeatureTrackSelectionDialog dlg = new FeatureTrackSelectionDialog(IGV.getInstance().getMainFrame(), nonJunctionTracks);
            dlg.setTitle("Select Gene Track");
            dlg.setVisible(true);
            if (dlg.getIsCancelled()) {
                return;
            }
            geneTrack = dlg.getSelectedTrack();
        }


        Collection<AlignmentTrack> alignmentTracks = IGV.getInstance().getAlignmentTracks().stream()
                .filter(t -> t.getExperimentType() == AlignmentTrack.ExperimentType.RNA).collect(Collectors.toList());
        ;
        if (alignmentTracks.size() > 1) {
            TrackSelectionDialog<AlignmentTrack> alDlg =
                    new TrackSelectionDialog<AlignmentTrack>(IGV.getInstance().getMainFrame(), TrackSelectionDialog.SelectionMode.MULTIPLE, alignmentTracks);
            alDlg.setTitle("Select Alignment Tracks");
            alDlg.setVisible(true);
            if (alDlg.getIsCancelled()) return;

            alignmentTracks = alDlg.getSelectedTracks();
        }

        SashimiPlot sashimiPlot = new SashimiPlot(FrameManager.getDefaultFrame(), alignmentTracks, geneTrack);
        //sashimiPlot.setShapeType(shapeType);
        sashimiPlot.setVisible(true);


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
            UIUtilities.invokeOnEventThread(() -> {
                component.revalidate();
                component.repaint();
            });

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
            RenderContext context = new RenderContext(this, (Graphics2D) g, frame, visibleRect);
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
