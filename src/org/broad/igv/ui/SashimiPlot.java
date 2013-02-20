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

import com.google.common.eventbus.Subscribe;
import org.broad.igv.feature.IExon;
import org.broad.igv.renderer.SashimiJunctionRenderer;
import org.broad.igv.sam.AlignmentDataManager;
import org.broad.igv.sam.SpliceJunctionFinderTrack;
import org.broad.igv.track.*;
import org.broad.igv.ui.event.ViewChange;
import org.broad.igv.ui.panel.*;

import javax.swing.*;
import javax.swing.event.PopupMenuEvent;
import javax.swing.event.PopupMenuListener;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.*;
import java.util.List;

/** Window for displaying sashimi style junction plot
 * See http://genes.mit.edu/burgelab/miso/docs/sashimi.html
 *
 * User: jacob
 * Date: 2013-Jan-11
 */
public class SashimiPlot extends JFrame{

    private SpliceJunctionFinderTrack spliceJunctionTrack;
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


    public SashimiPlot(ReferenceFrame iframe, SpliceJunctionFinderTrack track, FeatureTrack geneTrack){
        this.frame = new ReferenceFrame(iframe);
        this.frame.getEventBus().register(this);

        minOrigin = this.frame.getOrigin();
        maxEnd = this.frame.getEnd();
        int minZoom = this.frame.getZoom();

        initSize(frame.getWidthInPixels());

        BoxLayout boxLayout = new BoxLayout(getContentPane(), BoxLayout.Y_AXIS);
        getContentPane().setLayout(boxLayout);

        spliceJunctionTrack = new SpliceJunctionFinderTrack(track.getResourceLocator(), track.getName(),
               track.getDataManager(), track.getGenome());

        spliceJunctionTrack.setRendererClass(SashimiJunctionRenderer.class);
        TrackComponent<SpliceJunctionFinderTrack> trackComponent = new TrackComponent<SpliceJunctionFinderTrack>(frame, spliceJunctionTrack);

        setDataManager(track.getDataManager());
        getRenderer().setBackground(getBackground());

        SelectableFeatureTrack geneTrackClone = new SelectableFeatureTrack(geneTrack);
        geneTrackClone.setDisplayMode(Track.DisplayMode.SQUISHED);
        TrackComponent<SelectableFeatureTrack> geneComponent = new TrackComponent<SelectableFeatureTrack>(frame, geneTrackClone);
        //Hacky way of clearing packed features
        geneTrackClone.setVisibilityWindow(geneTrackClone.getVisibilityWindow());

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


        getContentPane().add(trackComponent);
        getContentPane().add(geneComponent);

        initMouseAdapters(trackComponent, geneComponent);

        validate();
    }

    private void setDataManager(AlignmentDataManager dataManager) {
        if(!dataManager.isShowSpliceJunctions()){
            dataManager.setShowSpliceJunctions(true);
            dataManager.clear();
        }
        getRenderer().setDataManager(dataManager);
    }

    private void initSize(int width) {
        setSize(width, 500);
    }

    public void setShapeType(SashimiJunctionRenderer.ShapeType shapeType) {
        ((SashimiJunctionRenderer) spliceJunctionTrack.getRenderer()).setShapeType(shapeType);
        repaint();
    }

    private void initMouseAdapters(TrackComponent<SpliceJunctionFinderTrack> trackComponent, TrackComponent<SelectableFeatureTrack> geneComponent) {
        JunctionTrackMouseAdapter ad1 = new JunctionTrackMouseAdapter(trackComponent);
        trackComponent.addMouseListener(ad1);
        trackComponent.addMouseMotionListener(ad1);

        GeneTrackMouseAdapter ad2 = new GeneTrackMouseAdapter(geneComponent);
        geneComponent.addMouseListener(ad2);
        geneComponent.addMouseMotionListener(ad2);
    }

    private SashimiJunctionRenderer getRenderer(){
        return (SashimiJunctionRenderer) spliceJunctionTrack.getRenderer();
    }

    @Subscribe
    public void respondViewResult(ViewChange.Result e) {
        repaint();
    }

    /**
     * Should consider using this elsewhere. Single component
     * which contains a single track
     */
    private static class TrackComponent<T extends Track> extends JComponent{

        private T track;
        private ReferenceFrame frame;

        public TrackComponent(ReferenceFrame frame, T track){
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

    private class JunctionTrackMouseAdapter extends TrackComponentMouseAdapter<SpliceJunctionFinderTrack>{

        JunctionTrackMouseAdapter(TrackComponent<SpliceJunctionFinderTrack> trackComponent){
            super(trackComponent);
        }

        @Override
        protected void handleDataClick(MouseEvent e) {
            //Show data of some sort?
        }

        @Override
        protected IGVPopupMenu getPopupMenu(MouseEvent e) {
            IGVPopupMenu menu = new IGVPopupMenu();

            JMenuItem maxDepthItem = new JMenuItem("Set Max Depth");
            maxDepthItem.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    String input = JOptionPane.showInputDialog("Set Maximum Depth", getRenderer().getMaxDepth());
                    try {
                        int newMaxDepth = Integer.parseInt(input);
                        getRenderer().setMaxDepth(newMaxDepth);
                        repaint();
                    } catch (NumberFormatException ex) {
                        JOptionPane.showMessageDialog(SashimiPlot.this, input + " is not an integer");
                    }
                }
            });

            menu.add(maxDepthItem);

            for(JMenuItem item: getRenderMenuItems(null, SashimiPlot.this)){
                menu.add(item);
            }
            return menu;
        }
    }

    private class GeneTrackMouseAdapter extends TrackComponentMouseAdapter<SelectableFeatureTrack>{

        GeneTrackMouseAdapter(TrackComponent<SelectableFeatureTrack> trackComponent){
            super(trackComponent);
        }

        @Override
        protected void handleDataClick(MouseEvent e) {
            trackComponent.track.handleDataClick(createTrackClickEvent(e));
            Set<IExon> selectedExon = trackComponent.track.getSelectedExons();
            getRenderer().setSelectedExons(selectedExon);
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

    private abstract class TrackComponentMouseAdapter<T extends Track> extends MouseAdapter{

        protected TrackComponent<T> trackComponent;
        protected PanTool currentTool;

        TrackComponentMouseAdapter(TrackComponent<T> trackComponent){
            this.trackComponent = trackComponent;
            currentTool = new PanTool(null);
            currentTool.setReferenceFrame(this.trackComponent.frame);
        }


        @Override
        public void mouseDragged(MouseEvent e) {
            double diff = e.getX() - currentTool.getLastMousePoint().getX();
            // diff > 0 means moving mouse to the right, which drags the frame towards the negative direction
            boolean hitBounds = SashimiPlot.this.frame.getOrigin() <= minOrigin && diff > 0;
            hitBounds |= SashimiPlot.this.frame.getEnd() >= maxEnd && diff < 0;
            if(!hitBounds){
                currentTool.mouseDragged(e);
                repaint();
            }
        }

        @Override
        public void mouseReleased(MouseEvent e) {
            currentTool.mouseReleased(e);
        }

        @Override
        public void mousePressed(MouseEvent e) {
            if(e.isPopupTrigger()){
                doPopupMenu(e);
            }else{
                currentTool.mousePressed(e);
                super.mousePressed(e);
            }

        }

        protected void doPopupMenu(MouseEvent e){
            IGVPopupMenu menu = getPopupMenu(e);
            if(menu != null) menu.show(trackComponent, e.getX(), e.getY());
        }

        protected TrackClickEvent createTrackClickEvent(MouseEvent e) {
            return new TrackClickEvent(e, trackComponent.frame);
        }

        @Override
        public void mouseClicked(MouseEvent e) {
            if(e.isPopupTrigger()){
                doPopupMenu(e);
                return;
            }

            currentTool.mouseClicked(e);
            handleDataClick(e);
        }

        /**
         * Essentially left click
         * @param e
         */
        protected abstract void handleDataClick(MouseEvent e);

        /**
         * Essentially right click
         * @param e
         * @return
         */
        protected abstract IGVPopupMenu getPopupMenu(MouseEvent e);
    }

    public static List<JMenuItem> getRenderMenuItems(final SpliceJunctionFinderTrack junctionFinderTrack, final SashimiPlot dialog){

        //JCheckBoxMenuItem setSplice = new JCheckBoxMenuItem("Splice Junction");
        //setSplice.addActionListener(getChangeClassListener(setSplice, SpliceJunctionRenderer.class, null));
        //setSplice.setSelected(SpliceJunctionFinderTrack.this.getRenderer().getClass().equals(SpliceJunctionRenderer.class));

        //setRenderingStyle.add(setSplice);

        List<JMenuItem> menuItemList = new ArrayList<JMenuItem>(3);

        Map<String, SashimiJunctionRenderer.ShapeType> renderTypes = new LinkedHashMap<String, SashimiJunctionRenderer.ShapeType>(3);
        renderTypes.put("Sashimi Ellipse", SashimiJunctionRenderer.ShapeType.ELLIPSE);
        renderTypes.put("Sashimi Circle", SashimiJunctionRenderer.ShapeType.CIRCLE);
        renderTypes.put("Sashimi Text", SashimiJunctionRenderer.ShapeType.TEXT);
        for(final Map.Entry<String, SashimiJunctionRenderer.ShapeType> entry: renderTypes.entrySet()){
            JMenuItem tmpSashimi = new JCheckBoxMenuItem(entry.getKey());
            tmpSashimi.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    getSashimiPlot(junctionFinderTrack, dialog, entry.getValue());
                }
            });
            menuItemList.add(tmpSashimi);
        }
        return menuItemList;
    }

    /**
     * Show SashimiPlot window, or change settings of {@code currentWindow}
     * @param sashimiPlot
     * @param shapeType
     */
    public static void getSashimiPlot(SpliceJunctionFinderTrack junctionFinderTrack, SashimiPlot sashimiPlot, SashimiJunctionRenderer.ShapeType shapeType) {
        if(sashimiPlot == null){
            FeatureTrackSelectionDialog dlg = new FeatureTrackSelectionDialog(IGV.getMainFrame());
            dlg.setVisible(true);
            if (dlg.getIsCancelled()) return;

            FeatureTrack geneTrack = dlg.getSelectedTrack();

            sashimiPlot = new SashimiPlot(FrameManager.getDefaultFrame(), junctionFinderTrack, geneTrack);
            sashimiPlot.setShapeType(shapeType);
            sashimiPlot.setVisible(true);
        }else{
            sashimiPlot.setShapeType(shapeType);
        }


    }

    private static class RepaintPopupMenuListener implements PopupMenuListener{

        Component component;

        RepaintPopupMenuListener(Component component){
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
}
