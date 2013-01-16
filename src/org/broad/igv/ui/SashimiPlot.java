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

import org.broad.igv.renderer.SashimiJunctionRenderer;
import org.broad.igv.sam.SpliceJunctionFinderTrack;
import org.broad.igv.track.*;
import org.broad.igv.ui.panel.FeatureTrackSelectionDialog;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.IGVPopupMenu;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.tribble.Feature;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/** Window for displaying sashimi style junction plot
 * See http://genes.mit.edu/burgelab/miso/docs/sashimi.html
 *
 * User: jacob
 * Date: 2013-Jan-11
 */
public class SashimiPlot extends JFrame{

    private SpliceJunctionFinderTrack spliceJunctionTrack;

    public SashimiPlot(ReferenceFrame frame, SpliceJunctionFinderTrack track, FeatureTrack geneTrack){
        initSize(frame.getWidthInPixels());
        BoxLayout boxLayout = new BoxLayout(getContentPane(), BoxLayout.Y_AXIS);
        getContentPane().setLayout(boxLayout);

        spliceJunctionTrack = new SpliceJunctionFinderTrack(track.getResourceLocator(), track.getName(),
               track.getDataManager(), track.getGenome());

        spliceJunctionTrack.setRendererClass(SashimiJunctionRenderer.class);
        TrackComponent<SpliceJunctionFinderTrack> trackComponent = new TrackComponent<SpliceJunctionFinderTrack>(frame, spliceJunctionTrack);

        FeatureTrack geneTrackClone = new FeatureTrack(geneTrack);
        TrackComponent<FeatureTrack> geneComponent = new TrackComponent<FeatureTrack>(frame, geneTrackClone);

        getContentPane().add(trackComponent);
        getContentPane().add(geneComponent);

        initMouseAdapters(trackComponent, geneComponent);

        validate();
    }

    private void initSize(int width) {
        setSize(width, 500);
    }

    public void setShapeType(SashimiJunctionRenderer.ShapeType shapeType) {
        ((SashimiJunctionRenderer) spliceJunctionTrack.getRenderer()).setShapeType(shapeType);
        repaint();
    }

    private void initMouseAdapters(TrackComponent<SpliceJunctionFinderTrack> trackComponent, TrackComponent<FeatureTrack> geneComponent) {
        trackComponent.addMouseListener(new JunctionTrackMouseAdapter(trackComponent));
        geneComponent.addMouseListener(new GeneTrackMouseAdapter(geneComponent));
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
            for(JMenuItem item: getRenderMenuItems(null, SashimiPlot.this)){
                menu.add(item);
            }
            return menu;
        }
    }

    private class GeneTrackMouseAdapter extends TrackComponentMouseAdapter<FeatureTrack>{

        GeneTrackMouseAdapter(TrackComponent<FeatureTrack> trackComponent){
            super(trackComponent);
        }

        @Override
        protected void handleDataClick(MouseEvent e) {
            trackComponent.track.handleDataClick(createTrackClickEvent(e));
            Feature selectedExon = trackComponent.track.getSelectedExon();
            ((SashimiJunctionRenderer) spliceJunctionTrack.getRenderer()).setSelectedExon(selectedExon);
            repaint();
        }

        @Override
        protected IGVPopupMenu getPopupMenu(MouseEvent e) {
            return null; //TODO
        }
    }

    private abstract static class TrackComponentMouseAdapter<T extends Track> extends MouseAdapter{

        protected TrackComponent<T> trackComponent;

        TrackComponentMouseAdapter(TrackComponent<T> trackComponent){
            this.trackComponent = trackComponent;
        }


        @Override
        public void mousePressed(MouseEvent e) {
            if(e.isPopupTrigger()){
                doPopupMenu(e);
            }else{
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
}
