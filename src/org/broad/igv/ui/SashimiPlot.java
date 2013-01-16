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

import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.renderer.SashimiJunctionRenderer;
import org.broad.igv.sam.SpliceJunctionFinderTrack;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.RenderContextImpl;
import org.broad.igv.track.Track;
import org.broad.igv.ui.panel.ReferenceFrame;

import javax.swing.*;
import java.awt.*;

/** Window for displaying sashimi style junction plot
 * See http://genes.mit.edu/burgelab/miso/docs/sashimi.html
 *
 * User: jacob
 * Date: 2013-Jan-11
 */
public class SashimiPlot extends JFrame{

    private SpliceJunctionFinderTrack spliceJunctionTrack;

    public SashimiPlot(ReferenceFrame frame, SpliceJunctionFinderTrack track){
        initSize(frame.getWidthInPixels());
        BoxLayout boxLayout = new BoxLayout(getContentPane(), BoxLayout.Y_AXIS);
        getContentPane().setLayout(boxLayout);

        spliceJunctionTrack = new SpliceJunctionFinderTrack(track.getResourceLocator(), track.getName(),
               track.getDataManager(), track.getGenome());

        spliceJunctionTrack.setRendererClass(SashimiJunctionRenderer.class);
        TrackComponent trackComponent = new TrackComponent(frame, spliceJunctionTrack);

        FeatureTrack geneTrack = GenomeManager.getInstance().getCurrentGenome().getGeneTrack();
        FeatureTrack geneTrackClone = new FeatureTrack(geneTrack);
        TrackComponent geneComponent = new TrackComponent(frame, geneTrackClone);

        getContentPane().add(trackComponent);
        getContentPane().add(geneComponent);

        validate();
    }

    private void initSize(int width) {
        setSize(width, 500);
    }

    public void setShapeType(SashimiJunctionRenderer.ShapeType shapeType) {
        ((SashimiJunctionRenderer) spliceJunctionTrack.getRenderer()).setShapeType(shapeType);
    }


    /**
     * Should consider using this elsewhere. Single component
     * which contains a single track
     */
    private static class TrackComponent extends JComponent{

        private Track track;
        private ReferenceFrame frame;

        public TrackComponent(ReferenceFrame frame, Track track){
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
}
