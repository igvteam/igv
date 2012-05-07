package org.broad.igv.hic;

import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.RenderContextImpl;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.Track;
import org.broad.igv.ui.panel.ReferenceFrame;

import javax.swing.*;
import java.awt.*;

/**
 * Created by IntelliJ IDEA.
 * User: neva
 * Date: 4/3/12
 * Time: 4:08 PM
 * To change this template use File | Settings | File Templates.
 */
public class TrackPanel extends JPanel {

    HiC hic;
    Track track;
    Genome genome;

    public TrackPanel(HiC hiC) {
        this.hic = hiC;
    }

    public void setTrack(Track track) {
        this.track = track;
    }

    //   public RenderContext(String genomeId, JComponent panel, Graphics2D graphics, ReferenceFrame referenceFrame, Rectangle visibleRect) {

    protected void paintComponent(Graphics graphics) {

        if (track == null) {
            return;
        }

        track.setHeight(getHeight());  // <= TODO move to setBounds
        if (hic.xContext != null) {
            RenderContext context = new HiCRenderContext(hic.xContext, this, (Graphics2D) graphics, getVisibleRect(), genome);
            track.render(context, getBounds());
        }

    }

    public void setGenome(Genome genome) {
        this.genome = genome;
    }
}
