package org.broad.igv.hic;

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


    Track track;
    ReferenceFrame referenceFrame;

    public TrackPanel() {
        // Hack
        referenceFrame = new ReferenceFrame("HiC");
     }

    public void setTrack(Track track) {
        this.track = track;
    }

    //   public RenderContext(String genomeId, JComponent panel, Graphics2D graphics, ReferenceFrame referenceFrame, Rectangle visibleRect) {

    protected void paintComponent(Graphics graphics) {

        track.setHeight(getHeight());

        RenderContext context = new RenderContextImpl("hg18", this, (Graphics2D) graphics, referenceFrame,
                getVisibleRect());
        track.render(context, getBounds());


    }

}
