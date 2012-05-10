package org.broad.igv.hic.track;

import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.hic.HiC;
import org.broad.igv.hic.HiCRenderContext;
import org.broad.igv.track.RenderContextImpl;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackLoader;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.ResourceLocator;

import javax.swing.*;
import java.awt.*;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: neva
 * Date: 4/3/12
 * Time: 4:08 PM
 * To change this template use File | Settings | File Templates.
 */
public class TrackPanel extends JPanel {

    HiC hic;
    Genome genome;

    public TrackPanel(HiC hiC) {
        this.hic = hiC;
        setAutoscrolls(true);
    }


    /**
     * Returns the current height of this component.
     * This method is preferable to writing
     * <code>component.getBounds().height</code>, or
     * <code>component.getSize().height</code> because it doesn't cause any
     * heap allocations.
     *
     * @return the current height of this component
     */
    @Override
    public int getHeight() {
        int h = 0;
        for (Track t : HiCTrackManager.getLoadedTracks()) {
            h += t.getHeight();
        }
        return h;
    }

    /**
     * If the <code>preferredSize</code> has been set to a
     * non-<code>null</code> value just returns it.
     * If the UI delegate's <code>getPreferredSize</code>
     * method returns a non <code>null</code> value then return that;
     * otherwise defer to the component's layout manager.
     *
     * @return the value of the <code>preferredSize</code> property
     * @see #setPreferredSize
     * @see javax.swing.plaf.ComponentUI
     */
    @Override
    public Dimension getPreferredSize() {
        return new Dimension(500, getHeight());
    }

    //   public RenderContext(String genomeId, JComponent panel, Graphics2D graphics, ReferenceFrame referenceFrame, Rectangle visibleRect) {

    protected void paintComponent(Graphics graphics) {

        java.util.List<Track> tracks = new ArrayList<Track>(HiCTrackManager.getLoadedTracks());
        if (tracks == null || tracks.isEmpty()) {
            return;
        }


        Rectangle rect = getBounds();

        for (Track track : tracks) {
            if (track.getHeight() > 0) {
                rect.height = track.getHeight();
                if (hic.xContext != null) {
                    RenderContext context = new HiCRenderContext(hic.xContext, this, (Graphics2D) graphics, rect, genome);
                    track.render(context, rect);
                }
                rect.y += rect.height;
            }
        }

    }

    public void setGenome(Genome genome) {
        this.genome = genome;
    }

    private void test() {
        try {
            String genomePath = "/Users/jrobinso/igv/genomes/hg19.genome";
            Genome genome = GenomeManager.getInstance().loadGenome(genomePath, null);

            String testURL = "http://www.broadinstitute.org/igvdata/encode/hg19/uwDnase/wgEncodeUwDnaseGm06990RawRep1.bigWig";
            java.util.List<Track> tracks = new ArrayList();
            (new TrackLoader()).loadBWFile(new ResourceLocator(testURL), tracks, genome);

            // eigenvectorTrack = new EigenvectorTrack("eigen", "Eigenvectors");
            Track track = tracks.get(0);
            track.setHeight(40);
            HiCTrackManager.addTrack(tracks.get(0));

        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }
}
