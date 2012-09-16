package org.broad.igv.hic.track;

import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.hic.HiC;
import org.broad.igv.hic.HiCRenderContext;
import org.broad.igv.renderer.GraphicUtils;
import org.broad.igv.track.*;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.panel.IGVPopupMenu;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.Pair;
import org.broad.igv.util.ResourceLocator;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;

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
    HiCTrack eigenvectorTrack;
    Collection<Pair<Rectangle, Track>> trackRectangles;

    public TrackPanel(HiC hiC) {
        this.hic = hiC;
        setAutoscrolls(true);
        trackRectangles = new ArrayList<Pair<Rectangle, Track>>();
        addMouseListener();
    }

    private void addMouseListener() {
        this.addMouseListener(new MouseAdapter() {
            @Override
            public void mouseReleased(MouseEvent mouseEvent) {
                if (mouseEvent.isPopupTrigger()) {
                    handlePopupEvent(mouseEvent);
                }
            }

            @Override
            public void mouseClicked(MouseEvent mouseEvent) {
                if (mouseEvent.isPopupTrigger()) {
                    handlePopupEvent(mouseEvent);
                }
            }

            @Override
            public void mousePressed(MouseEvent mouseEvent) {
                if (mouseEvent.isPopupTrigger()) {
                    handlePopupEvent(mouseEvent);
                }
             }

            private void handlePopupEvent(MouseEvent mouseEvent) {
                for (Pair<Rectangle, Track> p : trackRectangles) {
                    Rectangle r = p.getFirst();
                    if (r.contains(mouseEvent.getPoint())) {
//                        Collection<Track> selectedTracks = Arrays.asList(p.getSecond());
//                        TrackClickEvent te = new TrackClickEvent(mouseEvent, null);
//                        IGVPopupMenu menu = TrackMenuUtils.getPopupMenu(selectedTracks, "", te);
//                        menu.show(mouseEvent.getComponent(), mouseEvent.getX(), mouseEvent.getY());
                    }
                }
            }

        });
    }

    public void setEigenvectorTrack(HiCTrack eigenvectorTrack) {
        this.eigenvectorTrack = eigenvectorTrack;
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
        if (eigenvectorTrack != null) {
            h += eigenvectorTrack.getHeight();
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

    protected void paintComponent(Graphics graphics) {

        trackRectangles.clear();
        java.util.List<Track> tracks = new ArrayList<Track>(HiCTrackManager.getLoadedTracks());
        if ((tracks == null || tracks.isEmpty()) && eigenvectorTrack == null) {
            return;
        }


        Rectangle rect = getBounds();
        int y = rect.y;

        for (Track track : tracks) {
            if (track.getHeight() > 0) {
                int h = track.getHeight();
                Rectangle trackRectangle = new Rectangle(rect.x, y, rect.width, h);

                if (hic.xContext != null) {
                    RenderContext context = new HiCRenderContext(hic.xContext, this, (Graphics2D) graphics, trackRectangle, genome);
                    track.render(context, trackRectangle);
                    renderName(track, trackRectangle, graphics);
                    y += h;

                    trackRectangles.add(new Pair(trackRectangle, track));
                }


            }
        }
        if (eigenvectorTrack != null) {
            int h = rect.y + rect.height - y;
            Rectangle trackRectangle = new Rectangle(rect.x, y, rect.width, h);
            eigenvectorTrack.render((Graphics2D) graphics, hic.xContext, trackRectangle);
            trackRectangles.add(new Pair(trackRectangle, eigenvectorTrack));

        }

    }

    private void renderName(Track track, Rectangle rect, Graphics graphics) {
        Font font = FontManager.getFont(8);
        graphics.setFont(font);
        graphics.setColor(track.getColor());
        GraphicUtils.drawRightJustifiedText(track.getName(), rect.x + rect.width - 10, rect.y + 15, graphics);
    }


}
