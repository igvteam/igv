package org.broad.igv.hic.track;

import org.broad.igv.feature.genome.Genome;
import org.broad.igv.hic.HiC;
import org.broad.igv.hic.MainWindow;
import org.broad.igv.renderer.GraphicUtils;
import org.broad.igv.ui.FontManager;
import org.broad.igv.util.Pair;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.ArrayList;
import java.util.Collection;

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
    Collection<Pair<Rectangle, HiCTrack>> trackRectangles;

    public TrackPanel(HiC hiC) {
        this.hic = hiC;
        setAutoscrolls(true);
        trackRectangles = new ArrayList<Pair<Rectangle, HiCTrack>>();
        setBackground(new Color(238,238,238));
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
                for (Pair<Rectangle, HiCTrack> p : trackRectangles) {
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
        for (HiCTrack t : HiCTrackManager.getLoadedTracks()) {
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

        ((Graphics2D) graphics).setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

        trackRectangles.clear();
        java.util.List<HiCTrack> tracks = new ArrayList<HiCTrack>(HiCTrackManager.getLoadedTracks());
        if ((tracks == null || tracks.isEmpty()) && eigenvectorTrack == null) {
            return;
        }


        Rectangle rect = getBounds();
        graphics.setColor(getBackground());
        graphics.fillRect(rect.x, rect.y, rect.width, rect.height);
        int y = rect.y;

        for (HiCTrack hicTrack : tracks) {
            if (hicTrack.getHeight() > 0) {
                int h = hicTrack.getHeight();
                Rectangle trackRectangle = new Rectangle(rect.x, y, rect.width, h);

                if (hic.xContext != null) {
                    //  RenderContext context = new HiCRenderContext(hic.xContext, this, (Graphics2D) graphics, trackRectangle, genome);

                    hicTrack.render((Graphics2D) graphics, hic.xContext, trackRectangle);
                    renderName(hicTrack.getName(), hicTrack.getColor(), trackRectangle, graphics);
                    y += h;

                    trackRectangles.add(new Pair(trackRectangle, hicTrack));
                }


            }
        }
        if (eigenvectorTrack != null) {
            int h = rect.y + rect.height - y;
            Rectangle trackRectangle = new Rectangle(rect.x, y, rect.width, h);
            eigenvectorTrack.render((Graphics2D) graphics, hic.xContext, trackRectangle);

            renderName("Eigenvector", EigenvectorTrack.COLOR, trackRectangle, graphics);

            trackRectangles.add(new Pair(trackRectangle, eigenvectorTrack));


        }

        int cursorX = hic.getCursorX();
        if(cursorX > 0) {
            graphics.setColor(MainWindow.RULER_LINE_COLOR);
            graphics.drawLine(cursorX, 0, cursorX, getHeight());
        }

    }

    private void renderName(String name, Color color, Rectangle rect, Graphics graphics) {
        Font font = FontManager.getFont(8);
        graphics.setFont(font);
        graphics.setColor(color);
        GraphicUtils.drawRightJustifiedText(name, rect.x + rect.width - 10, rect.y + 15, graphics);
    }


}
