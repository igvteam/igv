package org.broad.igv.track;

import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.ResourceLocator;

import java.awt.*;
import java.util.List;
import java.util.stream.Collectors;


/**
 * Container for a collection of tracks.  Each subtrack is drawn at full height (content height), with no scrollbar.
 */

public class CompositeTrack extends AbstractTrack {

    List<Track> tracks;

    public CompositeTrack(ResourceLocator locator, List<Track> tracks) {
        super(locator);
        this.tracks = tracks;
    }

    @Override
    public boolean isReadyToPaint(ReferenceFrame frame) {
        return tracks.stream().allMatch(t -> t.isReadyToPaint(frame));
    }

    @Override
    public void load(ReferenceFrame frame) {
        for (Track t : tracks) {
            t.load(frame);
        }
    }

    @Override
    public void render(RenderContext context, Rectangle rect) {
        Rectangle r = new Rectangle(rect);
        for (Track t : tracks) {
            r.height = t.getContentHeight();
            t.render(context, r);
            r.y += r.height;
        }
    }

    @Override
    public int getContentHeight() {
        int h = tracks.stream().collect(Collectors.summingInt(t -> t.getContentHeight()));
        return h;
    }

    @Override
    public int getDefaultHeight() {
        int h = tracks.stream().collect(Collectors.summingInt(t -> t.getDefaultHeight()));
        return Math.min(300, h);
    }

    @Override
    public boolean handleDataClick(TrackClickEvent te) {
        Track t = getTrackAtY(te.getMouseEvent().getY());
        return t == null ? false : t.handleDataClick(te);
    }

    @Override
    public String getValueStringAt(String chr, double position, int mouseX, int mouseY, ReferenceFrame frame) {
        Track t = getTrackAtY(mouseY);
        return t == null ? "" : t.getValueStringAt(chr, position, mouseX, mouseY, frame);
    }

    protected Track getTrackAtY(int y) {
        int trackY = 0;
        for (Track t : tracks) {
            int h = t.getContentHeight();
            if(y > trackY && y <= trackY + h) {
                return t;
            }
            trackY += h;
        }
        return null;
    }
}
