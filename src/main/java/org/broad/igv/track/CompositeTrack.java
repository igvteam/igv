package org.broad.igv.track;

import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.ResourceLocator;

import java.awt.*;
import java.util.List;
import java.util.stream.Collectors;


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
        for(Track t : tracks) {
            t.load(frame);
        }
    }

    @Override
    public void render(RenderContext context, Rectangle rect) {

        Rectangle r = new Rectangle(rect);
        for(Track t : tracks) {
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
}
