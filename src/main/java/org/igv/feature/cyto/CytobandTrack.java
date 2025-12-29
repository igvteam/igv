package org.igv.feature.cyto;

import org.igv.Globals;
import org.igv.feature.CytoBandFileParser;
import org.igv.feature.Cytoband;
import org.igv.feature.genome.Genome;
import org.igv.renderer.CytobandRenderer;
import org.igv.track.AbstractTrack;
import org.igv.track.RenderContext;
import org.igv.track.TrackClickEvent;
import org.igv.ui.panel.IGVPopupMenu;
import org.igv.ui.panel.ReferenceFrame;
import org.igv.util.ResourceLocator;

import java.awt.*;
import java.io.BufferedReader;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class CytobandTrack extends AbstractTrack {

    Map<String, List<Cytoband>> chrCytoMap;
    CytobandRenderer renderer;

    public CytobandTrack(ResourceLocator locator, BufferedReader reader, Genome genome) {
        super(locator);
        this.setHeight(30);
        renderer = new CytobandRenderer(Globals.isDarkMode());
        load(reader, genome);
    }


    protected void load(BufferedReader reader, Genome genome) {
        Map<String, List<Cytoband>> map = CytoBandFileParser.loadData(reader);
        chrCytoMap = new HashMap<>();
        for (Map.Entry<String, List<Cytoband>> entry : map.entrySet()) {
            String chr = genome == null ? entry.getKey() : genome.getCanonicalChrName(entry.getKey());
            chrCytoMap.put(chr, entry.getValue());
        }
    }

    @Override
    public boolean isReadyToPaint(ReferenceFrame frame) {
        return chrCytoMap != null;
    }

    @Override
    public void load(ReferenceFrame frame) {
        // Nothing to do
    }

    @Override
    public void render(RenderContext context, Rectangle rect) {
        List<Cytoband> cytobands = chrCytoMap.get(context.getChr());
        if (cytobands != null) {
            Rectangle trackRect = context.getVisibleRect().intersection(rect);
            renderer.drawTrack(cytobands, context.getGraphics(), trackRect, context.getReferenceFrame());
        }
    }

    @Override
    public IGVPopupMenu getPopupMenu(final TrackClickEvent te) {
        return new IGVPopupMenu();
    }

}
