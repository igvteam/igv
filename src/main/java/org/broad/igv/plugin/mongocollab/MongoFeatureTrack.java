/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broad.igv.plugin.mongocollab;

import com.mongodb.DBCollection;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackClickEvent;
import org.broad.igv.track.TrackMenuUtils;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.action.SearchCommand;
import org.broad.igv.ui.panel.IGVPopupMenu;
import org.broad.igv.ui.panel.ReferenceFrame;
import htsjdk.tribble.Feature;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Arrays;
import java.util.List;

/**
 * @author jacob
 * @date 2013-Sep-27
 */
public class MongoFeatureTrack extends FeatureTrack {

    public MongoFeatureTrack(String id, String name, MongoFeatureSource source) {
        super(id, name, source);
    }

    @Override
    public void dispose() {
        super.dispose();
        SearchCommand.unregisterNamedFeatureSearcher((MongoFeatureSource) source);
    }

    /**
     * Return a menu entry for adding a new feature, or editing a feature the user clicked on.
     *
     * @param te
     * @param selFeat Selected feature. {@code null} indicates adding a new feature
     * @return
     */
    private JMenuItem createEditAnnotMenuEntry(TrackClickEvent te, final DBFeature.IGVFeat selFeat) {

        ReferenceFrame frame = te.getFrame();
        boolean hasFrame = frame != null;
        if (!hasFrame) return null;


        final boolean editing = selFeat != null;

        final String chr = editing ? selFeat.getChr() : frame.getChrName();
        final int start = editing ? selFeat.getStart() : (int) te.getChromosomePosition();
        final int end = editing ? selFeat.getEnd() : (int) Math.ceil(frame.getChromosomePosition(te.getMouseEvent().getX() + 1));

        String featDispName = editing && selFeat.getName() != null ? selFeat.getName() : "Feature";
        String name = editing ? "Edit " + featDispName : "Add annotation to " + MongoFeatureTrack.this.getName();
        JMenuItem item = new JMenuItem(name);
        item.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                Feature feat = editing ? selFeat.getDBFeature() : new BasicFeature(chr, start, end);
                FeatureAnnotDialog dialog = new FeatureAnnotDialog(IGV.getMainFrame(), getCollection(), feat);
                dialog.setVisible(true);
            }
        });
        return item;
    }

    private DBCollection getCollection() {
        return ((MongoFeatureSource) source).getCollection();
    }

    @Override
    public IGVPopupMenu getPopupMenu(TrackClickEvent te) {
        String title = getName();
        List<Track> tracks = Arrays.<Track>asList(this);
        IGVPopupMenu menu = TrackMenuUtils.getPopupMenu(tracks, title, te);

        //Add annotation or edit existing one
        Feature feat = getFeatureAtMousePosition(te);
        JMenuItem item = createEditAnnotMenuEntry(te, (DBFeature.IGVFeat) feat);
        if (item != null) menu.add(item);

        return menu;
    }
}
