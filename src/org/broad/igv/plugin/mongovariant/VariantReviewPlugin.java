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

package org.broad.igv.plugin.mongovariant;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.dev.api.IGVPlugin;
import org.broad.igv.dev.api.LoadHandler;
import org.broad.igv.track.*;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.PanelName;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.Utilities;
import org.broad.igv.variant.Variant;
import org.broad.igv.variant.VariantTrack;
import org.broad.igv.variant.vcf.VCFVariant;
import org.broadinstitute.gatk.tools.walkers.na12878kb.core.NA12878DBArgumentCollection;
import htsjdk.variant.variantcontext.VariantContext;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * User: jacob
 * Date: 2012-Dec-21
 */
public class VariantReviewPlugin implements IGVPlugin{

    private static Logger log = Logger.getLogger(VariantTrack.class);

    //TODO Experimental. Let user choose opinion and send info to DB
    private static final String SHOW_REVIEW_KEY = "SHOW_VARIANT_REVIEW";

    private static boolean hasReviewTrack = false;
    private static boolean showReviewOption = Boolean.parseBoolean(IGV.getInstance().getSession().getPersistent(SHOW_REVIEW_KEY, "false"));
    private static final String VARIANT_DB_EXT = "variant.db.json";
    private static ResourceLocator defaultLocator;

    @Override
    public void init(){
        if (showReviewOption) {
            //Test loading a class, won't work if running old Java
            if(!Globals.checkJavaVersion("1.7")){
                log.error("VariantReviewPlugin requires Java 7 or higher. This plugin will be disabled");
                return;
            }
            //Test that we can load this class, plugin loader will catch the exception and
            //not load plugin if it can't be loaded
            NA12878DBArgumentCollection col = new NA12878DBArgumentCollection();
            initMenuItems();

            TrackLoader.registerHandler(VARIANT_DB_EXT, new TrackLoadHandler());
        }
    }

    /**
     *
     * @param tracks
     * @return  true iff {@code tracks} is a single track of type {@code VariantTrack}
     */
    private boolean isSingleVariantTrack(Collection<Track> tracks){
        if (tracks.size() != 1) return false;

        Track track = tracks.iterator().next();

        return (track instanceof VariantTrack);
    }

    private static ResourceLocator loadVariantReviewTrack(String dbSpecPath,java.util.List<Track> newTracks, String trackName){
        ResourceLocator locator = new ResourceLocator(dbSpecPath);
        locator.setName(trackName);
        VariantReviewSource.loadVariantReview(locator, newTracks);
        return locator;
    }

    private void initMenuItems() {

        //Menu item for loading review track
        final JMenuItem loadReviewTrackItem = new JMenuItem("Load Review Track");

        loadReviewTrackItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                java.util.List<Track> newTracks = new ArrayList<Track>(1);
                ResourceLocator locator = loadVariantReviewTrack(VariantReviewPlugin.getDbSpecPath(), newTracks, "NA12878 KB");
                if(defaultLocator == null) defaultLocator = locator;
                IGV.getInstance().addTracks(newTracks, PanelName.DATA_PANEL);
                hasReviewTrack = true;
                loadReviewTrackItem.setEnabled(!hasReviewTrack);
                }
            });

        IGV.getInstance().addOtherToolMenu(loadReviewTrackItem);

        //Menu item for submitting new review
        TrackMenuUtils.addTrackMenuItemBuilder(new TrackMenuItemBuilder() {
            @Override
            public JMenuItem build(Collection<Track> selectedTracks, TrackClickEvent te) {
                if(!isSingleVariantTrack(selectedTracks)){
                    return null;
                }
                Track track = selectedTracks.iterator().next();
                final VariantTrack vTrack = (VariantTrack) track;
                final Variant variant = vTrack.getSelectedVariant(te);

                ResourceLocator trackLocator = vTrack.getResourceLocator();
                //If the track isn't connected to a database, use the default
                if(!trackLocator.getPath().toLowerCase().endsWith(VARIANT_DB_EXT) && defaultLocator != null){
                    trackLocator = defaultLocator;
                }
                String dbName = trackLocator.getName();
                final String dbSpecPath = trackLocator.getPath();
                JMenuItem addReviewMenuItem = new JMenuItem("Submit Review to " + dbName);
                addReviewMenuItem.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        VariantContext vc = VCFVariant.getVariantContext(variant);
                        (new VariantReviewDialog(IGV.getMainFrame(), vc, dbSpecPath)).setVisible(true);
                    }
                });
                addReviewMenuItem.setEnabled(variant != null);
                return addReviewMenuItem;
            }
        });
    }

    //All of this is is for loading from a db specified in the preferences / command line

    private static final String DB_PATH_KEY = "VARIANT_DB_PATH";
    public static final String DB_PATH_DEFAULT = NA12878DBArgumentCollection.DEFAULT_SPEC_PATH;

    private static final String PREFERENTIAL_SAMPLE_KEY = "PREFERENTIAL_SAMPLE";
    public static final String DEFAULT_PREFERENTIAL_SAMPLE = "NA12878";

    private static String preferentialSampleName;
    private static String dbSpecPath;

    public static String getPreferentialSampleName() {
        if(preferentialSampleName == null){
            preferentialSampleName = IGV.getPersistent(PREFERENTIAL_SAMPLE_KEY, DEFAULT_PREFERENTIAL_SAMPLE);
        }
        return preferentialSampleName;
    }

    public static String getDbSpecPath(){
        if(dbSpecPath == null){
            dbSpecPath = IGV.getPersistent(DB_PATH_KEY, DB_PATH_DEFAULT);
        }
        return dbSpecPath;
    }

    //------------//

    //For loading from a spec file
    private static class TrackLoadHandler implements LoadHandler{


        @Override
        public void load(String path, List<Track> newTracks) throws IOException {
            //TODO Put the track name in the dbSpec
            String name = Utilities.getFileNameFromURL(path);
            ResourceLocator locator = loadVariantReviewTrack(path, newTracks, name);
            if(defaultLocator == null) defaultLocator = locator;
        }
    }

}
