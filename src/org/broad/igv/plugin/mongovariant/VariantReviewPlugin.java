/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.plugin.mongovariant;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.dev.api.IGVPlugin;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackClickEvent;
import org.broad.igv.track.TrackMenuItemBuilder;
import org.broad.igv.track.TrackMenuUtils;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.PanelName;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.variant.Variant;
import org.broad.igv.variant.VariantTrack;
import org.broad.igv.variant.vcf.VCFVariant;
import org.broadinstitute.sting.gatk.walkers.na12878kb.core.NA12878DBArgumentCollection;
import org.broadinstitute.variant.variantcontext.VariantContext;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.Collection;

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

    @Override
    public void init(){
        if (showReviewOption) {
            //Test loading a class, won't work if running old Java
            if(!Globals.isVersionOrHigher("1.7")){
                log.error("VariantReviewPlugin requires Java 7 or higher. This plugin will be disabled");
                return;
            }
            //Test that we can load this class, plugin loader will catch the exception and
            //not load plugin if it can't be loaded
            NA12878DBArgumentCollection col = new NA12878DBArgumentCollection();
            initMenuItems();
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

    private void initMenuItems() {

        //Menu item for loading review track
        final JMenuItem loadReviewTrackItem = new JMenuItem("Load Review Track");
        loadReviewTrackItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                java.util.List<Track> newTracks = new ArrayList<Track>(1);
                ResourceLocator locator = new ResourceLocator(VariantReviewPlugin.getDbSpecPath());
                //TODO Put the track name in the dbSpec
                locator.setName("NA12878 KB");
                VariantReviewSource.loadVariantReview(locator, newTracks);
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
                //Not clear what to do with more than one track. Could be a bit nicer and only
                //require 1 variant track
                if(!isSingleVariantTrack(selectedTracks)){
                    return null;
                }
                Track track = selectedTracks.iterator().next();
                VariantTrack vTrack = (VariantTrack) track;
                final Variant variant = vTrack.getSelectedVariant(te);

                JMenuItem addReviewMenuItem = new JMenuItem("Submit Review to DB");
                addReviewMenuItem.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        VariantContext vc = VCFVariant.getVariantContext(variant);
                        (new VariantReviewDialog(IGV.getMainFrame(), vc)).setVisible(true);
                    }
                });
                addReviewMenuItem.setEnabled(variant != null);
                return addReviewMenuItem;
            }
        });
    }

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

}
