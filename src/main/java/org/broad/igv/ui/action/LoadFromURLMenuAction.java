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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.ui.action;

import org.apache.log4j.Logger;
import org.broad.igv.exceptions.HttpResponseException;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.google.GoogleUtils;
import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.IGVMenuBar;
import org.broad.igv.ui.util.LoadFromURLDialog;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.AmazonUtils;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.ResourceLocator;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import static org.broad.igv.util.AmazonUtils.isObjectAccessible;

/**
 * @author jrobinso
 */
public class LoadFromURLMenuAction extends MenuAction {

    static Logger log = Logger.getLogger(LoadFilesMenuAction.class);
    public static final String LOAD_FROM_DAS = "Load from DAS...";
    public static final String LOAD_FROM_URL = "Load from URL...";
    public static final String LOAD_FILE_AND_INDEX_FROM_URLS = "Load file and index from URLs...";
    public static final String LOAD_GENOME_FROM_URL = "Load Genome from URL...";
    private IGV igv;

    public LoadFromURLMenuAction(String label, int mnemonic, IGV igv) {
        super(label, null, mnemonic);
        this.igv = igv;
    }

    @Override
    public void actionPerformed(ActionEvent e) {

        JPanel ta = new JPanel();
        ta.setPreferredSize(new Dimension(600, 20));
        if (e.getActionCommand().equalsIgnoreCase(LOAD_FROM_URL)) {

            LoadFromURLDialog dlg = new LoadFromURLDialog(IGV.getMainFrame());
            dlg.setVisible(true);

            if (!dlg.isCanceled()) {

                String url = dlg.getFileURL();

                if (url != null && url.trim().length() > 0) {

                    url = mapURL(url.trim());

                    if (url.endsWith(".xml") || url.endsWith(".session")) {
                        try {
                            boolean merge = false;
                            String locus = null;
                            igv.doRestoreSession(url, locus, merge);
                        } catch (Exception ex) {
                            MessageUtils.showMessage("Error loading url: " + url + " (" + ex.toString() + ")");
                        }
                    } else {
                        try {
                            // If AWS support is active, check if objects are in accessible tiers via Load URL menu...
                            if (AmazonUtils.isAwsS3Path(url)) {
                                String bucket = AmazonUtils.getBucketFromS3URL(url);
                                String key = AmazonUtils.getKeyFromS3URL(url);

                                AmazonUtils.s3ObjectAccessResult res = isObjectAccessible(bucket, key);
                                if (!res.getObjAvailable()) { MessageUtils.showErrorMessage(res.getErrorReason(), null); return; }
                            }
                        } catch (NullPointerException npe) {
                            // User has not yet done Amazon->Login sequence
                            AmazonUtils.checkLogin();
                        }

                        ResourceLocator rl = new ResourceLocator(url.trim());

                        if (dlg.getIndexURL() != null) {
                            String indexUrl = dlg.getIndexURL().trim();

                            if (GoogleUtils.isGoogleCloud(indexUrl) || GoogleUtils.isGoogleDrive(indexUrl)) {
                                enableGoogleMenu();
                            }

                            rl.setIndexPath(indexUrl);
                        }
                        if (isHtsGet(url)) {
                            rl.setAttribute("htsget", true);
                        }
                        igv.loadTracks(Arrays.asList(rl));

                    }
                }
            }
        } else if ((e.getActionCommand().equalsIgnoreCase(LOAD_FROM_DAS))) {
            String url = JOptionPane.showInputDialog(IGV.getMainFrame(), ta, "Enter DAS feature source URL",
                    JOptionPane.QUESTION_MESSAGE);
            if (url != null && url.trim().length() > 0) {
                ResourceLocator rl = new ResourceLocator(url.trim());
                rl.setType("das");
                igv.loadTracks(Arrays.asList(rl));
            }
        } else if ((e.getActionCommand().equalsIgnoreCase(LOAD_GENOME_FROM_URL))) {

            String url = JOptionPane.showInputDialog(IGV.getMainFrame(), ta, "Enter URL to .genome or FASTA file",
                    JOptionPane.QUESTION_MESSAGE);
            if (url != null && url.trim().length() > 0) {
                try {
                    url = mapURL(url);
                    GenomeManager.getInstance().loadGenome(url.trim(), null);
                } catch (Exception e1) {
                    MessageUtils.showMessage("Error loading genome: " + e1.getMessage());
                }
            }
        }
    }
//  ToDo: need to figure out how to detect HtsGet resources
//    Also check : Supported htsget protocol version: vnd.ga4gh.htsget.v1.2.0may not be compatible with received content type: application/vnd.ga4gh.htsget.v0.2.0+json
    private boolean isHtsGet(String url) {
        return  url.contains("htsnexus.rnd.dnanex.us/v1/reads") ||
                url.contains("htsget.wtsi-npg-test.co.uk") ||
                url.contains("/reads/") ||
                url.contains("htsget.ga4gh.org");
    }

    private String mapURL(String url) {

        url = url.trim();

        if (GoogleUtils.isGoogleCloud(url) || GoogleUtils.isGoogleDrive(url)) {

            enableGoogleMenu();
            // if user is not currently logged in, attempt to
            // log in user
//            try {
//               OAuthUtils.getInstance().getProvider().doSecureLogin();
//            } catch (IOException e) {
//                log.error("Error connecting to OAuth: " + e.getMessage());
//            }

        }

        return url;
    }

    private void enableGoogleMenu() {

        if (!PreferencesManager.getPreferences().getAsBoolean(Constants.ENABLE_GOOGLE_MENU)) {
            PreferencesManager.getPreferences().put(Constants.ENABLE_GOOGLE_MENU, true);
            IGVMenuBar.getInstance().enableGoogleMenu(true);
        }
    }


    private boolean ping(String url) {
        InputStream is = null;
        try {
            Map<String, String> params = new HashMap();
            params.put("Range", "bytes=0-10");
            byte[] buffer = new byte[10];
            is = HttpUtils.getInstance().openConnectionStream(HttpUtils.createURL(url), params);
            is.read(buffer);
            is.close();
        } catch (HttpResponseException e1) {
            MessageUtils.showMessage(e1.getMessage());
            return false;
        } catch (IOException e) {
            log.error(e);

        } finally {
            if (is != null) try {
                is.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        return true;
    }
}

