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

import org.broad.igv.logging.*;
import org.broad.igv.exceptions.HttpResponseException;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.google.GoogleUtils;
import org.broad.igv.google.OAuthProvider;
import org.broad.igv.google.OAuthUtils;
import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.session.SessionReader;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.IGVMenuBar;
import org.broad.igv.ui.util.LoadFromURLDialog;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.AmazonUtils;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.LongRunningTask;
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

    static Logger log = LogManager.getLogger(LoadFilesMenuAction.class);
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

            LoadFromURLDialog dlg = new LoadFromURLDialog(IGV.getInstance().getMainFrame());
            dlg.setVisible(true);

            if (!dlg.isCanceled()) {

                final String inputURL = dlg.getFileURL();

                if (inputURL != null && inputURL.trim().length() > 0) {

                    final String url = mapURL(inputURL.trim());

                    if (url.startsWith("s3://")) {
                        checkAWSAccessbility(url);
                    }

                    if (SessionReader.isSessionFile(url)) {
                        try {
                            LongRunningTask.submit(() -> this.igv.loadSession(url, null));
                        } catch (Exception ex) {
                            MessageUtils.showMessage("Error loading url: " + url + " (" + ex.toString() + ")");
                        }
                    } else {
                        ResourceLocator rl = new ResourceLocator(url.trim());

                        if (dlg.getIndexURL() != null) {
                            String indexUrl = dlg.getIndexURL().trim();

                            if (GoogleUtils.isGoogleCloud(indexUrl) || GoogleUtils.isGoogleDrive(indexUrl)) {
                                enableGoogleMenu();
                            }

                            rl.setIndexPath(indexUrl);
                        }
                        igv.loadTracks(Arrays.asList(rl));

                    }
                }
            }
        } else if ((e.getActionCommand().equalsIgnoreCase(LOAD_FROM_DAS))) {
            String url = JOptionPane.showInputDialog(IGV.getInstance().getMainFrame(), ta, "Enter DAS feature source URL",
                    JOptionPane.QUESTION_MESSAGE);
            if (url != null && url.trim().length() > 0) {
                ResourceLocator rl = new ResourceLocator(url.trim());
                rl.setFormat("das");
                igv.loadTracks(Arrays.asList(rl));
            }
        } else if ((e.getActionCommand().equalsIgnoreCase(LOAD_GENOME_FROM_URL))) {
            String url = JOptionPane.showInputDialog(IGV.getInstance().getMainFrame(), ta, "Enter URL to .genome or FASTA file",
                    JOptionPane.QUESTION_MESSAGE);
            if (url != null && url.trim().length() > 0) {
                if (url.startsWith("s3://")) {
                    checkAWSAccessbility(url);
                } else if (url.startsWith("ftp://")) {
                    MessageUtils.showMessage("FTP protocol is not supported");
                }
                try {
                    url = mapURL(url);
                    GenomeManager.getInstance().loadGenome(url.trim(), null);
                } catch (Exception e1) {
                    MessageUtils.showMessage("Error loading genome: " + e1.getMessage());
                }

            }
        }
    }

    private String mapURL(String url) {

        url = url.trim();

        OAuthProvider oauthProvider = OAuthUtils.getInstance().getProvider();
        if (GoogleUtils.isGoogleCloud(url) || GoogleUtils.isGoogleDrive(url)) {

            enableGoogleMenu();
        }
        else if(oauthProvider != null && oauthProvider.appliesToUrl(url)){
            // if user is not currently logged in, attempt to
            // log in user
           try {
              oauthProvider.doSecureLogin();
           } catch (Exception e) {
                log.error("Error connecting to OAuth: " + e.getMessage());
           }

        }

        return url;
    }

    private void enableGoogleMenu() {

        if (!PreferencesManager.getPreferences().getAsBoolean(Constants.ENABLE_GOOGLE_MENU)) {
            PreferencesManager.getPreferences().put(Constants.ENABLE_GOOGLE_MENU, true);
            IGVMenuBar.getInstance().enableGoogleMenu(true);
        }
    }


    private void checkAWSAccessbility(String url) {
        try {
            // If AWS support is active, check if objects are in accessible tiers via Load URL menu...
            if (AmazonUtils.isAwsS3Path(url)) {
                String bucket = AmazonUtils.getBucketFromS3URL(url);
                String key = AmazonUtils.getKeyFromS3URL(url);
                AmazonUtils.s3ObjectAccessResult res = isObjectAccessible(bucket, key);
                if (!res.isObjectAvailable()) {
                    MessageUtils.showErrorMessage(res.getErrorReason(), null);
                    return;
                }
            }
        } catch (NullPointerException npe) {
            // User has not yet done Amazon->Login sequence
            AmazonUtils.checkLogin();
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

