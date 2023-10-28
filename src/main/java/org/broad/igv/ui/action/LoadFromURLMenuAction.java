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

import org.broad.igv.Globals;
import org.broad.igv.ext.ExtensionManager;
import org.broad.igv.ext.load.ILoadTracksFromUrlExtension;
import org.broad.igv.logging.*;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.util.GoogleUtils;
import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.session.SessionReader;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.IGVMenuBar;
import org.broad.igv.ui.util.LoadFromURLDialog;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.AmazonUtils;
import org.broad.igv.util.LongRunningTask;
import org.broad.igv.util.ResourceLocator;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.io.IOException;
import java.util.ArrayList;

import static org.broad.igv.util.AmazonUtils.isObjectAccessible;

/**
 * @author jrobinso
 */
public class LoadFromURLMenuAction extends MenuAction {

    static Logger log = LogManager.getLogger(LoadFilesMenuAction.class);
    public static final String LOAD_FROM_URL = "Load from URL...";
    public static final String LOAD_GENOME_FROM_URL = "Load Genome from URL...";
    public static final String LOAD_FROM_HTSGET = "Load from htsget Server...";
    private IGV igv;

    public LoadFromURLMenuAction(String label, int mnemonic, IGV igv) {
        super(label, null, mnemonic);
        this.igv = igv;
    }

    @Override
    public void actionPerformed(ActionEvent e) {

        JPanel ta = new JPanel();
        ta.setPreferredSize(new Dimension(600, 20));
        boolean isHtsGet = e.getActionCommand().equalsIgnoreCase(LOAD_FROM_HTSGET);
        if (e.getActionCommand().equalsIgnoreCase(LOAD_FROM_URL) || isHtsGet) {

            LoadFromURLDialog dlg = new LoadFromURLDialog(IGV.getInstance().getMainFrame(), isHtsGet);
            dlg.setVisible(true);

            if (!dlg.isCanceled()) {

                String inputURLs = dlg.getFileURL();
                if (inputURLs != null && inputURLs.trim().length() > 0) {

                    String[] inputs = Globals.whitespacePattern.split(inputURLs.trim());
                    checkURLs(inputs);
                    if (inputs.length == 1 && SessionReader.isSessionFile(inputs[0])) {
                        // Session URL
                        String url = inputs[0];
                        if (url.startsWith("s3://")) {
                            checkAWSAccessbility(url);
                        }
                        // extension code
                        ILoadTracksFromUrlExtension ext = (ILoadTracksFromUrlExtension) ExtensionManager.getExtentionFor(ILoadTracksFromUrlExtension.class, url);
                        if ( ext != null ) {
                            igv.loadTracks(ext.locatorsForUrl(url, dlg.getIndexURL()));
                        } else if (SessionReader.isSessionFile(url)) {

                            try {
                                LongRunningTask.submit(() -> this.igv.loadSession(url, null));
                            } catch (Exception ex) {
                                MessageUtils.showMessage("Error loading url: " + url + " (" + ex.toString() + ")");
                            }
                        }
                    } else {
                        // Files, possibly indexed
                        String[] indexes = null;
                        String indexURLs = dlg.getIndexURL();
                        if (indexURLs != null && indexURLs.trim().length() > 0) {
                            indexes = Globals.whitespacePattern.split(indexURLs.trim());
                            if (indexes.length != inputs.length) {
                                throw new RuntimeException("The number of Index URLs must equal the number of File URLs");
                            }
                            checkURLs(indexes);
                        }

                        ArrayList<ResourceLocator> locators = new ArrayList<>();
                        for (int i = 0; i < inputs.length; i++) {
                            String url = inputs[i];
                            ResourceLocator rl = new ResourceLocator(url.trim());
                            if (indexes != null) {
                                String indexUrl = indexes[i];
                                rl.setIndexPath(indexUrl);
                            }
                            if(isHtsGet) {
                                rl.setHtsget(true);
                            }
                            locators.add(rl);
                        }
                        igv.loadTracks(locators);
                    }
                }
            }
        } else if ((e.getActionCommand().equalsIgnoreCase(LOAD_GENOME_FROM_URL))) {

            String url = JOptionPane.showInputDialog(IGV.getInstance().getMainFrame(), ta, "Enter URL to .genome or FASTA file",
                    JOptionPane.QUESTION_MESSAGE);

            if (url != null && url.trim().length() > 0) {
                url = url.trim();
                try {
                    checkURLs(new String[]{url});
                    GenomeManager.getInstance().loadGenome(url, null);
                } catch (Exception e1) {
                    MessageUtils.showMessage("Error loading genome: " + e1.getMessage());
                }

            }
        }
    }

    private void checkURLs(String[] urls) {
        for (String url : urls) {
            if (url.startsWith("s3://")) {
                checkAWSAccessbility(url);
            } else if (url.startsWith("ftp://")) {
                MessageUtils.showMessage("FTP protocol is not supported");
            }
        }
    }
    public static String mapURL(String url) {

        url = url.trim();
        if (GoogleUtils.isGoogleDrive(url) || GoogleUtils.isGoogleDrive(url)) {
            enableGoogleMenu();
        }

        return url;
    }

    public static void enableGoogleMenu() {

        if (!PreferencesManager.getPreferences().getAsBoolean(Constants.ENABLE_GOOGLE_MENU)) {
            PreferencesManager.getPreferences().put(Constants.ENABLE_GOOGLE_MENU, true);
            try {
                IGVMenuBar.getInstance().enableGoogleMenu(true);
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
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
                }
            }
        } catch (NullPointerException npe) {
            // User has not yet done Amazon->Login sequence
            AmazonUtils.checkLogin();
        }
    }
}

