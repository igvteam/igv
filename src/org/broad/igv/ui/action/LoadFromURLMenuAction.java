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
import org.broad.igv.PreferenceManager;
import org.broad.igv.exceptions.HttpResponseException;
import org.broad.igv.ga4gh.OAuthUtils;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.IGVMenuBar;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.ResourceLocator;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.io.IOException;
import java.io.InputStream;
import java.io.UnsupportedEncodingException;
import java.net.URL;
import java.net.URLEncoder;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

/**
 * @author jrobinso
 */
public class LoadFromURLMenuAction extends MenuAction {

    static Logger log = Logger.getLogger(LoadFilesMenuAction.class);
    public static final String LOAD_FROM_DAS = "Load from DAS...";
    public static final String LOAD_FROM_URL = "Load from URL...";
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
            String url = JOptionPane.showInputDialog(IGV.getMainFrame(), ta, "Enter URL (http or ftp)", JOptionPane.QUESTION_MESSAGE);

            if (url != null && url.trim().length() > 0) {
                url = url.trim();
                if (url.startsWith("gs://")) {
                    enableGoogleMenu();
                    url = translateGoogleCloudURL(url);
                }

                if (OAuthUtils.isGoogleCloud(url)) {

                    // Access a few bytes as a means to check authorization
                    if(!ping(url)) return;

                    if (url.indexOf("alt=media") < 0) {
                        url = url + (url.indexOf('?') > 0 ? "&" : "?") + "alt=media";
                    }

                }


                if (url.endsWith(".xml") || url.endsWith(".session")) {
                    try {
                        boolean merge = false;
                        String locus = null;
                        igv.doRestoreSession(url, locus, merge);
                    } catch (Exception ex) {
                        MessageUtils.showMessage("Error loading url: " + url + " (" + ex.toString() + ")");
                    }
                } else {
                    ResourceLocator rl = new ResourceLocator(url.trim());
                    igv.loadTracks(Arrays.asList(rl));

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
                    igv.loadGenome(url.trim(), null, true);
                } catch (IOException e1) {
                    MessageUtils.showMessage("Error loading genome: " + e1.getMessage());
                }
            }
        }
    }

    private void enableGoogleMenu() {

        if (!PreferenceManager.getInstance().getAsBoolean(PreferenceManager.ENABLE_GOOGLE_MENU)) {
            PreferenceManager.getInstance().put(PreferenceManager.ENABLE_GOOGLE_MENU, true);
            IGVMenuBar.getInstance().enableGoogleMenu(true);
        }
    }


    /**
     * gs://igv-bam-test/NA12878.bam
     * https://www.googleapis.com/storage/v1/b/igv-bam-test/o/NA12878.bam
     *
     * @param gsUrl
     * @return
     */
    private String translateGoogleCloudURL(String gsUrl) {

        int i = gsUrl.indexOf('/', 5);
        if (i < 0) {
            log.error("Invalid gs url: " + gsUrl);
            return gsUrl;
        }

        String bucket = gsUrl.substring(5, i);
        String object = gsUrl.substring(i + 1);
        try {
            object = URLEncoder.encode(object, "UTF8");
        } catch (UnsupportedEncodingException e) {
            // This isn't going to happen
            log.error(e);
        }

        return "https://www.googleapis.com/storage/v1/b/" + bucket + "/o/" + object;

    }


    private boolean ping(String url) {
        InputStream is = null;
        try {
            Map<String, String> params = new HashMap();
            params.put("Range", "0-10");
            byte [] buffer = new byte[10];
            is = HttpUtils.getInstance().openConnectionStream(new URL(url), params);
            is.read(buffer);
            is.close();
        }  catch (HttpResponseException e1) {
            MessageUtils.showMessage(e1.getMessage());
            return false;
        }
        catch (IOException e) {
            log.error(e);

        } finally {
            if(is != null) try {
                is.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        return true;
    }
}

