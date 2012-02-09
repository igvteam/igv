/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTIES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package org.broad.igv.cbio;

import org.apache.commons.collections.Predicate;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.feature.FeatureDB;
import org.broad.igv.feature.NamedFeature;
import org.broad.igv.track.RegionScoreType;
import org.broad.igv.track.Track;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.Utilities;
import org.broad.tribble.Feature;
import org.w3c.dom.Document;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import javax.xml.parsers.ParserConfigurationException;
import java.io.IOException;
import java.io.InputStream;
import java.util.List;

/**
 * Class for taking a network of genes and annotating them
 * with additional information.
 * <p/>
 * User: jacob
 * Date: 2012/02/09
 */
public class NetworkAnnotator {
    private Logger logger = Logger.getLogger(NetworkAnnotator.class);

    private Document document;
    private Node graph;

    public static final String NODE_TAG = "node";
    public static final String KEY = "key";

    public static float collectScoreData(String name, List<Track> tracks, RegionScoreType type) {
        int zoom;

        List<NamedFeature> features = FeatureDB.getFeaturesList(name, Integer.MAX_VALUE);
        float totalScore = 0.0f;

        for (Feature feat : features) {
            double feat_score = 0;
            for (Track track : tracks) {
                float score = track.getRegionScore(feat.getChr(), feat.getStart(), feat.getEnd(), zoom = -1,
                        type, Globals.isHeadless() ? null : FrameManager.getDefaultFrame());
                feat_score += score > 0 ? score : 0;
            }
            totalScore += feat_score;
        }
        return totalScore / features.size();
    }

    /**
     * Download and store gene network from cbio
     *
     * @param path
     */
    public boolean loadNetwork(String path) {
        String error = null;
        try {
            InputStream cbioStream = ParsingUtils.openInputStream(new ResourceLocator(path));
            Document document = Utilities.createDOMDocumentFromXmlStream(cbioStream);
            this.document = document;
            this.graph = document.getElementsByTagName("graph").item(0);
        } catch (IOException e) {
            error = e.getMessage();
        } catch (ParserConfigurationException e) {
            error = e.getMessage();
        } catch (SAXException e) {
            error = e.getMessage();
        } finally {
            if (error != null) {
                logger.error(error);
                return false;
            } else {
                return true;
            }
        }
    }

    public NodeList getNodes() {
        return this.document.getElementsByTagName(NODE_TAG);
    }

    /**
     * Add the data specified by the score-types to our
     * network, using data from the tracks.
     * <p/>
     * TODO Our RegionScoreTypes do not map 1-1 with cbio data
     *
     * @param tracks
     * @param types
     */
    public void annotate(List<Track> tracks, List<RegionScoreType> types) {


    }

    private String getNodeKeyData(Node node, String key) {
        NodeList elements = node.getChildNodes();
        for (int ee = 0; ee < elements.getLength(); ee++) {
            Node el = elements.item(ee);
            try {
                NamedNodeMap map = el.getAttributes();
                Node label = map.getNamedItem(KEY);
                String textContent = label.getTextContent();
                if (textContent.equals(key)) {
                    return textContent;
                }
            } catch (NullPointerException e) {
                //In general these get hit due to newlines and such
                //We simply skip
                continue;
            }
        }
        return null;
    }

    public int filterNodes(Predicate predicate) {
        NodeList nodeList = this.getNodes();
        int count = 0;
        for (int nn = 0; nn < nodeList.getLength(); nn++) {
            Node node = nodeList.item(nn);
            if (!predicate.evaluate(node)) {
                this.graph.removeChild(node);
                count++;
            }
        }
        return count;
    }
}
