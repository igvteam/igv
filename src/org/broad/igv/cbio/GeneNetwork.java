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
import org.broad.igv.util.HttpUtils;
import org.broad.tribble.Feature;
import org.broad.tribble.readers.AsciiLineReader;
import org.broad.tribble.readers.LineReader;
import org.jgrapht.EdgeFactory;
import org.jgrapht.VertexFactory;
import org.jgrapht.graph.Pseudograph;

import java.io.IOException;
import java.io.InputStream;
import java.io.UnsupportedEncodingException;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLEncoder;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Fetch a gene network from cBio portal. WORK IN PROGRESS.
 * User: jacob
 * Date: 2012/01/31
 */
public class GeneNetwork extends Pseudograph<Vertex, BaseElement> {

    private Logger log = Logger.getLogger(GeneNetwork.class);

    static final String BASE_URL = "http://www.cbioportal.org/public-portal/webservice.do";
    static final String CMD = "getNetwork";
    static final String GENE_LIST_KEY = "gene_list";
    static final String BASIC_URL = BASE_URL + "?cmd=" + CMD + "&" + GENE_LIST_KEY + "=";
    /**
     * The middle column returned by CBIO is delimited by ":". These are the names, in order, of that data.
     */
    static final String[] COLUMN_NAMES = new String[]{"type", "datasource", "evidence_codes", "pubmed_ids", "entrez_id_a", "entrez_id_b"};
    static final int A_KEY = 0;
    static final int B_KEY = 2;
    static final int DATA_COL = 1;

    public static final String PERCENT_MUTATED = "PERCENT_MUTATED";

    protected KeyFactory edgeKeyFactory;
    protected KeyFactory vertexKeyFactory;
    protected List<KeyFactory> factoryList;

    public GeneNetwork() {
        this(BaseElement.class);
    }

    public GeneNetwork(EdgeFactory edgeFactory) {
        super(edgeFactory);
    }

    public GeneNetwork(Class clazz) {
        super(clazz);

        edgeKeyFactory = new KeyFactory("edge");
        vertexKeyFactory = new KeyFactory("node");
        factoryList = Arrays.asList(edgeKeyFactory, vertexKeyFactory);
    }

    public List<KeyFactory> getFactoryList() {
        return this.factoryList;
    }

    public boolean filterNodes(Predicate predicate) {
        Set<Vertex> rejected = new HashSet<Vertex>(this.vertexSet().size());
        for (Vertex v : this.vertexSet()) {
            if (!predicate.evaluate(v)) {
                rejected.add(v);
            }
        }
        return this.removeAllVertices(rejected);
    }

    public boolean filterEdges(Predicate predicate) {
        Set<BaseElement> rejected = new HashSet<BaseElement>(this.edgeSet().size());
        for (BaseElement e : this.edgeSet()) {
            if (!predicate.evaluate(e)) {
                rejected.add(e);
            }
        }
        return this.removeAllEdges(rejected);
    }

    public boolean pruneGraph() {
        Predicate unconnected = new Predicate<Vertex>() {
            public boolean evaluate(Vertex object) {
                return edgesOf(object).size() > 0;
            }
        };
        return this.filterNodes(unconnected);
    }

    public void collectMutationData(List<Track> tracks) {
        this.collectScoreData(tracks, RegionScoreType.MUTATION_COUNT, PERCENT_MUTATED);
    }

    public void collectScoreData(List<Track> tracks, RegionScoreType type, String data_label) {
        int zoom;

        Set<Vertex> rejected = new HashSet<Vertex>();
        for (Vertex v : this.vertexSet()) {
            List<NamedFeature> features = FeatureDB.getFeaturesList(v.getName(), Integer.MAX_VALUE);
            double total_tracks = 0;
            double total_score = 0;
            for (Feature feat : features) {
                for (Track track : tracks) {
                    float score = track.getRegionScore(feat.getChr(), feat.getStart(), feat.getEnd(), zoom = -1,
                            type, Globals.isHeadless() ? null : FrameManager.getDefaultFrame().getName());
                    total_tracks++;
                    total_score += score > 0 ? score : 0;
                }
            }
            double perc_score = 100.0 * total_score / total_tracks;
            v.put(data_label, perc_score);
        }
    }

    public int loadCBioLines(String[] gene_list) throws IOException {
        LineReader reader = getCBioLines(gene_list);
        int newVertices = 0;
        //Skip header
        String line = reader.readLine();


        while ((line = reader.readLine()) != null) {
            String[] tokens = line.split("\\t");
            if (tokens.length != 3) {
                log.error("Bad cbio line: " + line);
                continue;
            }
            String[] stringinfo = tokens[DATA_COL].split(":");
            BaseElement edge = new Edge(edgeKeyFactory);
            for (int ii = 0; ii < COLUMN_NAMES.length; ii++) {
                edge.put(COLUMN_NAMES[ii], stringinfo[ii]);
            }

            Vertex v1 = new Vertex(tokens[A_KEY], vertexKeyFactory);
            Vertex v2 = new Vertex(tokens[B_KEY], vertexKeyFactory);

            newVertices += this.addVertex(v1) ? 1 : 0;
            newVertices += this.addVertex(v2) ? 1 : 0;
            this.addEdge(v1, v2, edge);

        }

        return newVertices;
    }

    LineReader getCBioLines(String[] gene_list) throws IOException {

        try {
            String string_gl = URLEncoder.encode(gene_list[0], "UTF-8");
            for (int gi = 1; gi < gene_list.length; gi++) {
                string_gl += "," + URLEncoder.encode(gene_list[gi], "UTF-8");
            }

            String url = GeneNetwork.BASIC_URL + string_gl;
            //System.out.println(url);
            InputStream is = HttpUtils.getInstance().openConnectionStream(new URL(url));
            LineReader reader = new AsciiLineReader(is);
            return reader;
        } catch (UnsupportedEncodingException e) {
            throw new IllegalArgumentException("Bad argument in genelist: " + e.getMessage());
        } catch (MalformedURLException e) {
            //It's not a malformed URL. There's essentially no way it could be,
            //unless the encoding malfunctions but throws no exception
            return null;
        }
    }

    public static class SimpleVertexFactory implements VertexFactory<Vertex> {
        private KeyFactory keyFactory;
        private long vCounter = 0;

        public SimpleVertexFactory(KeyFactory keyFactory) {
            this.keyFactory = keyFactory;
        }

        public Vertex createVertex() {
            Vertex v = new Vertex("v" + vCounter, keyFactory);
            vCounter++;
            return v;
        }
    }

    public static class SimpleEdgeFactory implements EdgeFactory<Vertex, BaseElement> {
        private KeyFactory keyFactory;
        private long eCounter = 0;

        public SimpleEdgeFactory(KeyFactory keyFactory) {
            this.keyFactory = keyFactory;
        }

        public BaseElement createEdge(Vertex sourceVertex, Vertex targetVertex) {
            BaseElement e = new Edge(this.keyFactory);
            e.put("label", "" + eCounter);
            eCounter++;
            return e;
        }
    }

}
