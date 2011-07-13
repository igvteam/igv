/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package org.broad.igv.graph;

import org.broad.igv.feature.AbstractFeatureParser;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.Exon;
import org.broad.igv.feature.FeatureParser;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.util.ColorUtilities;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.tribble.Feature;
import org.broad.tribble.readers.AsciiLineReader;

import java.awt.*;
import java.io.IOException;
import java.util.*;
import java.util.List;

/**
 * @author jrobinso
 * @date Oct 19, 2010
 */
public class GeneUtils {

    static String file = "hg18_refGene.txt";
    static HashMap<String, List<BasicFeature>> transcripts;

    public static Graph getGraphFor(String gene, Genome genome) throws IOException {
        if (transcripts == null) {
            loadGenes(genome);
        }
        return createGraph(transcripts.get(gene));
    }


    // TODO -- account for strand

    private static Graph createGraph(List<BasicFeature> transcripts) {

        Graph graph = new Graph();

        // get unique exons.  Uniqueness is defined by locus
        HashSet<String> locusStrings = new HashSet();
        ArrayList<Exon> uniqExons = new ArrayList();
        for (BasicFeature t : transcripts) {
            for (Exon exon : t.getExons()) {
                String locus = exon.getLocusString();
                if (!locusStrings.contains(locus)) {
                    uniqExons.add(exon);
                    locusStrings.add(locus);
                }
            }
        }

        // Sort exons by start position
        Collections.sort(uniqExons, new Comparator<Exon>() {
            public int compare(Exon exon1, Exon exon2) {
                return exon1.getStart() - exon2.getStart();
            }
        });


        // Create nodes
        HashMap<String, Node> nodes = new HashMap();
        for (Exon exon : uniqExons) {
            String key = exon.getLocusString();
            nodes.put(key, new Node(exon.getStart()));
        }

        // Create edges and subgraphs
        int idx = 1;
        for (BasicFeature t : transcripts) {
            SubGraph sg = new SubGraph();
            Color color = ColorUtilities.randomColor(idx, 0.5f);
            List<Exon> exons = t.getExons();
            for (int i = 0; i < exons.size() - 1; i++) {
                Node n1 = nodes.get(exons.get(i).getLocusString());
                Node n2 = nodes.get(exons.get(i + 1).getLocusString());
                sg.addEdge(n1, n2, color);
                graph.addNode(n1);
                graph.addNode(n2);
            }
            idx++;

            graph.addSubGraph(sg);
        }


        graph.updateLayout();

        return graph;
    }

    private static void loadGenes(Genome genome) throws IOException {

        transcripts = new HashMap();

        FeatureParser fp = AbstractFeatureParser.getInstanceFor(file, genome);
        AsciiLineReader reader = ParsingUtils.openAsciiReader(new ResourceLocator(file));
        List<Feature> features = fp.loadFeatures(reader);

        for (Feature f : features) {


            // This cast is not safe, not for production use
            BasicFeature bf = (BasicFeature) f;
            List<BasicFeature> tlist = transcripts.get(bf.getName().toUpperCase());
            if (tlist == null) {
                tlist = new ArrayList();
                transcripts.put(bf.getName().toUpperCase(), tlist);
            }
            tlist.add(bf);
        }


    }

}
