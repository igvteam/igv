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

import org.broad.igv.util.HttpUtils;
import org.broad.tribble.readers.AsciiLineReader;
import org.broad.tribble.readers.LineReader;
import org.jgrapht.EdgeFactory;
import org.jgrapht.Graph;
import org.jgrapht.graph.Pseudograph;

import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.net.URLEncoder;

/**
 * Fetch a gene network from cBio portal. WORK IN PROGRESS.
 * User: jacob
 * Date: 2012/01/31
 */
public class GeneNetwork {

    private static final String BASE_URL = "http://www.cbioportal.org/public-portal/webservice.do";
    private static final String CMD = "getNetwork";
    private static final String GENE_LIST_KEY = "gene_list";
    private static final String BASIC_URL = BASE_URL + "?cmd=" + CMD + "&" + GENE_LIST_KEY + "=";

    public static void main(String[] args) throws IOException {

        //String[] gene_list = new String[]{"TP53", "BRCA1", "KRAS"};
        String filepath = "test/data/out/test_cbio.graphml";
        String[] gene_list = new String[]{"EGFR", "BRCA1"};
        String string_gl = URLEncoder.encode(gene_list[0], "UTF-8");
        for (int gi = 1; gi < gene_list.length; gi++) {
            string_gl += "," + URLEncoder.encode(gene_list[gi], "UTF-8");
        }

        Graph<String, BaseEdge> graph = new Pseudograph<String, BaseEdge>(BaseEdge.class);

        String url = BASIC_URL + string_gl;
        InputStream is = HttpUtils.getInstance().openConnectionStream(new URL(url));
        LineReader reader = new AsciiLineReader(is);
        //Skip header
        String line = reader.readLine();

        while ((line = reader.readLine()) != null) {
            String[] tokens = line.split("\\t");
            if (tokens.length != 3) {
                System.out.println("Bad line: " + line);
                continue;
            }
            String[] stringinfo = tokens[1].split(":");
            //TODO Validate line
            String[] names = new String[]{"type", "pubmed_id", "datasource", "experimental_type"};
            BaseEdge edge = new Edge();
            for (int ii = 0; ii < names.length; ii++) {
                edge.put(names[ii], new GraphMLData(stringinfo[ii]));
            }

            graph.addVertex(tokens[0]);
            graph.addVertex(tokens[2]);
            graph.addEdge(tokens[0], tokens[2], edge);
        }

        GraphMLExporter<String, BaseEdge> exporter = new GraphMLExporter<String, BaseEdge>(String.class, BaseEdge.class);
        exporter.exportGraph(filepath, graph);
    }

    public class testEdgeFactory implements EdgeFactory {

        public Object createEdge(Object o, Object o1) {
            return null; //TODO
        }
    }


}
