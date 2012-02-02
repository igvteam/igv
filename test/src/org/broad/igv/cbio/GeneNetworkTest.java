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

import org.broad.tribble.readers.LineReader;
import org.jgrapht.Graph;
import org.jgrapht.graph.Pseudograph;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import static junit.framework.Assert.*;

/**
 * User: jacob
 * Date: 2012/02/02
 */
public class GeneNetworkTest {

    Graph<Vertex, BaseElement> graph;

    @Before
    public void setUp() throws Exception {
        graph = new Pseudograph<Vertex, BaseElement>(BaseElement.class);
    }

    @After
    public void tearDown() throws Exception {
        graph = null;
    }

    public void testGetCBioData() throws Exception {

    }


    public static void main(String[] args) throws IOException {

        //String[] gene_list = new String[]{"TP53", "BRCA1", "KRAS"};
        String filepath = "test/data/out/test_cbio.graphml";
        String[] gene_list = new String[]{"EGFR", "BRCA1"};

        LineReader reader = new GeneNetwork().getCBioData(gene_list);
        //Skip header
        String line = reader.readLine();

        Graph<Vertex, BaseElement> graph = new Pseudograph<Vertex, BaseElement>(BaseElement.class);
        KeyFactory edgeKeyFactory = new KeyFactory("edge");
        KeyFactory vertexKeyFactory = new KeyFactory("node");
        List<KeyFactory> factoryList = Arrays.asList(edgeKeyFactory, vertexKeyFactory);

        while ((line = reader.readLine()) != null) {
            String[] tokens = line.split("\\t");
            if (tokens.length != 3) {
                System.out.println("Bad line: " + line);
                continue;
            }
            String[] stringinfo = tokens[1].split(":");
            //TODO Validate line
            //Line labels, IN ORDER
            String[] names = new String[]{"type", "datasource", "evidence_codes", "pubmed_ids", "entrez_id_a", "entrez_id_b"};
            BaseElement edge = new Edge(edgeKeyFactory);
            for (int ii = 0; ii < names.length; ii++) {
                edge.put(names[ii], stringinfo[ii]);
            }

            Vertex v1 = new Vertex(tokens[0], vertexKeyFactory);
            Vertex v2 = new Vertex(tokens[2], vertexKeyFactory);

            graph.addVertex(v1);
            graph.addVertex(v2);
            graph.addEdge(v1, v2, edge);
        }

        GraphMLExporter exporter = new GraphMLExporter();
        exporter.exportGraph(filepath, graph, factoryList);
    }

    /**
     * Test that graph is consistent (no duplicate vertices or edges).
     * Two vertices should be considered equal if their names are
     * equal. We give them different attributes
     *
     * @throws Exception
     */
    @Test
    public void testGraphVertices() throws Exception {
        KeyFactory factory = new KeyFactory("testGraphVertices");
        Random random = new Random(132294492231L);

        int numv = 20;
        String[] vertex_labels = new String[numv];
        String[] seed_attribs = new String[]{"attr1", "score", "squizzle", "thenadierwalzoftreachery"};
        String[] test_attribs = new String[]{"attr1", "score", "squizzle", "thenadierwalzoftreachery"};

        for (int ii = 0; ii < numv; ii++) {
            vertex_labels[ii] = "v" + ii;
            Vertex v = new Vertex(vertex_labels[ii], factory);
            assertTrue(graph.addVertex(v));

            for (String attr : seed_attribs) {
                v.put(attr, random.nextDouble());
            }
        }

        assertEquals(numv, graph.vertexSet().size());

        for (String label : vertex_labels) {
            Vertex v = new Vertex(label, factory);
            for (String attr : test_attribs) {
                v.put(attr, random.nextDouble());
            }
            assertFalse(graph.addVertex(v));
        }

        assertEquals(numv, graph.vertexSet().size());


    }

}
