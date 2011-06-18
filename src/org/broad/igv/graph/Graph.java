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

import java.awt.*;
import java.util.*;
import java.util.List;

/**
 * @author jrobinso
 * @date Oct 12, 2010
 */
public class Graph {

    int maxX = Integer.MIN_VALUE;
    int minX = Integer.MAX_VALUE;

    Set<Node> nodes = new HashSet();
    List<SubGraph> subGraphs = new ArrayList();


    public void addNode(Node node) {
        nodes.add(node);
    }

    public Set<Node> getNodes() {
        return nodes;
    }


    /**
     * This method is responsible for updating the logical layout out the graph,  that is setting the "y" cell for each
     * node in the graph. This (the Y cell position) is the only free parameter.  It should be set to minimize
     * crossings of edges.
     */
    public void updateLayout() {
        alg2();


    }

    private void alg2() {

        NodeGrid grid = new NodeGrid();

        // The dumbest (nearly) possible algorithm, ignores edges altogether

        // Sort subgraphs by size
        Collections.sort(subGraphs, new Comparator<SubGraph>() {
            public int compare(SubGraph sg1, SubGraph sg2) {
                return sg2.getSize() - sg1.getSize();
            }
        });

        int level = 0;
        for (SubGraph sg : subGraphs) {
            for (Edge edge : sg.edges) {
                Node n1 = edge.getParent();
                if (n1.getCellY() == Integer.MIN_VALUE) {
                    grid.addNode(n1, level);
                }
                Node n2 = edge.getChild();
                if (n2.getCellY() == Integer.MIN_VALUE) {
                    grid.addNode(n2, level);
                }
            }
            level++;
        }

    }

    private void alg1() {
        // The dumbest (nearly) possible algorithm, ignores edges altogether

        // 1.  Segregate nodes by X position
        Map<Integer, List<Node>> nodeMap = new HashMap();
        for (Node node : nodes) {
            Integer cellX = node.getCellX();
            List<Node> nodeList = nodeMap.get(cellX);
            if (nodeList == null) {
                nodeList = new ArrayList();
                nodeMap.put(cellX, nodeList);
            }
            nodeList.add(node);
        }

        // 2.  Find the deepest stack of nodes (maxY)
        int maxY = 0;
        for (List<Node> nodeList : nodeMap.values()) {
            maxY = Math.max(nodeList.size(), maxY);
        }

        // 3. Simply assign the Y position based on the list order
        for (List<Node> nodeList : nodeMap.values()) {

            int initY = (maxY - nodeList.size()) / 2;

            for (int listPosition = 0; listPosition < nodeList.size(); listPosition++) {
                Node node = nodeList.get(listPosition);
                node.setCellY(initY + listPosition);
            }
        }
    }

    public void addSubGraph(SubGraph sg) {
        this.subGraphs.add(sg);
        maxX = Math.max(maxX, sg.getMaxX());
        minX = Math.min(minX, sg.getMinX());
    }

    public int getMaxX() {
        return maxX;
    }

    public int getMinX() {
        return minX;
    }

    public List<SubGraph> getSubGraphs() {
        return subGraphs;
    }


    class NodeGrid {

        Map<Integer, List<Node>> nodeMap = new HashMap();

        public void addNode(Node node, int level) {
            int x = node.getCellX();
            List<Node> nodeList = nodeMap.get(x);
            if(nodeList == null) {
                nodeList = new ArrayList();
                nodeMap.put(x, nodeList);
            }
            int y = nodeList.size(); //Math.max(level, nodeList.size());
            node.setCellY(y);
            nodeList.add(node);
        }


    }
}
