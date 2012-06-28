/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.feature;

import org.jgrapht.graph.DefaultDirectedGraph;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

/**
 * Takes a set of features, which are assumed to alternative
 * versions of the same feature. For example, multiple
 * alternative splicings of a single gene. Generates the simplest
 * graph for these features
 * <p/>
 * User: jacob
 * Date: 2012/03/28
 */
public class AlternativeSpliceGraph<T> extends DefaultDirectedGraph<IExon, Object> {

    private IExon lastExon = null;

    /**
     * For storing information about each exon.
     */
    private Map<IExon, T> parameters = new HashMap<IExon, T>();

    public AlternativeSpliceGraph() {
        super(Object.class);
    }


    public AlternativeSpliceGraph(Collection<? extends IGVFeature> features) {
        this();
        addFeatures(features);
    }

    /**
     * Add exon to graph, and parameter to map.
     * Note that even if exon already exists in
     * graph, the parameter will be overwritten
     *
     * @param exon
     * @param parameter
     * @return
     */
    public boolean put(IExon exon, T parameter) {
        parameters.put(exon, parameter);
        return addExon(exon);
    }

    public T getParameter(IExon exon) {
        return parameters.get(exon);
    }

    public boolean hasParameter(IExon exon) {
        return parameters.containsKey(exon);
    }

    /**
     * If an exon overlaps exactly a different exon,
     * we add appropriate edges.
     * <p/>
     * If we find a new exon, we add a node for it.
     * <p/>
     * If an exon partially overlaps, we
     * treat that as a different exon.
     * TODO Treat overlapping regions as same
     *
     * @param exon The exon to add
     */
    public boolean addExon(IExon exon) {

        //Should be a no-op if exon already there
        boolean added = addVertex(exon);

        if (lastExon != null) {
            if (getEdge(lastExon, exon) == null) {
                addEdge(lastExon, exon);
            }
        }
        lastExon = exon;

        return added;
    }

    private void addFeature(IGVFeature feature) {
        startFeature();
        for (Exon exon : feature.getExons()) {
            addExon(exon);
        }
    }


    /**
     * Convenience method for calling
     * {@link #addFeature} on each.
     *
     * @param features
     */
    private void addFeatures(Collection<? extends IGVFeature> features) {
        for (IGVFeature feature : features) {
            addFeature(feature);
        }
    }

    public void startFeature() {
        lastExon = null;
    }

}
