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
