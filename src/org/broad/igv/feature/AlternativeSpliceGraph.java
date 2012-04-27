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

import java.lang.reflect.InvocationHandler;
import java.lang.reflect.Method;
import java.lang.reflect.Proxy;
import java.util.Collection;

/**
 * Takes a set of features, which are assumed to alternative
 * versions of the same feature. For example, multiple
 * alternative splicings of a single gene. Generates the simplest
 * graph for these features
 * <p/>
 * User: jacob
 * Date: 2012/03/28
 */
public class AlternativeSpliceGraph extends DefaultDirectedGraph<IGVFeature, Object> {

    public AlternativeSpliceGraph(Collection<? extends IGVFeature> features) {
        super(Object.class);
        createFeatureGraph(features);
    }

    private void createFeatureGraph(Collection<? extends IGVFeature> features) {

        /**
         * We loop through features to build the graph.
         * If an exon overlaps exactly a different exon,
         * we add appropriate edges.
         *
         * If we find a new exon, we add a node for it.
         *
         * If an exon partially overlaps, we
         * treat that as a different exon.
         * TODO Treat overlapping regions as same
         *
         */
        for (IGVFeature feature : features) {

            IExon lastExon = null;
            for (Exon exon : feature.getExons()) {

                InvocationHandler handler = new ExonLocHandler(exon);
                IExon eProx = (IExon) Proxy.newProxyInstance(IExon.class.getClassLoader(),
                        new Class[]{IExon.class},
                        handler);

                //Should be a no-op if exon already there
                addVertex(eProx);

                if (lastExon != null) {
                    addEdge(lastExon, eProx);
                }
                lastExon = eProx;
            }
        }


    }

    private class ExonLocHandler implements InvocationHandler {

        private IExon parent;
        private int hashCode = 0;

        public ExonLocHandler(IExon parent) {
            this.parent = parent;
        }

        private boolean equals(IExon parent, Object inother) {
            if (inother == null || !(inother instanceof IExon)) {
                return false;
            }
            IExon other = (IExon) inother;
            boolean eq = parent.getChr().equals(other.getChr());
            eq &= parent.getStart() == other.getStart();
            eq &= parent.getEnd() == other.getEnd();
            eq &= parent.getCdStart() == other.getCdStart();
            eq &= parent.getCdEnd() == other.getCdEnd();
            eq &= parent.getStrand() == other.getStrand();
            return eq;
        }

        private int hashCode(IExon parent) {
            if (hashCode != 0) {
                return hashCode;
            }

            String conc = parent.getChr() + parent.getStrand().toString() + parent.getStart();
            conc += parent.getEnd();
            conc += parent.getCdStart();
            conc += parent.getCdEnd();
            int hc = conc.hashCode();

            if (hc == 0) {
                hc = 1;
            }
            hashCode = hc;
            return hc;
        }

        @Override
        public Object invoke(Object proxy, Method method, Object[] args) throws Throwable {
            if (method.getName().equals("hashCode")) {
                return hashCode(parent);
            } else if (method.getName().equals("equals")) {
                return equals(parent, args[0]);
            } else {
                return method.invoke(parent, args);
            }
        }
    }
}
