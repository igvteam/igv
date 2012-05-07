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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.session;

import org.broad.igv.renderer.*;

/**
 * Maps between renderer type name and renderer class.  Used for saving and restoring sessions.
 *
 * @author eflakes
 */
public class RendererFactory {


    static public enum RendererType {

        BAR_CHART,
        BASIC_FEATURE,
        FEATURE_DENSITY,
        GENE_TRACK,
        GISTIC_TRACK,
        HEATMAP,
        MUTATION,
        SCATTER_PLOT,
        LINE_PLOT
    }

    static Class defaultRendererClass = BarChartRenderer.class;

    static public Class getRendererClass(String rendererTypeName) {

        String typeName = rendererTypeName.toUpperCase();

        // Try know type names
        if (typeName.equals(RendererType.BAR_CHART.name()) || typeName.equals("BAR")) {
            return BarChartRenderer.class;
        } else if (typeName.equals(RendererType.BASIC_FEATURE.name())) {
            return IGVFeatureRenderer.class;
        } else if (typeName.equals(RendererType.FEATURE_DENSITY.name())) {
            return FeatureDensityRenderer.class;
        } else if (typeName.equals(RendererType.GENE_TRACK.name())) {
            return GeneTrackRenderer.class;
        } else if (typeName.equals(RendererType.GISTIC_TRACK.name())) {
            return GisticTrackRenderer.class;
        } else if (typeName.equals(RendererType.HEATMAP.name())) {
            return HeatmapRenderer.class;
        } else if (typeName.equals(RendererType.MUTATION.name())) {
            return MutationRenderer.class;
        } else if (typeName.equals(RendererType.SCATTER_PLOT.name()) ||
                typeName.toUpperCase().equals("POINTS")) {
            return PointsRenderer.class;
        } else if (typeName.equals(RendererType.LINE_PLOT.name()) ||
                typeName.toUpperCase().equals("LINE")) {
            return LineplotRenderer.class;
        }
        return null;

    }

    static public RendererType getRenderType(Renderer renderer) {

        Class rendererClass = renderer.getClass();

        RendererType rendererType = null;

        if (rendererClass.equals(BarChartRenderer.class)) {
            rendererType = RendererType.BAR_CHART;
        } else if (rendererClass.equals(IGVFeatureRenderer.class)) {
            rendererType = RendererType.BASIC_FEATURE;
        } else if (rendererClass.equals(FeatureDensityRenderer.class)) {
            rendererType = RendererType.FEATURE_DENSITY;
        } else if (rendererClass.equals(GeneTrackRenderer.class)) {
            rendererType = RendererType.GENE_TRACK;
        } else if (rendererClass.equals(GisticTrackRenderer.class)) {
            rendererType = RendererType.GISTIC_TRACK;
        } else if (rendererClass.equals(HeatmapRenderer.class)) {
            rendererType = RendererType.HEATMAP;
        } else if (rendererClass.equals(MutationRenderer.class)) {
            rendererType = RendererType.MUTATION;
        } else if (rendererClass.equals(PointsRenderer.class)) {
            rendererType = RendererType.SCATTER_PLOT;
        } else if (rendererClass.equals(LineplotRenderer.class)) {
            rendererType = RendererType.LINE_PLOT;
        }
        return rendererType;
    }
}
