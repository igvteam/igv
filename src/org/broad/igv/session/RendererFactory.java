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
