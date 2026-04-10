package org.igv.session;

import org.igv.feature.mut.MutationRenderer;
import org.igv.renderer.*;

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
        LINE_PLOT,
        DYNSEQ
    }

    static Class defaultRendererClass = BarChartRenderer.class;

    static public Class getRendererClass(String rendererTypeName) {

        String typeName = rendererTypeName.toUpperCase();

        // Try know type names
        if (typeName.equals(RendererType.BAR_CHART.name()) || typeName.equalsIgnoreCase("BAR")) {
            return BarChartRenderer.class;
        } else if (typeName.equalsIgnoreCase(RendererType.BASIC_FEATURE.name())) {
            return IGVFeatureRenderer.class;
        } else if (typeName.equalsIgnoreCase(RendererType.FEATURE_DENSITY.name())) {
            return FeatureDensityRenderer.class;
        } else if (typeName.equalsIgnoreCase(RendererType.GENE_TRACK.name())) {
            return GeneTrackRenderer.class;
        } else if (typeName.equals(RendererType.GISTIC_TRACK.name())) {
            return GisticTrackRenderer.class;
        } else if (typeName.equalsIgnoreCase(RendererType.HEATMAP.name())) {
            return HeatmapRenderer.class;
        } else if (typeName.equalsIgnoreCase(RendererType.MUTATION.name())) {
            return MutationRenderer.class;
        } else if (typeName.equals(RendererType.SCATTER_PLOT.name()) ||
                typeName.equalsIgnoreCase("POINTS")) {
            return PointsRenderer.class;
        } else if (typeName.equals(RendererType.LINE_PLOT.name()) ||
                typeName.equalsIgnoreCase("LINE")) {
            return LineplotRenderer.class;
        } else if (typeName.equalsIgnoreCase(RendererType.DYNSEQ.name()))  {
            return DynSeqRenderer.class;
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
        } else if (rendererClass.equals(DynSeqRenderer.class)) {
            rendererType = RendererType.DYNSEQ;
        }
        return rendererType;
    }
}
