package org.igv.sashimi;

//~--- non-JDK imports --------------------------------------------------------

import htsjdk.tribble.Feature;
import org.igv.Globals;
import org.igv.feature.IExon;
import org.igv.feature.IGVFeature;
import org.igv.feature.SpliceJunctionFeature;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.prefs.PreferencesManager;
import org.igv.renderer.DataRange;
import org.igv.renderer.GraphicUtils;
import org.igv.renderer.IGVFeatureRenderer;
import org.igv.sam.AlignmentDataManager;
import org.igv.sam.AlignmentInterval;
import org.igv.sam.CoverageTrack;
import org.igv.track.RenderContext;
import org.igv.track.Track;
import org.igv.ui.FontManager;
import org.igv.ui.color.ColorUtilities;

import java.awt.*;
import java.awt.geom.Ellipse2D;
import java.awt.geom.GeneralPath;
import java.awt.geom.Rectangle2D;
import java.util.*;
import java.util.List;

import static org.igv.prefs.Constants.SASHIMI_SHOW_COVERAGE;


/**
 * Renderer for splice junctions. Draws a filled-in arc for each junction, with the width of the
 * arc representing the depth of coverage. If coverage information is present for the flanking
 * regions, draws that, too; otherwise indicates flanking regions with rectangles
 */
public class SashimiJunctionRenderer extends IGVFeatureRenderer {

    private static Logger log = LogManager.getLogger(SashimiJunctionRenderer.class);

    public enum ShapeType {
        CIRCLE,
        ELLIPSE,
        TEXT,
        NONE
    }


    // Static -- property is shared by all instances
    private static ShapeType shapeType = ShapeType.TEXT;

    Color ARC_COLOR_POS = new Color(150, 50, 50, 140); //transparent dull red

    private Color color;

    //central horizontal line color
    private Color colorCenterline;

    //maximum depth that can be displayed, due to track height limitations. Junctions with
    //this depth and deeper will all look the same
    protected int DEFAULT_MAX_DEPTH = 50;
    protected int maxDepth = DEFAULT_MAX_DEPTH;
    private Set<IExon> selectedExons;

    private CoverageTrack coverageTrack = null;
    private AlignmentDataManager dataManager = null;

    private Color background;

    /**
     * We want the features to alternate above and below, but don't want
     * the arcs to switch around when zooming /panning. So we store the above/below
     * status from the first rendering, and keep using that. This won't necessarily persist
     * between window openings, don't care.
     */
    private Map<Feature, Boolean> drawFeatureAbove = null;

    public SashimiJunctionRenderer() {
        this.darkMode = Globals.isDarkMode();
        this.color = darkMode ? ColorUtilities.modifyAlpha(Color.RED, 140) : ARC_COLOR_POS;
        this.colorCenterline = darkMode ?
                ColorUtilities.modifyAlpha(Color.WHITE, 100) :
                new Color(0, 0, 0, 100);
    }

    /**
     * If there are multiple arcs with the same start/end positions (e.g. different strands)
     * want to make sure they don't overlap
     */
    //private HashBasedTable<Integer, Integer, Feature> featureStartEndTable = null;
    public void setSelectedExons(Set<IExon> selectedExons) {
        this.selectedExons = selectedExons;
    }

    /**
     * Set the data manager
     *
     * @param dataManager
     */
    public void setDataManager(AlignmentDataManager dataManager) {
        this.dataManager = dataManager;
    }

    public AlignmentDataManager getDataManager() {
        return this.dataManager;
    }

    public void setBackground(Color background) {
        this.background = background;
    }

    public CoverageTrack getCoverageTrack() {
        return coverageTrack;
    }

    public int getMaxDepth() {
        return maxDepth;
    }

    public void setMaxDepth(int maxDepth) {
        this.maxDepth = maxDepth;
    }

    public static void setShapeType(ShapeType st) {
        shapeType = st;
    }

    public static ShapeType getShapeType() {
        return shapeType;
    }

    public Color getColor() {
        return color;
    }

    public void setColor(Color color) {
        this.color = color;
        if (this.coverageTrack != null) this.coverageTrack.setColor(color);
    }

    public void setCoverageTrack(CoverageTrack covTrack) {
        this.coverageTrack = new CoverageTrack(covTrack);
        this.coverageTrack.setColor(color);
        //Don't want to color SNPs, so we just set an impossibly high threshold
        coverageTrack.setSnpThreshold(2.0f);

        coverageTrack.setAutoScale(true);
        coverageTrack.setGlobalAutoScale(false);
    }


    /**
     * Note:  assumption is that featureList is sorted by pStart position.
     *
     * @param featureList
     * @param context
     * @param trackRectangle
     * @param track
     */
    @Override
    public void render(List<IGVFeature> featureList,
                       RenderContext context,
                       Rectangle trackRectangle,
                       Track track) {

        this.setColor(track.getColor());

        Rectangle coverageRectangle = new Rectangle(trackRectangle);

        final boolean showCoverage = PreferencesManager.getPreferences().getAsBoolean(SASHIMI_SHOW_COVERAGE);

        if (coverageTrack != null && showCoverage) {
            //Only want the coverage track to go so high so that the arcs still have room
            int newHeight = coverageRectangle.height / 2;
            int newY = coverageRectangle.y + coverageRectangle.height / 2 - newHeight;
            coverageRectangle.setBounds(coverageRectangle.x, newY, coverageRectangle.width, newHeight);

            coverageTrack.render(context, coverageRectangle);
        }

        double origin = context.getOrigin();
        double locScale = context.getScale();

        if ((featureList != null) && !featureList.isEmpty()) {

            // Create a graphics object to draw feature names.  Graphics are not cached
            // by font, only by color, so its necessary to create a new one to prevent
            // affecting other tracks.
            Font font = FontManager.getFont(track.getFontSize());
            Graphics2D fontGraphics = (Graphics2D) context.getGraphic2DForColor(Color.BLACK).create();
            fontGraphics.setFont(font);

            // Track coordinates
            double trackRectangleX = trackRectangle.getX();
            double trackRectangleMaxX = trackRectangle.getMaxX();

            //Draw track name
            fontGraphics.drawString(track.getName(), (int) (trackRectangleMaxX * 0.85), (int) (trackRectangle.getY() + font.getSize()));

            // Draw the lines that represent the bounds of
            // a feature's region
            Set<IExon> locselectedExons = selectedExons;

            if (drawFeatureAbove == null) drawFeatureAbove = new HashMap<Feature, Boolean>(featureList.size());
            //if(featureStartEndTable == null) featureStartEndTable = HashBasedTable.create();
            boolean lastDrawAbove = true;
            for (IGVFeature feature : featureList) {
                SpliceJunctionFeature junctionFeature = (SpliceJunctionFeature) feature;
                int junctionStart = junctionFeature.getJunctionStart();
                int junctionEnd = junctionFeature.getJunctionEnd();

                Boolean drawAbove = drawFeatureAbove.get(junctionFeature);
                if (drawAbove == null) {
                    drawAbove = !lastDrawAbove;
                    drawFeatureAbove.put(junctionFeature, drawAbove);
                }

                //Only show arcs for the selected feature, if applicable
                if (locselectedExons != null && locselectedExons.size() > 0) {
                    boolean inSelected = false;
                    for (IExon selectedExon : locselectedExons) {
                        if ((junctionStart >= selectedExon.getStart() && junctionStart <= selectedExon.getEnd())
                                || (junctionEnd >= selectedExon.getStart() && junctionEnd <= selectedExon.getEnd())) {
                            inSelected = true;
                            break;
                        }
                    }

                    if (!inSelected) continue;
                }

                double virtualPixelJunctionStart = Math.round((junctionStart - origin) / locScale);
                double virtualPixelJunctionEnd = Math.round((junctionEnd - origin) / locScale);

                // If the any part of the feature fits in the track rectangle draw it
                if ((virtualPixelJunctionEnd >= trackRectangleX) && (virtualPixelJunctionStart <= trackRectangleMaxX)) {

                    int depth = junctionFeature.getJunctionDepth();
                    Color color = feature.getColor();

                    //Calculate the height of the coverage track. This doesn't need to be exact,
                    //but we want the arcs to start at roughly the top of the coverage
                    int pixelYstart = 0;
                    int pixelYend = 0;

                    if (coverageTrack != null && showCoverage) {
                        pixelYstart = getYOffset(coverageRectangle, coverageTrack.getDataRange(), getCoverage(junctionStart));
                        pixelYend = getYOffset(coverageRectangle, coverageTrack.getDataRange(), getCoverage(junctionEnd));
                    }

                    drawFeature(pixelYstart, pixelYend,
                            (int) virtualPixelJunctionStart, (int) virtualPixelJunctionEnd, depth,
                            trackRectangle, context, color, drawAbove);
                    lastDrawAbove = drawAbove;
                }
            }

            //draw a central horizontal line
            Graphics2D g2D = context.getGraphic2DForColor(colorCenterline);
            g2D.drawLine((int) trackRectangleX, (int) trackRectangle.getCenterY(),
                    (int) trackRectangleMaxX, (int) trackRectangle.getCenterY());


        }
    }

    /**
     * Get the coverage around this approximate genome position. We actually look
     * for a maximum around a certain window. This is intended for plotting, we just
     * want the arc to look like it's coming from the top
     *
     * @param genomePos
     * @return
     */
    private int getCoverage(int genomePos) {

        if (dataManager == null) return 0;
        Collection<AlignmentInterval> intervals = dataManager.getLoadedIntervals();
        if (intervals == null) return 0;

        int buffer = 4;
        int coverage = 0;
        for (AlignmentInterval interval : intervals) {
            if (interval.overlaps(interval.getChr(), genomePos - buffer, genomePos + buffer)) {
                for (int loc = genomePos - buffer; loc < genomePos + buffer; loc++) {
                    coverage = Math.max(coverage, interval.getTotalCount(loc));
                }
                return coverage;
            }
        }

        return 0;
    }

    private int getYOffset(Rectangle rect, DataRange range, int totalCount) {

        double maxRange = range.isLog() ? Math.log10(range.getMaximum()) : range.getMaximum();
        double tmp = range.isLog() ? Math.log10(totalCount) / maxRange : totalCount / maxRange;
        int height = (int) (tmp * rect.height);

        height = Math.min(height, rect.height - 1);
        return -height;
    }

    /**
     * Draw a filled arc representing a single feature. The thickness and height of the arc are proportional to the
     * depth of coverage.  Some of this gets a bit arcane -- the result of lots of visual tweaking.
     *
     * @param pixelYStartOffset  the y offset from center line that the arc should start at
     * @param pixelYEndOffset    thy y offset from center line that the arc should end at
     * @param pixelJunctionStart the starting position of the junction, whether on-screen or not
     * @param pixelJunctionEnd   the ending position of the junction, whether on-screen or not
     * @param depth              coverage depth
     * @param trackRectangle
     * @param context
     * @param featureColor       the color specified for this feature.  May be null.
     * @param drawAbove          Whether to draw the arc above or below the centerline
     */
    protected void drawFeature(int pixelYStartOffset, int pixelYEndOffset,
                               int pixelJunctionStart, int pixelJunctionEnd, int depth,
                               Rectangle trackRectangle, RenderContext context, Color featureColor,
                               boolean drawAbove) {

        //If the feature color is specified, use it, except that we set our own alpha depending on whether
        //the feature is highlighted.  Otherwise default based on strand and highlight.
        Color color;
        if (featureColor != null) {
            int r = featureColor.getRed();
            int g = featureColor.getGreen();
            int b = featureColor.getBlue();
            int alpha = 140;
            color = new Color(r, g, b, alpha);
        } else {
            color = this.color;
        }

        int length = pixelJunctionEnd - pixelJunctionStart;
        int minArcHeight = (trackRectangle.height - 1) / 8;
        //We adjust the height slightly by length of junction, just so arcs don't overlap as much
        int arcHeight = minArcHeight + (length > 0 ? (int) (Globals.log2(length)) : 0);

        int minY = (int) trackRectangle.getCenterY() + Math.min(pixelYStartOffset - arcHeight, pixelYEndOffset - arcHeight);
        //Check if arc goes too high. All arcs going below have the same height,
        //so no need to check case-by-case
        if (drawAbove && minY < trackRectangle.getMinY()) {
            pixelYStartOffset = 0;
            pixelYEndOffset = 0;
        }

        //float depthProportionOfMax = Math.min(1, depth / maxDepth);
        int effDepth = Math.min(maxDepth, depth);

        //We adjust up or down depending on whether drawing up or down
        int yPosModifier = drawAbove ? -1 : 1;

        int arcBeginY = (int) trackRectangle.getCenterY() + yPosModifier + (drawAbove ? pixelYStartOffset - 2 : 0);
        int arcEndY = (int) trackRectangle.getCenterY() + yPosModifier + (drawAbove ? pixelYEndOffset - 2 : 0);

        Graphics2D g2D = context.getGraphic2DForColor(color);
        g2D.setBackground(background);

        double minStrokeSize = 0.1f;
        double maxStrokeSize = 3.0f;

        double strokeSize = maxStrokeSize;

        //Setting maxDepth to 1 should just max everything out, but 1/0 tends
        //to make things crash
        if (maxDepth >= 2) {
            double scale = (maxStrokeSize - minStrokeSize) / Math.log(maxDepth);
            strokeSize = scale * Math.log(effDepth) + minStrokeSize;
        }

        Stroke stroke = new BasicStroke((float) strokeSize);
        g2D.setStroke(stroke);

        if (pixelJunctionStart == pixelJunctionEnd) {
            // Junction is less than a pixel wide, draw a vertical line.
            int lineEndY = arcBeginY + yPosModifier * arcHeight;
            g2D.drawLine(pixelJunctionStart, arcBeginY, pixelJunctionStart, lineEndY);
        } else {
            //We use corners of a square as control points because why not
            //The control point is never actually reached
            int arcControlPeakY = arcBeginY + yPosModifier * arcHeight;

            GeneralPath arcPath = new GeneralPath();
            arcPath.moveTo(pixelJunctionStart, arcBeginY);
            arcPath.curveTo(pixelJunctionStart, arcControlPeakY,
                    pixelJunctionEnd, arcControlPeakY,
                    pixelJunctionEnd, arcEndY);
            g2D.draw(arcPath);
        }


        float midX = ((float) pixelJunctionStart + (float) pixelJunctionEnd) / 2;
        double actArcPeakY = arcBeginY + yPosModifier * Math.pow(0.5, 3) * (6) * arcHeight;
        int maxPossibleArcHeight = (trackRectangle.height - 1) / 4;
        float depthProportionOfMax = Math.min(1, (float) depth / maxDepth);
        float maxPossibleShapeHeight = maxPossibleArcHeight / 2.0f;
        Shape shape = null;

        switch (shapeType) {
            case CIRCLE:
                shape = createDepthCircle(maxPossibleShapeHeight, depthProportionOfMax, midX, actArcPeakY);
                break;
            case ELLIPSE:
                shape = createDepthEllipse(maxPossibleShapeHeight, depthProportionOfMax, midX, actArcPeakY);
                break;
            case TEXT:
                String text = "" + depth;
                Rectangle2D textBounds = g2D.getFontMetrics().getStringBounds(text, g2D);

                float floatX = (float) (midX - textBounds.getWidth() / 2);
                float floatY = (float) actArcPeakY + (float) textBounds.getHeight() / 2;

                //Clear background so we aren't drawing numbers over arcs
                int rectHeight = (int) textBounds.getHeight();

                int intX = (int) floatX;
                int intY = (int) floatY - rectHeight;
                int w = (int) textBounds.getWidth();
                int h = rectHeight;
                g2D.clearRect(intX, intY, w, h);
                GraphicUtils.drawCenteredText(text, intX, intY, w, h, g2D);

                break;
        }

        if (shape != null) {
            g2D.draw(shape);
            g2D.fill(shape);
        }
    }

    private Shape createDepthEllipse(double maxPossibleShapeHeight, double depthProportionOfMax, double arcMidX, double actArcPeakY) {
        double w = 5f;
        double x = arcMidX - w / 2;

        double h = maxPossibleShapeHeight * depthProportionOfMax;

        //The ellipse is always specified from the top left corner
        double y = actArcPeakY - h / 2;

        return new Ellipse2D.Double(x, y, w, h);
    }

    private Shape createDepthCircle(double maxPossibleShapeHeight, double depthProportionOfMax, double arcMidX, double actArcPeakY) {

        double h = maxPossibleShapeHeight * Math.sqrt(depthProportionOfMax);
        double w = h;
        double x = arcMidX - w / 2;

        //The ellipse is always specified from the top left corner
        double y = actArcPeakY - h / 2;

        return new Ellipse2D.Double(x, y, w, h);
    }


}
