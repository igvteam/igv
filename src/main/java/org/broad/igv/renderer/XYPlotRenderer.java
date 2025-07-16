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
package org.broad.igv.renderer;

//~--- non-JDK imports --------------------------------------------------------

import org.broad.igv.Globals;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.prefs.IGVPreferences;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.Track;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.UIConstants;

import java.awt.*;
import java.text.DecimalFormat;
import java.util.List;

import static org.broad.igv.prefs.Constants.*;

/**
 * @author jrobinso
 */
public abstract class XYPlotRenderer extends DataRenderer {

    private double marginFraction = 0.2;


    protected void drawDataPoint(Color graphColor, int dx, int pX, int baseY, int pY,
                                 RenderContext context) {
        context.getGraphic2DForColor(graphColor).fillRect(pX, pY, dx, 2);

    }

    /**
     * Render the track in the given rectangle.
     *
     * @param track
     * @param locusScores
     * @param context
     * @param arect
     */
    public synchronized void renderScores(Track track, List<LocusScore> locusScores, RenderContext context, Rectangle arect) {


        Rectangle adjustedRect = calculateDrawingRect(arect);
        double origin = context.getOrigin();
        double locScale = context.getScale();

        Color posColor = track.getColor();
        Color negColor = track.getAltColor();

        // Get the Y axis definition, consisting of minimum, maximum, and base value.  Often
        // the base value is == min value which is == 0.

        DataRange dataRange = track.getDataRange();
        float maxValue = dataRange.getMaximum();
        float baseValue = dataRange.getBaseline();
        float minValue = dataRange.getMinimum();
        boolean isLog = dataRange.isLog();

        if (isLog) {
            minValue = (float) (minValue == 0 ? 0 : Math.log10(minValue));
            maxValue = (float) Math.log10(maxValue);
        }


        // Calculate the Y scale factor.

        double delta = (maxValue - minValue);
        double yScaleFactor = adjustedRect.getHeight() / delta;

        // Calculate the Y position in pixels of the base value.  Clip to bounds of rectangle
        double baseDelta = maxValue - baseValue;
        int baseY = (int) (adjustedRect.getY() + baseDelta * yScaleFactor);
        if (baseY < adjustedRect.y) {
            baseY = adjustedRect.y;
        } else if (baseY > adjustedRect.y + adjustedRect.height) {
            baseY = adjustedRect.y + adjustedRect.height;
        }

        for (LocusScore score : locusScores) {

            // Note -- don't cast these to an int until the range is checked.
            // could get an overflow.
            double pX = ((score.getStart() - origin) / locScale);
            double dx = Math.ceil((Math.max(1, score.getEnd() - score.getStart())) / locScale) + 1;
            if ((pX + dx < 0)) {
                continue;
            } else if (pX > adjustedRect.getMaxX()) {
                break;
            }

            float dataY = score.getScore();
            if (isLog && dataY <= 0) {
                continue;
            }

            if (!Float.isNaN(dataY)) {

                // Compute the pixel y location.  Clip to bounds of rectangle.
                double dy = isLog ? Math.log10(dataY) - baseValue : (dataY - baseValue);
                int pY = baseY - (int) (dy * yScaleFactor);
                if (pY < adjustedRect.y) {
                    pY = adjustedRect.y;
                } else if (pY > adjustedRect.y + adjustedRect.height) {
                    pY = adjustedRect.y + adjustedRect.height;
                }

                Color color = (dataY >= baseValue) ? posColor : negColor;
                drawDataPoint(color, (int) dx, (int) pX, baseY, pY, context);

            }

        }


    }

    static DecimalFormat formatter = new DecimalFormat();

    /**
     * Method description
     *
     * @param track
     * @param context
     * @param arect
     */
    @Override
    public void renderAxis(Track track, RenderContext context, Rectangle arect) {

        if (context.multiframe) {
            return;
        }

        super.renderAxis(track, context, arect);

        Rectangle drawingRect = calculateDrawingRect(arect);

        IGVPreferences prefs = PreferencesManager.getPreferences();

        Color labelColor = prefs.getAsBoolean(CHART_COLOR_TRACK_NAME) ? track.getColor() : Color.black;
        Graphics2D labelGraphics = context.getGraphic2DForColor(labelColor);

        labelGraphics.setFont(FontManager.getFont(8));

        if (prefs.getAsBoolean(CHART_DRAW_TRACK_NAME)) {

            // Only attempt if track height is > 25 pixels
            if (arect.getHeight() > 25) {
                Rectangle labelRect = new Rectangle(arect.x, arect.y + 10, arect.width, 10);
                labelGraphics.setFont(FontManager.getFont(10));
                GraphicUtils.drawCenteredText(track.getName(), labelRect, labelGraphics);
            }
        }

        if (prefs.getAsBoolean(CHART_DRAW_Y_AXIS)) {

            Rectangle axisRect = new Rectangle(arect.x, arect.y + 1, AXIS_AREA_WIDTH, arect.height);

            DataRange axisDefinition = track.getDataRange();
            float maxValue = axisDefinition.getMaximum();
            float baseValue = axisDefinition.getBaseline();
            float minValue = axisDefinition.getMinimum();
            
            // Bottom (minimum tick mark)
            int pY = computeYPixelValue(drawingRect, axisDefinition, minValue);

            labelGraphics.drawLine(axisRect.x + AXIS_AREA_WIDTH - 10, pY, axisRect.x + AXIS_AREA_WIDTH - 5, pY);
            GraphicUtils.drawRightJustifiedText(formatter.format(minValue), axisRect.x + AXIS_AREA_WIDTH - 15, pY, labelGraphics);

            // Top (maximum tick mark)
            int topPY = computeYPixelValue(drawingRect, axisDefinition, maxValue);

            labelGraphics.drawLine(axisRect.x + AXIS_AREA_WIDTH - 10, topPY,
                    axisRect.x + AXIS_AREA_WIDTH - 5, topPY);
            GraphicUtils.drawRightJustifiedText(formatter.format(maxValue),
                    axisRect.x + AXIS_AREA_WIDTH - 15, topPY + 4, labelGraphics);

            // Connect top and bottom
            labelGraphics.drawLine(axisRect.x + AXIS_AREA_WIDTH - 10, topPY,
                    axisRect.x + AXIS_AREA_WIDTH - 10, pY);

            // Middle tick mark.  Draw only if room
            int midPY = computeYPixelValue(drawingRect, axisDefinition, baseValue);

            if ((midPY < pY - 15) && (midPY > topPY + 15)) {
                labelGraphics.drawLine(axisRect.x + AXIS_AREA_WIDTH - 10, midPY,
                        axisRect.x + AXIS_AREA_WIDTH - 5, midPY);
                GraphicUtils.drawRightJustifiedText(formatter.format(baseValue),
                        axisRect.x + AXIS_AREA_WIDTH - 15, midPY + 4, labelGraphics);
            }

        } else if (track.isShowDataRange() && arect.height > 20) {
            drawScale(track.getDataRange(), context, arect);
        }
    }

    @Override
    public void renderBorder(Track track, RenderContext context, Rectangle arect) {

        Rectangle adjustedRect = calculateDrawingRect(arect);

        // Draw boundaries if there is room
        if (adjustedRect.getHeight() >= 10) {

            ///TrackProperties pros = track.getProperties();


            // midline

            DataRange axisDefinition = track.getDataRange();
            float maxValue = axisDefinition.getMaximum();
            float baseValue = axisDefinition.getBaseline();
            float minValue = axisDefinition.getMinimum();


            double maxX = adjustedRect.getMaxX();
            double x = adjustedRect.getX();
            double y = adjustedRect.getY();

            if ((baseValue > minValue) && (baseValue < maxValue)) {
                int baseY = computeYPixelValue(adjustedRect, axisDefinition, baseValue);

                getBaselineGraphics(context).drawLine((int) x, baseY, (int) maxX, baseY);
            }

            IGVPreferences prefs = PreferencesManager.getPreferences();

            Color altColor = track.getAltColor();
            Color borderColor = (prefs.getAsBoolean(CHART_COLOR_BORDERS) && altColor != null && altColor.equals(track.getColor()))
                    ? track.getColor() : Color.lightGray;
            Graphics2D borderGraphics = context.getGraphic2DForColor(borderColor);

            // Draw the baseline -- todo, this is a wig track option?
            double zeroValue = axisDefinition.getBaseline();
            int zeroY = computeYPixelValue(adjustedRect, axisDefinition, zeroValue);
            borderGraphics.drawLine(adjustedRect.x, zeroY, adjustedRect.x + adjustedRect.width, zeroY);

            // Optionally draw "Y" line  (UCSC track line option)
            if (track.isDrawYLine()) {
                Graphics2D yLineGraphics = context.getGraphic2DForColor(Color.gray);
                int yLine = computeYPixelValue(adjustedRect, axisDefinition, track.getYLine());
                GraphicUtils.drawDashedLine(borderGraphics, adjustedRect.x, yLine, adjustedRect.x + adjustedRect.width, yLine);
            }


            // If the chart has + and - numbers draw both borders or none. This
            // needs documented somewhere.
            boolean drawBorders = true;

            if (minValue * maxValue < 0) {
                drawBorders = prefs.getAsBoolean(CHART_DRAW_BOTTOM_BORDER) &&
                        prefs.getAsBoolean(CHART_DRAW_TOP_BORDER);
            }

            if (drawBorders && prefs.getAsBoolean(CHART_DRAW_TOP_BORDER)) {
                borderGraphics.drawLine(adjustedRect.x, adjustedRect.y,
                        adjustedRect.x + adjustedRect.width, adjustedRect.y);
            }

            if (drawBorders && prefs.getAsBoolean(CHART_DRAW_BOTTOM_BORDER)) {
                borderGraphics.drawLine(adjustedRect.x, adjustedRect.y + adjustedRect.height,
                        adjustedRect.x + adjustedRect.width,
                        adjustedRect.y + adjustedRect.height);
            }
        }
        /*
        (CHART_DRAW_TOP_BORDER));
        prefs.setDrawBottomBorder(getBooleanPreference(CHART_DRAW_BOTTOM_BORDER));
        prefs.setColorBorders(getBooleanPreference(CHART_COLOR_BORDERS));
        prefs.setDrawAxis(getBooleanPreference(CHART_DRAW_Y_AXIS));
        prefs.setDrawTrackName(getBooleanPreference(CHART_DRAW_TRACK_NAME));
        prefs.setColorTrackName(getBooleanPreference(CHART_COLOR_TRACK_NAME));
        prefs.setAutoscale(getBooleanPreference(CHART_AUTOSCALE));
        prefs.setShowDataRange(getBooleanPreference(CHART_SHOW_DATA_RANGE));
         */
    }

    /**
     * Get a graphics object for the baseline.
     * TODO -- make the line style settable by the user
     *
     * @param context
     * @return
     */
    private static Graphics2D getBaselineGraphics(RenderContext context) {
        Graphics2D baselineGraphics;
        baselineGraphics = (Graphics2D) context.getGraphic2DForColor(Color.lightGray).create();
        return baselineGraphics;
    }

    /**
     * Method description
     *
     * @return
     */
    public String getDisplayName() {
        return "Scatter Plot";
    }

    protected int computeYPixelValue(Rectangle drawingRect, DataRange axisDefinition, double dataY) {

        double maxValue = axisDefinition.getMaximum();
        double minValue = axisDefinition.getMinimum();

        double yScaleFactor = drawingRect.getHeight() / (maxValue - minValue);

        // Compute the pixel y location.  Clip to bounds of rectangle.
        // The distince in pixels frmo the data value to the axis maximum
        double delta = (maxValue - dataY) * yScaleFactor;
        double pY = drawingRect.getY() + delta;

        return (int) Math.max(drawingRect.getMinY(), Math.min(drawingRect.getMaxY(), pY));
    }

    protected Rectangle calculateDrawingRect(Rectangle arect) {

        double buffer = Math.min(arect.getHeight() * marginFraction, 10);
        Rectangle adjustedRect = new Rectangle(arect);
        adjustedRect.y = (int) (arect.getY() + buffer);
        adjustedRect.height = (int) (arect.height - (adjustedRect.y - arect.getY()));


        return adjustedRect;
    }

    public void setMarginFraction(double marginFraction) {
        this.marginFraction = marginFraction;
    }
}
