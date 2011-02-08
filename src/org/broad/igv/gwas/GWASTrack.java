package org.broad.igv.gwas;

import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.renderer.GraphicUtils;
import org.broad.igv.renderer.Renderer;
import org.broad.igv.track.AbstractTrack;
import org.broad.igv.track.RegionScoreType;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.WindowFunction;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.ChromosomeColors;
import org.broad.igv.util.ColorUtilities;
import org.broad.igv.util.ResourceLocator;

import java.awt.*;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

//import org.broad.igv.ui.IGVModel;

/**
 * Created by IntelliJ IDEA.
 * User: jussi
 * Date: Nov 23, 2009
 * Time: 4:56:40 PM
 * To change this template use File | Settings | File Templates.
 */
public class GWASTrack extends AbstractTrack {

    private GWASData gData;
    private static final Logger log = Logger.getLogger(GWASTrack.class);

    private static final int AXIS_AREA_WIDTH = 60;
    protected static Color axisLineColor = new Color(255, 180, 180);
    private double trackMinY;
    private double maxY;
    private double scale;
    private GWASParser parser;

    String getDisplayName() {
        return displayName;
    }

    public void setDisplayName(String displayName) {
        this.displayName = displayName;
    }

    private String displayName = "GWAS Track";


    public GWASTrack(String name) {
        super(name);
    }

    /**
     * Constructor for a new GWAS track
     *
     * @param locator
     * @param id
     * @param name
     * @param gData
     * @param parser
     */
    public GWASTrack(ResourceLocator locator, String id, String name, GWASData gData, GWASParser parser) {


        super(locator, id, name);

        PreferenceManager prefs = PreferenceManager.getInstance();
        super.setHeight(prefs.getAsInt(PreferenceManager.GWAS_TRACK_HEIGHT));

        this.gData = gData;

        // Set range from 0 to highest value rounded to greater integer
        int maxValue = (int) Math.ceil(gData.getMaxValue());
        super.setDataRange(new DataRange(0, (int) (maxValue / 2), maxValue));

        this.parser = parser;

    }

    /**
     * Render GWAS data
     *
     * @param context
     * @param arect
     */
    public void render(RenderContext context, Rectangle arect) {

        PreferenceManager prefs = PreferenceManager.getInstance();

        this.trackMinY = arect.getMinY();
        Rectangle adjustedRect = calculateDrawingRect(arect);
        renderAxis(context, arect);
        this.maxY = adjustedRect.getMaxY();
        this.scale = context.getScale();
        int bufferX = (int) adjustedRect.getMaxX();
        int bufferY = (int) adjustedRect.getMaxY();
        Color[][] drawBuffer = new Color[bufferX + 1][bufferY + 1];
        double origin = context.getOrigin();
        double locScale = context.getScale();


        // Get the Y axis definition, consisting of minimum, maximum, and base value.  Often
        // the base value is == min value which is == 0.

        DataRange axisDefinition = getDataRange();
        float maxValue = axisDefinition.getMaximum();
        float baseValue = axisDefinition.getBaseline();
        float minValue = axisDefinition.getMinimum();

        // Calculate the Y scale factor.
        double yScaleFactor = adjustedRect.getHeight() / (maxValue - minValue);

        // Calculate the Y position in pixels of the base value.
        int baseY = (int) (adjustedRect.getY() + (maxValue - baseValue) * yScaleFactor);

        //int lastPx = 0;
        String chrName = context.getChr();
        List<String> chrList = new ArrayList();
        if (chrName.equals("All")) {
            for (String key : gData.getLocations().keySet()) {
                chrList.add(key);
            }

        } else {
            chrList.add(chrName);

        }
        double dx = Math.ceil(1 / locScale) + 1;
        double rangeMaxValue = Math.ceil(gData.getMaxValue());
        int minPointSize = prefs.getAsInt(PreferenceManager.GWAS_MIN_POINT_SIZE);
        int maxPointSize = prefs.getAsInt(PreferenceManager.GWAS_MAX_POINT_SIZE);


        for (String chr : chrList) {
            if (this.gData.getLocations().containsKey(chr) && this.gData.getValues().containsKey(chr)) {

                int size = this.gData.getLocations().get(chr).size();
                for (int j = 0; j < size; j++) {
                    int start;


                    if (chrName.equals("All"))

                        start = GenomeManager.getInstance().getCurrentGenome().getGenomeCoordinate(chr, this.gData.getLocations().get(chr).get(j));

                    else


                        start = this.gData.getLocations().get(chr).get(j);

                    Color drawColor = ColorUtilities.convertRGBStringToColor(prefs.get(PreferenceManager.GWAS_PRIMARY_COLOR));
                    if (prefs.getAsBoolean(PreferenceManager.GWAS_USE_CHR_COLORS))
                        drawColor = ChromosomeColors.getColor(chr);
                    else if (prefs.getAsBoolean(PreferenceManager.GWAS_ALTERNATING_COLORS)) {

                        Object[] keys = this.gData.getLocations().keySet().toArray();

                        int chrCounter = 0;
                        //int lineCounter = 0;
                        while (chrCounter < keys.length && !keys[chrCounter].toString().equals(chr)) {
                            chrCounter++;

                        }
                        if (chrCounter % 2 == 0)
                            drawColor = ColorUtilities.convertRGBStringToColor(prefs.get(PreferenceManager.GWAS_SECONDARY_COLOR));

                    }
                    // Note -- don't cast these to an int until the range is checked.
                    // could get an overflow.
                    double pX = ((start - origin) / locScale);
                    // Todo -- No need to calculate end for the genomic region?
                    //double dx = Math.ceil((Math.max(1, start - start)) / locScale) + 1;
                    if ((pX + dx < 0)) {
                        continue;
                    } else if (pX > adjustedRect.getMaxX()) {
                        break;
                    }

                    float dataY = this.gData.getValues().get(chr).get(j);

                    if (!Float.isNaN(dataY)) {


                        // Compute the pixel y location.  Clip to bounds of rectangle.
                        int pY = (int)
                                Math.min(adjustedRect.getMaxY(),
                                        adjustedRect.getY()
                                                + (maxValue - dataY) * yScaleFactor);

                        int xPointSize = (int) Math.ceil(dataY / rangeMaxValue * maxPointSize);
                        // Min x size is the width of a single nucleotide
                        int xMinPointSize = (int) (1 / locScale);
                        // If min x size is smaller than minimun point size, use min point size
                        if (xMinPointSize < minPointSize)
                            xMinPointSize = minPointSize;

                        if (xPointSize < xMinPointSize)
                            xPointSize = xMinPointSize;
                        // Scale y size based on the used on the data value and max point size
                        int yPointSize = (int) Math.ceil(dataY / rangeMaxValue * maxPointSize);
                        if (yPointSize < minPointSize)
                            yPointSize = minPointSize;

                        int tmpX = (int) pX;
                        if (tmpX < 0)
                            tmpX = 0;
                        // Calculate offsets that can be used to center larger data points around their location
                        int xCenterOffset = (int) xPointSize / 2;
                        int yCenterOffset = (int) yPointSize / 2;

                        // Loop through all the pixels of the data point
                        for (int drawX = 0; drawX < xPointSize; drawX++) {
                            for (int drawY = 0; drawY < yPointSize; drawY++) {

                                int tmpDrawX = tmpX + drawX - xCenterOffset;
                                int tmpDrawY = pY + drawY - yCenterOffset;
                                // Check if pixels are inside array, if so, mark them in the buffer
                                if (tmpDrawX >= 0 && tmpDrawX <= adjustedRect.getMaxX() && tmpDrawY >= 0 && tmpDrawY <= adjustedRect.getMaxY()) {


                                    drawBuffer[tmpDrawX][tmpDrawY] = drawColor;
                                }
                            }
                        }
                    }

                }
            }


        }
        int tmpDx = (int) Math.ceil(1 / locScale) + 1;
        for (int x = 0; x < bufferX; x++) {
            for (int y = 0; y < bufferY; y++) {
                if (drawBuffer[x][y] != null)
                    drawDataPoint(drawBuffer[x][y], tmpDx, x, baseY, y, context);
            }


        }

    }


    void drawDataPoint(Color graphColor, int dx, int pX, int baseY, int pY,
                       RenderContext context) {
        context.getGraphic2DForColor(graphColor).fillRect(pX, pY, 1, 1);

    }

    private static final DecimalFormat formatter = new DecimalFormat();


    void renderAxis(RenderContext context, Rectangle arect) {


        Rectangle drawingRect = calculateDrawingRect(arect);

        Color labelColor = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.CHART_COLOR_TRACK_NAME) ? getColor() : Color.black;
        Graphics2D labelGraphics = context.getGraphic2DForColor(labelColor);
        labelGraphics.setFont(FontManager.getScalableFont(8));

        if (PreferenceManager.getInstance().getAsBoolean(PreferenceManager.CHART_DRAW_TRACK_NAME)) {
            // Only attempt if track height is > 25 pixels
            if (arect.getHeight() > 25) {
                Rectangle labelRect = new Rectangle(arect.x, arect.y + 10, arect.width, 10);
                labelGraphics.setFont(FontManager.getScalableFont(10));
                GraphicUtils.drawCenteredText(getDisplayName(), labelRect, labelGraphics);
            }
        }

        if (PreferenceManager.getInstance().getAsBoolean(PreferenceManager.CHART_DRAW_Y_AXIS)) {
            Rectangle axisRect = new Rectangle(arect.x, arect.y + 1, AXIS_AREA_WIDTH, arect.height);
            DataRange axisDefinition = getDataRange();
            float maxValue = axisDefinition.getMaximum();
            float minValue = axisDefinition.getMinimum();


            // Bottom (minimum tick mark)
            int pY = computeYPixelValue(drawingRect, axisDefinition, minValue);

            labelGraphics.drawLine(axisRect.x + AXIS_AREA_WIDTH - 10, pY,
                    axisRect.x + AXIS_AREA_WIDTH - 5, pY);
            GraphicUtils.drawRightJustifiedText(formatter.format(minValue),
                    axisRect.x + AXIS_AREA_WIDTH - 15, pY, labelGraphics);

            // Top (maximum tick mark)
            int topPY = computeYPixelValue(drawingRect, axisDefinition, maxValue);

            labelGraphics.drawLine(axisRect.x + AXIS_AREA_WIDTH - 10, topPY,
                    axisRect.x + AXIS_AREA_WIDTH - 5, topPY);
            GraphicUtils.drawRightJustifiedText(formatter.format(maxValue),
                    axisRect.x + AXIS_AREA_WIDTH - 15, topPY + 4, labelGraphics);

            // Connect top and bottom
            labelGraphics.drawLine(axisRect.x + AXIS_AREA_WIDTH - 10, topPY,
                    axisRect.x + AXIS_AREA_WIDTH - 10, pY);

            // Draw middle tick marks from top to bottom, draw only if space available
            int lastY = pY;
            for (int i = (int) minValue + 1; i < (int) maxValue; i++) {

                int midPY = computeYPixelValue(drawingRect, axisDefinition, i);
                if ((midPY < lastY - 15) && (midPY < pY - 15) && (midPY > topPY + 15)) {
                    labelGraphics.drawLine(axisRect.x + AXIS_AREA_WIDTH - 10, midPY,
                            axisRect.x + AXIS_AREA_WIDTH - 5, midPY);
                    GraphicUtils.drawRightJustifiedText(formatter.format(i),
                            axisRect.x + AXIS_AREA_WIDTH - 15, midPY + 4, labelGraphics);
                    lastY = midPY;
                }
            }
        }
    }


    int computeYPixelValue(Rectangle drawingRect, DataRange axisDefinition, double dataY) {

        double maxValue = axisDefinition.getMaximum();
        double minValue = axisDefinition.getMinimum();

        double yScaleFactor = drawingRect.getHeight() / (maxValue - minValue);

        // Compute the pixel y location.  Clip to bounds of rectangle.
        // The distance in pixels from the data value to the axis maximum
        double delta = (maxValue - dataY) * yScaleFactor;
        double pY = drawingRect.getY() + delta;

        return (int) Math.max(drawingRect.getMinY(), Math.min(drawingRect.getMaxY(), pY));
    }


    Rectangle calculateDrawingRect(Rectangle arect) {

        double buffer = Math.min(arect.getHeight() * 0.2, 10);
        Rectangle adjustedRect = new Rectangle(arect);
        adjustedRect.y = (int) (arect.getY() + buffer);
        adjustedRect.height = (int) (arect.height - (adjustedRect.y - arect.getY()));


        return adjustedRect;
    }


    /**
     * Find index of data point closest to given chromosomal location and y-coordinate
     *
     * @param chr
     * @param y
     * @param location
     * @param maxDistance
     * @return
     */
    public int findIndex(String chr, int y, int location, int maxDistance) {


        // Calculate offset from track location by other tracks
        y = y - (int) this.trackMinY;

        // Estimate values near y coordinate
        double valueEstimate = (1 - y / this.maxY) * this.getDataRange().getMaximum();
        // Percentage threshold for searching values on y-axis
        double threshold = 0.1;

        // Based on threshold, set max and min y values to search for values
        double topValue = valueEstimate + threshold * this.getDataRange().getMaximum();
        double bottomValue = valueEstimate - threshold * this.getDataRange().getMaximum();

        if (bottomValue < 0)
            bottomValue = 0;


        // Find data point based on the given coordinates and search parameters
        return this.gData.getNearestIndexByLocation(chr, location, bottomValue, topValue, maxDistance);

    }

    /**
     * Get description for a data point at given chromosome and index. If description is not cached, re-populate cache.
     *
     * @param chr
     * @param index
     * @return
     */

    public String getDescription(String chr, int index) {

        String textValue = "";

        float value = this.gData.getValues().get(chr).get(index);
        int hitLocation = this.gData.getLocations().get(chr).get(index);
        int rowIndex = gData.getCumulativeChrLocation(chr) + index;

        textValue += chr + ": " + hitLocation + "<br>";
        textValue += "Value: " + value + "<br>";
        textValue += "-----<br>";

        try {

            // Look for data point description from cache
            String tmpDescription = gData.getDescriptionCache().getDescriptionString(chr, hitLocation);

            // If no description found, populate cache with the description
            if (tmpDescription == null) {
                // Calculate starting row based on the cache size, i.e. cache descriptions before and after estimated hit location
                int tmpRow = rowIndex - (gData.getDescriptionCache().getMaxSize() / 2);
                if (tmpRow < 0)
                    tmpRow = 0;

                this.gData = parser.parseDescriptions(gData, chr, hitLocation, tmpRow);
                tmpDescription = gData.getDescriptionCache().getDescriptionString(chr, hitLocation);

            }

            // Add fetched description
            textValue += tmpDescription;


        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
        return textValue;

    }


    public String getValueStringAt(String chr, double position, int y, ReferenceFrame frame) {


        String textValue = "";

        int location = (int) position;

        // Set maximum search distance to be the amount of nucleotides corresponding to 2 pixels on the screen
        int maxDistance = (int) (this.scale) * 2;

        // Convert from All view chr coordinates
        if (chr.equals("All")) {

            Genome.ChromosomeCoordinate chrCoordinate = GenomeManager.getInstance().getCurrentGenome().getChromosomeCoordinate(location);
            chr = chrCoordinate.getChr();
            location = chrCoordinate.getCoordinate();
            maxDistance = maxDistance * 1000;
        }

        int index = findIndex(chr, y, location, maxDistance);

        // If there is a data point at the given location, fetch description
        if (index != -1)
            textValue += getDescription(chr, index);

        return textValue;
    }

    // Needed to compile

    public GWASTrack(ResourceLocator dataResourceLocator, String id, String name) {
        super(dataResourceLocator, id, name);
    }

    public void setStatType(WindowFunction type) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    public WindowFunction getWindowFunction() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public void setRendererClass(Class rc) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    public Renderer getRenderer() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public boolean isLogNormalized() {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }


    public float getRegionScore(String chr, int start, int end, int zoom, RegionScoreType type) {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
