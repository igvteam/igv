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
    private static Logger log = Logger.getLogger(GWASTrack.class);

    private static final int AXIS_AREA_WIDTH = 60;
    protected static Color axisLineColor = new Color(255, 180, 180);
    //  private RenderContext latestContext;
    private double maxY;
    private double scale;
    private GWASParser parser;
    private GWASPreferences preferences = new GWASPreferences();

    public String getDisplayName() {
        return displayName;
    }

    public void setDisplayName(String displayName) {
        this.displayName = displayName;
    }

    private String displayName = "GWAS Track";


    public GWASTrack(String name) {
        super(name);
    }

    public GWASTrack(ResourceLocator locator, String id, String name, GWASData gData, GWASParser parser) {
        super(locator, id, name);
        super.setHeight(140);
        //log.info("DATA RANGE: " + this.getDataRange().getMinimum() + ".." + this.getDataRange().getBaseline() + ".." +this.getDataRange().getMaximum());

        this.gData = gData;
        int maxValue = (int) Math.ceil(gData.getMaxValue());

        // Set range from 0 to highest value rounded to greater integer
        super.setDataRange(new DataRange(0, (int) maxValue / 2, maxValue));
        //log.info("DATA RANGE: " + this.getDataRange().getMinimum() + ".." + this.getDataRange().getBaseline() + ".." +this.getDataRange().getMaximum());
        this.parser = parser;

    }

    public void render(RenderContext context, Rectangle arect) {
        //log.info("GWASTrack render called");
        //log.info("context chr: " + context.getChr());
        // this.latestContext = context;


        //Graphics2D noDataGraphics = context.getGraphic2DForColor(IGVConstants.NO_DATA_COLOR);


        Rectangle adjustedRect = calculateDrawingRect(arect);
        renderAxis(context, arect);

        this.maxY = adjustedRect.getMaxY();
        this.scale = context.getScale();

        int bufferX = (int) adjustedRect.getMaxX();
        int bufferY = (int) adjustedRect.getMaxY();

        Color[][] drawBuffer = new Color[bufferX + 1][bufferY + 1];


        //log.info("maxY: " + adjustedRect.getMaxY() + " Y: " + adjustedRect.getY() + " maxx: " + adjustedRect.getMaxX() + " X: " + adjustedRect.getX());

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
        int minPointSize = this.preferences.getMinPointSize();
        int maxPointSize = this.preferences.getMaxPointSize();


        for (String chr : chrList) {
            if (this.gData.getLocations().containsKey(chr) && this.gData.getValues().containsKey(chr)) {

                int size = this.gData.getLocations().get(chr).size();
                for (int j = 0; j < size; j++) {
                    //log.info(locationList.get(j));
                    int start;


                    if (chrName.equals("All"))

                        //start = IGVModel.getInstance().getViewContext().getGenome().getGenomeCoordinate(chr, this.gData.getLocations().get(chr).get(j));
                        start = GenomeManager.getInstance().getCurrentGenome().getGenomeCoordinate(chr, this.gData.getLocations().get(chr).get(j));

                    else


                        start = this.gData.getLocations().get(chr).get(j);

                    //log.info("start: " + start);
                    Color drawColor = this.preferences.getPrimaryColor();
                    if (this.preferences.isChrColors())
                        drawColor = ChromosomeColors.getColor(chr);
                    else if (this.preferences.isAlternatingColors()) {

                        Object[] keys = this.gData.getLocations().keySet().toArray();

                        int chrCounter = 0;
                        int lineCounter = 0;
                        while (chrCounter < keys.length && !keys[chrCounter].toString().equals(chr)) {
                            chrCounter++;

                        }
                        if (chrCounter % 2 == 0)
                            drawColor = this.preferences.getSecondaryColor();

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

                        // Commented out because it draws of the chart y-values
                        // Compute the pixel y location.  Clip to bounds of rectangle.
                        /*
                       int pY = (int) Math.max(adjustedRect.getMinY(),
                               Math.min(adjustedRect.getMaxY(),
                                       adjustedRect.getY()
                                               + (maxValue - dataY) * yScaleFactor));
                        */
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
                    //drawDataPoint(color, (int) dx, (int) pX, baseY, pY, context);
                    //log.info(color + " " + " " + (int) dx  + " " +  (int) pX  + " " +  baseY  + " " +  pY  + " " );

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
        //context.getGraphic2DForColor(graphColor).fillRect(pX, pY, dx, 2);
        context.getGraphic2DForColor(graphColor).fillRect(pX, pY, 1, 1);

    }

    private static final DecimalFormat formatter = new DecimalFormat();


    public void renderAxis(RenderContext context, Rectangle arect) {

        // For now disable axes for all chromosome view
        /*
           if (context.getChr().equals(IGVConstants.CHR_ALL)) {
               return;
           }
        */
        //super.renderAxis(track, context, arect);

        Rectangle drawingRect = calculateDrawingRect(arect);

        //PreferenceManager.ChartPreferences prefs = PreferenceManager.getInstance().getChartPreferences();

        // Color labelColor = prefs.isColorTrackName() ? getColor() : Color.black;

        Color labelColor = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.CHART_COLOR_TRACK_NAME) ? getColor() : Color.black;


        Graphics2D labelGraphics = context.getGraphic2DForColor(labelColor);

        labelGraphics.setFont(FontManager.getScalableFont(8));

        //if (prefs.isDrawTrackName()) {
        if (PreferenceManager.getInstance().getAsBoolean(PreferenceManager.CHART_DRAW_TRACK_NAME)) {
            // Only attempt if track height is > 25 pixels
            if (arect.getHeight() > 25) {
                Rectangle labelRect = new Rectangle(arect.x, arect.y + 10, arect.width, 10);
                labelGraphics.setFont(FontManager.getScalableFont(10));
                GraphicUtils.drawCenteredText(getDisplayName(), labelRect, labelGraphics);
            }
        }

        //if (prefs.isDrawAxis()) {

        if (PreferenceManager.getInstance().getAsBoolean(PreferenceManager.CHART_DRAW_Y_AXIS)) {

            Rectangle axisRect = new Rectangle(arect.x, arect.y + 1, AXIS_AREA_WIDTH, arect.height);


            DataRange axisDefinition = getDataRange();
            float maxValue = axisDefinition.getMaximum();
            //float baseValue = axisDefinition.getBaseline();
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


    /**
     * Get a grapphics object for the baseline.
     * TODO -- make the line style settable by the user
     *
     * @param context
     * @return
     */
    private static Graphics2D getBaselineGraphics(RenderContext context) {
        Graphics2D baselineGraphics;
        Stroke thindashed = new BasicStroke(1.0f,    // line width

                /* cap style */
                BasicStroke.CAP_BUTT,

                /* join style, miter limit */
                BasicStroke.JOIN_BEVEL, 1.0f,

                /* the dash pattern */
                new float[]{8.0f, 3.0f, 2.0f, 3.0f},

                /* the dash phase */
                0.0f);    /* on 8, off 3, on 2, off 3 */

        baselineGraphics = (Graphics2D) context.getGraphic2DForColor(Color.lightGray).create();

        // baselineGraphics.setStroke(thindashed);
        return baselineGraphics;
    }


    int computeYPixelValue(Rectangle drawingRect, DataRange axisDefinition, double dataY) {

        double maxValue = axisDefinition.getMaximum();
        double minValue = axisDefinition.getMinimum();

        double yScaleFactor = drawingRect.getHeight() / (maxValue - minValue);

        // Compute the pixel y location.  Clip to bounds of rectangle.
        // The distince in pixels frmo the data value to the axis maximum
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


    public String getValueStringAt(String chr, double position, int y, ReferenceFrame frame) {


        String textValue = "";

        // Estimate values near y coordinate
        double valueEstimate = (1 - y / this.maxY) * this.getDataRange().getMaximum();
        // Percentage threshold for searching values on y-axis
        double threshold = 0.1;

        // Based on threshold, set max and min y values to search for values
        double topValue = valueEstimate + threshold * this.getDataRange().getMaximum();
        double bottomValue = valueEstimate - threshold * this.getDataRange().getMaximum();

        if (bottomValue < 0)
            bottomValue = 0;

        // Set maximum search distance to be the amount of nucleotides corresponding to 2 pixels on the screen
        int maxDistance = (int) (this.scale) * 2;

        int location = (int) position;

        // Convert from All view chr coordinates
        if (chr.equals("All")) {

            Genome.ChromosomeCoordinate chrCoordinate = GenomeManager.getInstance().getCurrentGenome().getChromosomeCoordinate(location);
            chr = chrCoordinate.getChr();
            location = chrCoordinate.getCoordinate();
            maxDistance = maxDistance * 1000;
        }

        // Find data point based on the given coordinates and search parameters
        int index = gData.getNearestIndexByLocation(chr, location, bottomValue, topValue, maxDistance);

        float value = -1;
        int hitLocation = -1;
        int rowIndex = -1;
        if (index != -1) {

            value = this.gData.getValues().get(chr).get(index);
            hitLocation = this.gData.getLocations().get(chr).get(index);
            rowIndex = gData.getCumulativeChrLocation(chr) + index;

            textValue += chr + ": " + hitLocation + "<br>";
            textValue += "Value: " + value + "<br>";
            textValue += "-----<br>";

            try {


                String tmpDescription = gData.getDescriptionCache().getDescriptionString(chr, hitLocation);
                if (tmpDescription == null) {
                    // Calculate starting row based on the cache size, i.e. cache descriptions before and after estimated hit location
                    int tmpRow = rowIndex - (gData.getDescriptionCache().getMaxSize() / 2);
                    if (rowIndex < 0) {
                        rowIndex = 0;
                    }
                    this.gData = parser.parseDescriptions(gData, chr, hitLocation, tmpRow);
                    tmpDescription = gData.getDescriptionCache().getDescriptionString(chr, hitLocation);

                }
                textValue += tmpDescription;


            } catch (IOException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }


        }

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
