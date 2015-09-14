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

package org.broad.igv.gwas;

import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.genome.ChromosomeCoordinate;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.renderer.GraphicUtils;
import org.broad.igv.session.IGVSessionReader;
import org.broad.igv.session.SessionXmlAdapters;
import org.broad.igv.session.SubtlyImportant;
import org.broad.igv.track.*;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.ui.panel.IGVPopupMenu;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.util.UIUtilities;
import org.broad.igv.util.ChromosomeColors;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.collections.DoubleArrayList;
import org.broad.igv.util.collections.IntArrayList;

import javax.swing.*;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlType;
import javax.xml.bind.annotation.adapters.XmlJavaTypeAdapter;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;

/**
 * @author  jussi
 * @since  Nov 23, 2009
 */
@XmlType(factoryMethod = "getNextTrack")
public class GWASTrack extends AbstractTrack {

    // Color properties
    @XmlAttribute private int minPointSize;
    @XmlAttribute private int maxPointSize;

    @XmlAttribute private boolean useChrColors;
    @XmlAttribute private boolean singleColor;
    @XmlAttribute private boolean alternatingColors;

    @XmlJavaTypeAdapter(SessionXmlAdapters.Color.class)
    @XmlAttribute private Color primaryColor;
    @XmlJavaTypeAdapter(SessionXmlAdapters.Color.class)
    @XmlAttribute private Color secondaryColor;

    private GWASData gData;
    private static final Logger log = Logger.getLogger(GWASTrack.class);

    private static final int AXIS_AREA_WIDTH = 60;
    @XmlAttribute private double trackMinY;
    @XmlAttribute private double maxY;
    @XmlAttribute private double scale;
    private GWASParser parser;
    private static final DecimalFormat formatter = new DecimalFormat();
    //private String displayName = "GWAS Track";
    @XmlAttribute private String displayName = null;
    @XmlAttribute private boolean drawYAxis = true;
    private boolean showAxis = true;

    String getDisplayName() {
        return displayName;
    }

    public void setDisplayName(String displayName) {
        this.displayName = displayName;
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

        // Set range from 0 to highest value rounded to greater integer
        int maxValue = (int) Math.ceil(gData.getMaxValue());
        super.setDataRange(new DataRange(0, (maxValue / 2), maxValue));


        // Get default values
        super.setHeight(prefs.getAsInt(PreferenceManager.GWAS_TRACK_HEIGHT));

        this.minPointSize = prefs.getAsInt(PreferenceManager.GWAS_MIN_POINT_SIZE);
        this.maxPointSize = prefs.getAsInt(PreferenceManager.GWAS_MAX_POINT_SIZE);

        this.primaryColor = ColorUtilities.stringToColor(prefs.get(PreferenceManager.GWAS_PRIMARY_COLOR));
        this.secondaryColor = ColorUtilities.stringToColor(prefs.get(PreferenceManager.GWAS_SECONDARY_COLOR));
        this.singleColor = prefs.getAsBoolean(PreferenceManager.GWAS_SINGLE_COLOR);
        this.alternatingColors = prefs.getAsBoolean(PreferenceManager.GWAS_ALTERNATING_COLORS);
        this.useChrColors = prefs.getAsBoolean(PreferenceManager.GWAS_USE_CHR_COLORS);
        this.showAxis = prefs.getAsBoolean(PreferenceManager.GWAS_SHOW_AXIS);

        this.gData = gData;
        this.parser = parser;


    }


    public void render(RenderContext context, Rectangle arect) {

        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        this.trackMinY = arect.getMinY();
        Rectangle adjustedRect = calculateDrawingRect(arect);
        double adjustedRectMaxX = adjustedRect.getMaxX();
        double adjustedRectMaxY = adjustedRect.getMaxY();
        double adjustedRectY = adjustedRect.getY();
        this.maxY = adjustedRectMaxY;
        this.scale = context.getScale();
        int bufferX = (int) adjustedRectMaxX;
        int bufferY = (int) adjustedRectMaxY;
        Color[][] drawBuffer = new Color[bufferX + 1][bufferY + 1];
        double origin = context.getOrigin();
        double locScale = context.getScale();

        // Get the Y axis definition, consisting of minimum, maximum, and base value.  Often
        // the base value is == min value which is == 0.

        DataRange axisDefinition = this.getDataRange();
        float maxValue = axisDefinition.getMaximum();
        float minValue = axisDefinition.getMinimum();

        // Calculate the Y scale factor.
        double yScaleFactor = adjustedRect.getHeight() / (maxValue - minValue);

        //int lastPx = 0;
        String chrName = context.getChr();
        ArrayList<String> chrList = new ArrayList();
        if (chrName.equals("All")) {
            for (String key : gData.getLocations().keySet()) {
                chrList.add(key);
            }

        } else {
            chrList.add(chrName);

        }
        double dx = Math.ceil(1 / locScale) + 1;
        double rangeMaxValue = Math.ceil(gData.getMaxValue());

        double pointSizeScale = rangeMaxValue / maxPointSize;

        Color drawColor = this.primaryColor;
        Object[] chrs = this.gData.getLocations().keySet().toArray();

        int xMinPointSize = (int) (1 / locScale);

        // Loop through data points, chromosome by chromosome

        for (String chr : chrList) {
            if (this.gData.getLocations().containsKey(chr) && this.gData.getValues().containsKey(chr)) {


                // Choose a color for the chromosome

                // Use specific color for each chromosome
                if (this.useChrColors)
                    drawColor = ChromosomeColors.getColor(chr);

                    // Use two alternating colors for chromosomes
                else if (this.alternatingColors) {
                    int chrCounter = 0;
                    while (chrCounter < chrs.length && !chrs[chrCounter].toString().equals(chr))
                        chrCounter++;

                    if (chrCounter % 2 == 0)
                        drawColor = this.secondaryColor;
                    else
                        drawColor = this.primaryColor;

                }

                IntArrayList locations = this.gData.getLocations().get(chr);
                DoubleArrayList values = this.gData.getValues().get(chr);

                int size = locations.size();

                // Loop through data points in a chromosome
                for (int j = 0; j < size; j++) {

                    // Get location, e.g. start for the data point
                    int start;
                    if (chrName.equals("All"))
                        start = genome.getGenomeCoordinate(chr, locations.get(j));
                    else
                        start = locations.get(j);

                    // Based on location, calculate X-coordinate, or break if outside of the view
                    double pX = ((start - origin) / locScale);
                    if ((pX + dx < 0))
                        continue;
                    else if (pX > adjustedRectMaxX)
                        break;

                    // Based on value of the data point, calculate Y-coordinate
                    double dataY = values.get(j);

                    if (!Double.isNaN(dataY)) {

                        int xPointSize = (int) Math.ceil(dataY / pointSizeScale);

                        // Scale y size based on the used range, data value and max point size
                        int yPointSize = xPointSize;
                        if (yPointSize < minPointSize)
                            yPointSize = minPointSize;

                        // If x minimum size is smaller than point minimum size, use minimum point size
                        if (xMinPointSize < minPointSize)
                            xMinPointSize = minPointSize;

                        if (xPointSize < xMinPointSize)
                            xPointSize = xMinPointSize;

                        // Point sizes divided by two to center locations of large points
                        int x = (int) pX - (xPointSize / 2);
                        int y = ((int) Math.min(adjustedRectMaxY, adjustedRectY + (maxValue - dataY) * yScaleFactor)) - (yPointSize / 2);

                        int maxDrawX = x + xPointSize;
                        int maxDrawY = y + yPointSize;
                        if (x < 0)
                            x = 0;
                        if (y < 0)
                            y = 0;

                        if (maxDrawX > adjustedRectMaxX)
                            maxDrawX = (int) adjustedRectMaxX;
                        if (maxDrawY > adjustedRectMaxY)
                            maxDrawY = (int) adjustedRectMaxY;

                        // Loop through all the pixels of the data point and fill in drawing buffer
                        for (int drawX = x; drawX < maxDrawX; drawX++)
                            for (int drawY = y; drawY < maxDrawY; drawY++)
                                drawBuffer[drawX][drawY] = drawColor;
                    }
                }
            }
        }

        // Draw the pixels from the drawing buffer to the canvas

        Graphics2D g = context.getGraphics();
        Color color;
        Color prevColor = null;

        for (int x = 0; x < bufferX; x++)
            for (int y = 0; y < bufferY; y++) {
                color = drawBuffer[x][y];
                if (color != null) {
                    if (!color.equals(prevColor)) {
                        g.setColor(color);
                        prevColor = color;
                    }
                    g.fillRect(x, y, 1, 1);
                }
            }

        // Draw the legend axis
        if (showAxis) {
            this.renderAxis(context, arect);
        }

    }


    void renderAxis(RenderContext context, Rectangle arect) {

        Rectangle drawingRect = calculateDrawingRect(arect);

        Color labelColor = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.CHART_COLOR_TRACK_NAME) ? getColor() : Color.black;
        Graphics2D labelGraphics = context.getGraphic2DForColor(labelColor);
        //Graphics2D labelGraphics = context.getGraphic2DForColor(Color.black);

        labelGraphics.setFont(FontManager.getFont(8));

        String tmpDisplayName = this.getDisplayName();
        if (tmpDisplayName != null && tmpDisplayName.length() > 0 && PreferenceManager.getInstance().getAsBoolean(PreferenceManager.CHART_DRAW_TRACK_NAME)) {
            // Only attempt if track height is > 25 pixels
            if (arect.getHeight() > 25) {
                Rectangle labelRect = new Rectangle(arect.x, arect.y + 10, arect.width, 10);
                labelGraphics.setFont(FontManager.getFont(10));
                GraphicUtils.drawCenteredText(tmpDisplayName, labelRect, labelGraphics);
            }
        }

        //if (PreferenceManager.getInstance().getAsBoolean(PreferenceManager.CHART_DRAW_Y_AXIS)) {
        if (this.drawYAxis) {
            labelGraphics = context.getGraphic2DForColor(Color.black);
            labelGraphics.setFont(FontManager.getFont(11));


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
    int findIndex(String chr, int y, int location, int maxDistance) {


        // Calculate offset from track location by other tracks
        y = y - (int) this.trackMinY;
        double tmpMaxY = this.maxY - this.trackMinY;

        // Estimate values near y coordinate
        double valueEstimate = (1 - y / tmpMaxY) * this.getDataRange().getMaximum();
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

    String getDescription(String chr, int index) {

        String textValue = "";

        double value = this.gData.getValues().get(chr).get(index);
        int hitLocation = this.gData.getLocations().get(chr).get(index);
        int rowIndex = gData.getCumulativeChrLocation(chr) + index;

        textValue += chr + ": " + hitLocation + "<br>";
        textValue += "Value: " + value + "<br>";
        textValue += "-----<br>";

        try {

            // Look for data point description from cache
            //String tmpDescription = gData.getDescriptionCache().getDescriptionString(chr, hitLocation);
            String tmpDescription = gData.getDescriptionCache().getDescriptionString(chr, hitLocation, value);


            // If no description found, populate cache with the description
            if (tmpDescription == null) {
                // Calculate starting row based on the cache size, i.e. cache descriptions before and after estimated hit location
                int tmpRow = rowIndex - (gData.getDescriptionCache().getMaxSize() / 2);
                if (tmpRow < 0)
                    tmpRow = 0;

                this.gData = parser.parseDescriptions(gData, chr, hitLocation, tmpRow);
                tmpDescription = gData.getDescriptionCache().getDescriptionString(chr, hitLocation, value);

            }

            // Add fetched description
            textValue += tmpDescription;


        } catch (IOException e) {
            log.error(e.getMessage(), e);
        }
        return textValue;

    }


    public String getValueStringAt(String chr, double position, int y, ReferenceFrame frame) {

        int location = (int) position;

        // Set maximum search distance to be the amount of nucleotides corresponding to 2 pixels on the screen
        int maxDistance = (int) (this.scale) * 2;

        // Convert from All view chr coordinates
        if (chr.equals("All")) {

            ChromosomeCoordinate chrCoordinate = GenomeManager.getInstance().getCurrentGenome().getChromosomeCoordinate(location);
            chr = chrCoordinate.getChr();
            location = chrCoordinate.getCoordinate();
            maxDistance = maxDistance * 1000;
        }

        int index = findIndex(chr, y, location, maxDistance);

        // If there is a data point at the given location, fetch description
        return index >= 0 ? getDescription(chr, index) : null;

    }

    /**
     * Override to return a specialized popup menu
     *
     * @return
     */
    @Override
    public IGVPopupMenu getPopupMenu(TrackClickEvent te) {

        IGVPopupMenu popupMenu = new IGVPopupMenu();

        JLabel popupTitle = new JLabel("  " + getName(), JLabel.CENTER);

        Font newFont = popupMenu.getFont().deriveFont(Font.BOLD, 12);
        popupTitle.setFont(newFont);
        popupMenu.add(popupTitle);

        Collection<Track> tmpTracks = new ArrayList();
        tmpTracks.add(this);


        popupMenu.addSeparator();


        popupMenu.add(TrackMenuUtils.getTrackRenameItem(tmpTracks));
        popupMenu.add(TrackMenuUtils.getRemoveMenuItem(tmpTracks));
        popupMenu.addSeparator();

        popupMenu.add(TrackMenuUtils.getDataRangeItem(tmpTracks));
        popupMenu.add(TrackMenuUtils.getChangeTrackHeightItem(tmpTracks));
        popupMenu.addSeparator();

        popupMenu.add(new JLabel("Color scheme", JLabel.LEFT));

        popupMenu.add(addChrColorItem(popupMenu));
        popupMenu.add(addSingleColorItem(popupMenu));
        popupMenu.add(addAlternatingColorItem(popupMenu));
        popupMenu.addSeparator();

        popupMenu.add(addPrimaryColorItem(popupMenu));
        popupMenu.add(addSecondaryColorItem(popupMenu));

        popupMenu.addSeparator();
        popupMenu.add(addMinPointSizeItem(popupMenu));
        popupMenu.add(addMaxPointSizeItem(popupMenu));

        popupMenu.addSeparator();
        addShowAxisItem(popupMenu);

        return popupMenu;
    }


    /**
     * @param menu
     * @return
     */
    public JMenuItem addShowAxisItem(JPopupMenu menu) {
        final JCheckBoxMenuItem axisItem = new JCheckBoxMenuItem("Show axis");
        axisItem.setSelected(showAxis);
        axisItem.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                showAxis = axisItem.isSelected();
                PreferenceManager.getInstance().put(PreferenceManager.GWAS_SHOW_AXIS, String.valueOf(showAxis));
                IGV.getInstance().repaintDataPanels();

            }
        });
        menu.add(axisItem);
        return axisItem;
    }

    public JMenuItem addChrColorItem(JPopupMenu menu) {
        JMenuItem colorItem = new JCheckBoxMenuItem("Chromosome color", useChrColors);
        colorItem.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                singleColor = false;
                useChrColors = true;
                alternatingColors = false;
                updateColorPreferences();
                IGV.getInstance().repaintDataPanels();

            }
        });
        menu.add(colorItem);

        return colorItem;
    }

    public JMenuItem addAlternatingColorItem(JPopupMenu menu) {
        JMenuItem colorItem = new JCheckBoxMenuItem("Alternating color", alternatingColors);
        colorItem.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                singleColor = false;
                useChrColors = false;
                alternatingColors = true;
                updateColorPreferences();
                IGV.getInstance().repaintDataPanels();
            }
        });
        menu.add(colorItem);

        return colorItem;
    }

    public JMenuItem addSingleColorItem(JPopupMenu menu) {
        JMenuItem colorItem = new JCheckBoxMenuItem("Single color", singleColor);

        colorItem.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                singleColor = true;
                useChrColors = false;
                alternatingColors = false;
                updateColorPreferences();
                IGV.getInstance().repaintDataPanels();

            }
        });
        menu.add(colorItem);

        return colorItem;
    }


    private void updateColorPreferences() {
        PreferenceManager.getInstance().put(PreferenceManager.GWAS_SINGLE_COLOR, String.valueOf(singleColor));
        PreferenceManager.getInstance().put(PreferenceManager.GWAS_USE_CHR_COLORS, String.valueOf(useChrColors));
        PreferenceManager.getInstance().put(PreferenceManager.GWAS_ALTERNATING_COLORS, String.valueOf(alternatingColors));
    }


    public JMenuItem addPrimaryColorItem(JPopupMenu menu) {
        JMenuItem colorItem = new JMenuItem("Set primary color...");
        colorItem.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                Color color = UIUtilities.showColorChooserDialog("Set primary color", primaryColor);
                if (color != null) {
                    primaryColor = color;
                    String colorString = ColorUtilities.colorToString(primaryColor);
                    PreferenceManager.getInstance().put(PreferenceManager.GWAS_PRIMARY_COLOR, colorString);
                    IGV.getInstance().repaintDataPanels();
                }
            }
        });
        menu.add(colorItem);
        return colorItem;
    }


    public JMenuItem addSecondaryColorItem(JPopupMenu menu) {
        JMenuItem colorItem = new JMenuItem("Set alternating color...");
        colorItem.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                Color color = UIUtilities.showColorChooserDialog("Set alternating color", secondaryColor);
                if (color != null) {
                    secondaryColor = color;
                    String colorString = ColorUtilities.colorToString(secondaryColor);
                    PreferenceManager.getInstance().put(PreferenceManager.GWAS_SECONDARY_COLOR, colorString);
                    IGV.getInstance().repaintDataPanels();
                }

            }
        });
        menu.add(colorItem);

        return colorItem;
    }


    public JMenuItem addMinPointSizeItem(JPopupMenu menu) {
        JMenuItem menuItem = new JMenuItem("Set minimum point size...");
        menuItem.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                changeMinPointSizeValue();
            }
        });
        menu.add(menuItem);
        return menuItem;
    }

    public JMenuItem addMaxPointSizeItem(JPopupMenu menu) {
        JMenuItem menuItem = new JMenuItem("Set maximum point size...");
        menuItem.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                changeMaxPointSizeValue();
            }
        });
        menu.add(menuItem);
        return menuItem;
    }

    private void changeMinPointSizeValue() {

        int value = this.minPointSize;
        String tmpValue = JOptionPane.showInputDialog(
                IGV.getMainFrame(), "Minimum point size in pixels (1-20):", String.valueOf(value));
        if (!(tmpValue == null) || !tmpValue.trim().equals("")) {

            try {
                value = Integer.parseInt(tmpValue);
                if (value < 1 || value > 20) {
                    JOptionPane.showMessageDialog(IGV.getMainFrame(),
                            "Minimum point size must be an integer number between 1 and 20.");
                } else {
                    if (value > this.maxPointSize) {
                        this.maxPointSize = value;
                    }
                    this.minPointSize = value;
                    updatePointSizePreferences();
                    IGV.getInstance().repaintDataPanels();
                }

            } catch (NumberFormatException numberFormatException) {
                JOptionPane.showMessageDialog(IGV.getMainFrame(),
                        "Point size must be an integer number.");
            }
        }
    }

    private void changeMaxPointSizeValue() {

        int value = this.maxPointSize;
        String tmpValue = JOptionPane.showInputDialog(
                IGV.getMainFrame(), "Maximum point size in pixels (1-20):",
                String.valueOf(value));

        if (!(tmpValue == null) || !tmpValue.trim().equals("")) {

            try {
                value = Integer.parseInt(tmpValue);
                if (value < 1 || value > 20) {
                    JOptionPane.showMessageDialog(IGV.getMainFrame(),
                            "Maximum point size must be an integer number between 1 and 20.");
                } else {
                    if (value < this.minPointSize) {
                        this.minPointSize = value;
                    }
                    this.maxPointSize = value;
                    updatePointSizePreferences();
                    IGV.getInstance().repaintDataPanels();
                }

            } catch (NumberFormatException numberFormatException) {
                JOptionPane.showMessageDialog(IGV.getMainFrame(),
                        "Point size must be an integer number.");
            }
        }
    }

    private void updatePointSizePreferences() {
        PreferenceManager.getInstance().put(PreferenceManager.GWAS_MIN_POINT_SIZE, String.valueOf(minPointSize));
        PreferenceManager.getInstance().put(PreferenceManager.GWAS_MAX_POINT_SIZE, String.valueOf(maxPointSize));


    }

    public boolean isLogNormalized() {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @SubtlyImportant
    private static GWASTrack getNextTrack(){
        return (GWASTrack) IGVSessionReader.getNextTrack();
    }
}
