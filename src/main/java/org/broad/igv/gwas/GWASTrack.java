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
import org.broad.igv.Globals;
import org.broad.igv.event.RepaintEvent;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.IGVPreferences;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.renderer.GraphicUtils;
import org.broad.igv.track.*;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.ui.panel.IGVPopupMenu;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.util.UIUtilities;
import org.broad.igv.util.ChromosomeColors;
import org.broad.igv.util.ResourceLocator;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.text.DecimalFormat;
import java.util.List;
import java.util.*;
import java.util.stream.Collectors;

/**
 * @author jussi
 * @since Nov 23, 2009
 */
public class GWASTrack extends AbstractTrack {

    private static final Logger log = Logger.getLogger(GWASTrack.class);
    private static final int AXIS_AREA_WIDTH = 60;
    private static final DecimalFormat formatter = new DecimalFormat();

    private int minPointSize;
    private int maxPointSize;
    private boolean useChrColors;
    private boolean singleColor;
    private boolean alternatingColors;
    private Color primaryColor;
    private Color secondaryColor;
    private double trackMinY;
    private double maxY;
    private double scale;
    private boolean drawYAxis = true;
    private boolean showAxis = true;
    double maxValue = -1;
    private Map<String, List<GWASFeature>> gData;
    Genome genome;
    private String[] columns;
    IGV igv;

    /**
     * Constructor for a new GWAS track
     *
     * @param locator
     * @param id
     * @param name
     * @param gData
     */
    public GWASTrack(ResourceLocator locator,
                     String id,
                     String name,
                     Map<String, List<GWASFeature>> gData,
                     String[] columns,
                     Genome genome) {
        super(locator, id, name);

        this.genome = genome;
        this.igv = IGV.getInstance(); // TODO replace with parameter

        IGVPreferences prefs = PreferencesManager.getPreferences();

        // Set range from 0 to highest value rounded to greater integer
        setMaxValue(gData);
        int mv = (int) Math.ceil(maxValue);
        super.setDataRange(new DataRange(0, (mv / 2), mv));

        // Get default values
        super.setHeight(prefs.getAsInt(Constants.GWAS_TRACK_HEIGHT));
        this.minPointSize = prefs.getAsInt(Constants.GWAS_MIN_POINT_SIZE);
        this.maxPointSize = prefs.getAsInt(Constants.GWAS_MAX_POINT_SIZE);
        this.primaryColor = ColorUtilities.stringToColor(prefs.get(Constants.GWAS_PRIMARY_COLOR));
        this.secondaryColor = ColorUtilities.stringToColor(prefs.get(Constants.GWAS_SECONDARY_COLOR));
        this.singleColor = prefs.getAsBoolean(Constants.GWAS_SINGLE_COLOR);
        this.alternatingColors = prefs.getAsBoolean(Constants.GWAS_ALTERNATING_COLORS);
        this.useChrColors = prefs.getAsBoolean(Constants.GWAS_USE_CHR_COLORS);
        this.showAxis = prefs.getAsBoolean(Constants.GWAS_SHOW_AXIS);
        this.gData = gData;
        this.columns = columns;
    }

    public GWASTrack() {
    }

    private void setMaxValue(Map<String, List<GWASFeature>> gData) {
        for (List<GWASFeature> features : gData.values()) {
            for (GWASFeature f : features) {
                if (f.value > maxValue) maxValue = f.value;
            }
        }
    }

    @Override
    public boolean isNumeric() {
        return true;
    }

    @Override
    public boolean isReadyToPaint(ReferenceFrame frame) {
        return true;  // Track is initialized with all data
    }

    @Override
    public void load(ReferenceFrame frame) {
        // Nothing to do, track is initialized with all data
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
        List<String> chrList;
        if (chrName.equals("All")) {
            chrList = genome.getLongChromosomeNames();
        } else {
            chrList = Arrays.asList(chrName);
        }
        double dx = Math.ceil(1 / locScale) + 1;
        double rangeMaxValue = Math.ceil(maxValue);

        double pointSizeScale = rangeMaxValue / maxPointSize;

        Color drawColor = this.primaryColor;

        int xMinPointSize = (int) (1 / locScale);

        // Loop through data points, chromosome by chromosome

        int chrCounter = 0;
        for (String chr : chrList) {

            List<GWASFeature> featureList = gData.get(chr);

            if (featureList != null) {

                // Choose a color for the chromosome
                // Use specific color for each chromosome
                if (this.useChrColors)
                    drawColor = ChromosomeColors.getColor(chr);

                    // Use two alternating colors for chromosomes
                else if (this.alternatingColors) {
                    if (chrCounter % 2 == 0)
                        drawColor = this.secondaryColor;
                    else
                        drawColor = this.primaryColor;
                }

                // Loop through data points in a chromosome
                Graphics2D g = context.getGraphics();
                for (GWASFeature feature : featureList) {

                    // Get location, e.g. start for the data point
                    int start;
                    if (chrName.equalsIgnoreCase("all"))
                        start = genome.getGenomeCoordinate(feature.chr, feature.position);
                    else
                        start = feature.position;

                    // Based on location, calculate X-coordinate, or break if outside of the view
                    double pX = ((start - origin) / locScale);
                    if ((pX + dx < 0))
                        continue;
                    else if (pX > adjustedRectMaxX)
                        break;

                    // Based on value of the data point, calculate Y-coordinate
                    double dataY = feature.value;

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

                        g.setColor(drawColor);
                        g.fillRect(x, y, maxDrawX - x, maxDrawY - y);
                        feature.pixelX = (x + maxDrawX) / 2;
                        feature.pixelY = (y + maxDrawY) / 2;
                    }
                }
            }
            chrCounter++;
        }

        // Draw the legend axis
        if (showAxis) {
            this.renderAxis(context, arect);
        }

    }


    void renderAxis(RenderContext context, Rectangle arect) {

        Rectangle drawingRect = calculateDrawingRect(arect);

        Color labelColor = PreferencesManager.getPreferences().getAsBoolean(Constants.CHART_COLOR_TRACK_NAME) ? getColor() : Color.black;
        Graphics2D labelGraphics = context.getGraphic2DForColor(labelColor);
        //Graphics2D labelGraphics = context.getGraphic2DForColor(Color.black);

        labelGraphics.setFont(FontManager.getFont(8));

        String tmpDisplayName = this.getDisplayName();
        if (tmpDisplayName != null && tmpDisplayName.length() > 0 && PreferencesManager.getPreferences().getAsBoolean(Constants.CHART_DRAW_TRACK_NAME)) {
            // Only attempt if track height is > 25 pixels
            if (arect.getHeight() > 25) {
                Rectangle labelRect = new Rectangle(arect.x, arect.y + 10, arect.width, 10);
                labelGraphics.setFont(FontManager.getFont(10));
                GraphicUtils.drawCenteredText(tmpDisplayName, labelRect, labelGraphics);
            }
        }

        //if (IGVPreferences.getInstance().getAsBoolean(IGVPreferences.CHART_DRAW_Y_AXIS)) {
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
     */
    GWASFeature findFeature(String chr, int mouseX, int mouseY) {
        final int pixelThreshold = 3;
        List<GWASFeature> featureList;
        if ("all".equalsIgnoreCase(chr)) {
            featureList = new ArrayList<>();
            for (List<GWASFeature> chrFeatures : this.gData.values()) {
                for (GWASFeature f : chrFeatures) {
                    if (Math.abs(mouseX - f.pixelX) <= pixelThreshold && Math.abs(mouseY - f.pixelY) <= pixelThreshold) {
                        featureList.add(f);
                    }
                }
            }
        } else {
            featureList = this.gData.get(chr);
        }
        if (featureList == null) {
            return null;
        } else {
            List<GWASFeature> closeFeatures = featureList.stream()
                    .filter(f -> Math.abs(mouseX - f.pixelX) <= pixelThreshold && Math.abs(mouseY - f.pixelY) <= pixelThreshold)
                    .collect(Collectors.toList());
            if (closeFeatures.isEmpty()) {
                return null;
            } else {
                closeFeatures.sort((o1, o2) -> {
                    double d1 = Math.sqrt((mouseX - o1.pixelX) * (mouseX - o1.pixelX) + (mouseY - o1.pixelY) * (mouseY - o1.pixelY));
                    double d2 = Math.sqrt((mouseX - o2.pixelX) * (mouseX - o2.pixelX) + (mouseY - o2.pixelY) * (mouseY - o2.pixelY));
                    if (d1 == d2) return 0;
                    else if (d1 > d2) return 1;
                    else return -1;
                });
                return closeFeatures.get(0);
            }
        }
    }

    /**
     * Get description for a data point at given chromosome and index. If description is not cached, re-populate cache.
     */

    public String getDescriptionString(GWASFeature feature) {

        String description = feature.line;
        String descriptionString = null;

        if (description != null) {
            descriptionString = "";
            int headersSize = this.columns.length;
            String[] tokens = Globals.singleTabMultiSpacePattern.split(description);

            for (int i = 0; i < headersSize; i++) {
                String tmpHeaderToken = this.columns[i];
                if (tmpHeaderToken != null)
                    descriptionString += tmpHeaderToken + ": " + tokens[i] + "<br>";
            }

        }
        return descriptionString;

    }


    public String getValueStringAt(String chr, double position, int mouseX, int mouseY, ReferenceFrame frame) {


        GWASFeature feature = findFeature(chr, mouseX, mouseY);

        // If there is a data point at the given location, fetch description
        return feature != null ? getDescriptionString(feature) : null;

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
        axisItem.addActionListener(e -> {
            showAxis = axisItem.isSelected();
            PreferencesManager.getPreferences().put(Constants.GWAS_SHOW_AXIS, String.valueOf(showAxis));
            igv.postEvent(new RepaintEvent(GWASTrack.this));
        });
        menu.add(axisItem);
        return axisItem;
    }

    public JMenuItem addChrColorItem(JPopupMenu menu) {
        JMenuItem colorItem = new JCheckBoxMenuItem("Chromosome color", useChrColors);
        colorItem.addActionListener(e -> {
            singleColor = false;
            useChrColors = true;
            alternatingColors = false;
            updateColorPreferences();
            igv.postEvent(new RepaintEvent(GWASTrack.this));
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
                igv.postEvent(new RepaintEvent(GWASTrack.this));
            }
        });
        menu.add(colorItem);

        return colorItem;
    }

    public JMenuItem addSingleColorItem(JPopupMenu menu) {
        JMenuItem colorItem = new JCheckBoxMenuItem("Single color", singleColor);
        colorItem.addActionListener(e -> {
            singleColor = true;
            useChrColors = false;
            alternatingColors = false;
            updateColorPreferences();
            igv.postEvent(new RepaintEvent(GWASTrack.this));
        });
        menu.add(colorItem);
        return colorItem;
    }

    private void updateColorPreferences() {
        PreferencesManager.getPreferences().put(Constants.GWAS_SINGLE_COLOR, String.valueOf(singleColor));
        PreferencesManager.getPreferences().put(Constants.GWAS_USE_CHR_COLORS, String.valueOf(useChrColors));
        PreferencesManager.getPreferences().put(Constants.GWAS_ALTERNATING_COLORS, String.valueOf(alternatingColors));
    }

    public JMenuItem addPrimaryColorItem(JPopupMenu menu) {
        JMenuItem colorItem = new JMenuItem("Set primary color...");
        colorItem.addActionListener(e -> {
            Color color = UIUtilities.showColorChooserDialog("Set primary color", primaryColor);
            if (color != null) {
                primaryColor = color;
                String colorString = ColorUtilities.colorToString(primaryColor);
                PreferencesManager.getPreferences().put(Constants.GWAS_PRIMARY_COLOR, colorString);
                igv.postEvent(new RepaintEvent(GWASTrack.this));
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
                    PreferencesManager.getPreferences().put(Constants.GWAS_SECONDARY_COLOR, colorString);
                    igv.postEvent(new RepaintEvent(GWASTrack.this));
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
                    igv.postEvent(new RepaintEvent(GWASTrack.this));
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
                    igv.postEvent(new RepaintEvent(GWASTrack.this));
                }

            } catch (NumberFormatException numberFormatException) {
                JOptionPane.showMessageDialog(IGV.getMainFrame(),
                        "Point size must be an integer number.");
            }
        }
    }

    private void updatePointSizePreferences() {
        PreferencesManager.getPreferences().put(Constants.GWAS_MIN_POINT_SIZE, String.valueOf(minPointSize));
        PreferencesManager.getPreferences().put(Constants.GWAS_MAX_POINT_SIZE, String.valueOf(maxPointSize));


    }

    public boolean isLogNormalized() {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public void marshalXML(Document document, Element element) {

        super.marshalXML(document, element);

        element.setAttribute("maxPointSize", String.valueOf(maxPointSize));
        element.setAttribute("minPointSize", String.valueOf(minPointSize));
        element.setAttribute("scale", String.valueOf(scale));
        element.setAttribute("trackMinY", String.valueOf(trackMinY));
        element.setAttribute("maxY", String.valueOf(maxY));
        element.setAttribute("useChrColors", String.valueOf(useChrColors));
        element.setAttribute("singleColor", String.valueOf(singleColor));
        element.setAttribute("drawYAxis", String.valueOf(drawYAxis));
        element.setAttribute("alternatingColors", String.valueOf(alternatingColors));
        if (primaryColor != null) {
            element.setAttribute("primaryColor", ColorUtilities.colorToString(primaryColor));
        }
        if (secondaryColor != null) {
            element.setAttribute("secondaryColor", ColorUtilities.colorToString(secondaryColor));
        }


    }

    @Override
    public void unmarshalXML(Element element, Integer version) {

        super.unmarshalXML(element, version);

        maxPointSize = Integer.parseInt(element.getAttribute("maxPointSize"));
        minPointSize = Integer.parseInt(element.getAttribute("minPointSize"));
        scale = Double.parseDouble(element.getAttribute("scale"));
        trackMinY = Double.parseDouble(element.getAttribute("trackMinY"));
        maxY = Double.parseDouble(element.getAttribute("maxY"));
        useChrColors = Boolean.parseBoolean(element.getAttribute("useChrColors"));
        singleColor = Boolean.parseBoolean(element.getAttribute("singleColor"));
        drawYAxis = Boolean.parseBoolean(element.getAttribute("drawYAxis"));
        alternatingColors = Boolean.parseBoolean(element.getAttribute("alternatingColors"));
        if (element.hasAttribute("primaryColor")) {
            primaryColor = ColorUtilities.stringToColor(element.getAttribute("primaryColor"));
        }
        if (element.hasAttribute("secondaryColor ")) {
            secondaryColor = ColorUtilities.stringToColor(element.getAttribute("secondaryColor"));
        }
    }
}
