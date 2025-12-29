package org.broad.igv.gwas;

import org.broad.igv.logging.*;
import org.broad.igv.Globals;
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
import java.util.regex.Pattern;
import java.util.stream.Collectors;

/**
 * @author jussi
 * @since Nov 23, 2009
 */
public class GWASTrack extends AbstractTrack {

    private static final Logger log = LogManager.getLogger(GWASTrack.class);
    private static final int AXIS_AREA_WIDTH = 60;
    private static final DecimalFormat formatter = new DecimalFormat();

    private int minPointSize;
    private int maxPointSize;
    private boolean useChrColors;
    private boolean singleColor;
    private boolean alternatingColors;
    private Color primaryColor;
    private Color secondaryColor;
    private double maxY;
    private boolean drawYAxis = true;
    private boolean showAxis = true;
    double maxValue = -1;
    private GWASData gData;
    Genome genome;
    private String[] columns;

    private Pattern delimiter;
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
                     GWASData gData,
                     String[] columns,
                     Pattern delimiter,
                     Genome genome) {
        super(locator, id, name);

        this.delimiter = delimiter;
        this.genome = genome;
        this.igv = IGV.getInstance(); // TODO replace with parameter

        IGVPreferences prefs = PreferencesManager.getPreferences();

        // Set inital range from 0 to highest value rounded to greater integer, or 25.  This assumes pvalues
        float maxValue = Math.min(25, (float) Math.ceil(gData.getMaxValue()));
        float minValue = (float) Math.floor(gData.getMinValue());
        super.setDataRange(new DataRange(minValue, maxValue));

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

    public void render(RenderContext context, Rectangle rect) {

        Genome genome = GenomeManager.getInstance().getCurrentGenome();

        double origin = context.getOrigin();
        double locScale = context.getScale();

        Graphics2D g = null;

        try {
            g = (Graphics2D) context.getGraphics().create();
            g.setClip(rect);

            Rectangle plotRect = calculateDrawingRect(rect);

            DataRange dataRange = this.getDataRange();
            float maxValue = dataRange.getMaximum();
            float minValue = dataRange.getMinimum();
            double pointSizeScaleFactor = ((double) (this.maxPointSize - this.minPointSize)) / (maxValue - minValue);

            String chrName = context.getChr();
            List<String> chrList;
            if (chrName.equals("All")) {
                chrList = genome.getLongChromosomeNames();
            } else {
                chrList = Arrays.asList(chrName);
            }

            Color drawColor = this.primaryColor;


            // Loop through data points, chromosome by chromosome

            int chrCounter = 0;
            for (String chr : chrList) {

                List<GWASFeature> featureList = gData.get(chr);

                if (featureList != null) {

                    if (this.useChrColors)
                        drawColor = ChromosomeColors.getColor(chr);
                    else if (this.alternatingColors) {
                        if (chrCounter % 2 == 0)
                            drawColor = this.secondaryColor;
                        else
                            drawColor = this.primaryColor;
                    }

                    // Loop through data points
                    for (GWASFeature feature : featureList) {

                        // Get location, e.g. start for the data point
                        int start;
                        if (chrName.equalsIgnoreCase("all"))
                            start = genome.getGenomeCoordinate(feature.chr, feature.position) - 1;
                        else
                            start = feature.position - 1;


                        // Based on value of the data point, calculate Y-coordinate
                        double dataY = feature.value;

                        if (!Double.isNaN(dataY)) {

                            int pixelSize = Math.min(this.maxPointSize, this.minPointSize + (int) (pointSizeScaleFactor * (dataY - minValue)));
                            int x = (int) ((start - origin + 0.5) / locScale) - pixelSize / 2;
                            int y = computeYPixelValue(plotRect, this.dataRange, dataY) - pixelSize / 2;
                            g.setColor(drawColor);
                            g.fillOval(x, y, pixelSize, pixelSize);
                            feature.pixelX = x;
                            feature.pixelY = y;
                            feature.pixelSize = pixelSize;
                        }
                    }
                }
                chrCounter++;
            }

            // Draw the legend axis
            if (showAxis) {
                this.renderAxis(context, plotRect);
            }

        } finally {
            if (g != null) g.dispose();
        }

    }


    void renderAxis(RenderContext context, Rectangle plotRect) {

        Graphics2D g = null;

        try {
            g = (Graphics2D) context.getGraphics().create();

            g.setColor(Color.black);
            g.setFont(FontManager.getFont(11));

            Rectangle axisRect = new Rectangle(plotRect.x, plotRect.y + 1, AXIS_AREA_WIDTH, plotRect.height);
            DataRange axisDefinition = getDataRange();
            float maxValue = axisDefinition.getMaximum();
            float minValue = axisDefinition.getMinimum();

            // Bottom (minimum tick mark)
            int pY = computeYPixelValue(plotRect, axisDefinition, minValue);

            g.drawLine(axisRect.x + AXIS_AREA_WIDTH - 10, pY, axisRect.x + AXIS_AREA_WIDTH - 5, pY);
            GraphicUtils.drawRightJustifiedText(formatter.format(minValue), axisRect.x + AXIS_AREA_WIDTH - 15, pY, g);

            // Top (maximum tick mark)
            int topPY = computeYPixelValue(plotRect, axisDefinition, maxValue);

            g.drawLine(axisRect.x + AXIS_AREA_WIDTH - 10, topPY, axisRect.x + AXIS_AREA_WIDTH - 5, topPY);
            GraphicUtils.drawRightJustifiedText(formatter.format(maxValue), axisRect.x + AXIS_AREA_WIDTH - 15, topPY + 4, g);

            // Connect top and bottom
            g.drawLine(axisRect.x + AXIS_AREA_WIDTH - 10, topPY, axisRect.x + AXIS_AREA_WIDTH - 10, pY);

            // Draw middle tick marks from top to bottom, draw only if space available
            int lastY = pY;
            for (int i = (int) minValue + 1; i < (int) maxValue; i++) {

                int midPY = computeYPixelValue(plotRect, axisDefinition, i);
                if ((midPY < lastY - 15) && (midPY < pY - 15) && (midPY > topPY + 15)) {
                    g.drawLine(axisRect.x + AXIS_AREA_WIDTH - 10, midPY, axisRect.x + AXIS_AREA_WIDTH - 5, midPY);
                    GraphicUtils.drawRightJustifiedText(formatter.format(i), axisRect.x + AXIS_AREA_WIDTH - 15, midPY + 4, g);
                    lastY = midPY;
                }
            }
        } finally {
            if (g != null) g.dispose();
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

        List<GWASFeature> closeFeatures;
        Collection<String> chromosomes = "all".equalsIgnoreCase(chr) ? this.gData.keySet() : Arrays.asList(chr);
        closeFeatures = new ArrayList<>();
        for (String c : chromosomes) {
            List<GWASFeature> chrFeatures = this.gData.get(c);
            for (GWASFeature f : chrFeatures) {
                int pixelSize = Math.max(3, f.pixelSize);
                if (Math.abs(mouseX - f.pixelX) <= pixelSize && Math.abs(mouseY - f.pixelY) <= pixelSize) {
                    closeFeatures.add(f);
                }
            }
        }
        if (closeFeatures.isEmpty()) {
            return null;
        } else if (closeFeatures.size() == 1) {
            return closeFeatures.get(0);
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


    /**
     * Get description for a data point at given chromosome and index. If description is not cached, re-populate cache.
     */

    public String getDescriptionString(GWASFeature feature) {

        String description = feature.line;
        String descriptionString = null;

        if (description != null) {
            descriptionString = "";
            int headersSize = this.columns.length;
            String[] tokens = this.delimiter.split(description);

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
            repaint();
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
            repaint();
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
                repaint();
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
            repaint();
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
                repaint();
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
                    repaint();
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
                IGV.getInstance().getMainFrame(), "Minimum point size in pixels (1-20):", String.valueOf(value));
        if (!(tmpValue == null) || !tmpValue.trim().equals("")) {

            try {
                value = Integer.parseInt(tmpValue);
                if (value < 1 || value > 20) {
                    JOptionPane.showMessageDialog(IGV.getInstance().getMainFrame(),
                            "Minimum point size must be an integer number between 1 and 20.");
                } else {
                    if (value > this.maxPointSize) {
                        this.maxPointSize = value;
                    }
                    this.minPointSize = value;
                    updatePointSizePreferences();
                    repaint();
                }
            } catch (NumberFormatException numberFormatException) {
                JOptionPane.showMessageDialog(IGV.getInstance().getMainFrame(),
                        "Point size must be an integer number.");
            }
        }
    }

    private void changeMaxPointSizeValue() {

        int value = this.maxPointSize;
        String tmpValue = JOptionPane.showInputDialog(
                IGV.getInstance().getMainFrame(), "Maximum point size in pixels (1-20):",
                String.valueOf(value));

        if (!(tmpValue == null) || !tmpValue.trim().equals("")) {

            try {
                value = Integer.parseInt(tmpValue);
                if (value < 1 || value > 20) {
                    JOptionPane.showMessageDialog(IGV.getInstance().getMainFrame(),
                            "Maximum point size must be an integer number between 1 and 20.");
                } else {
                    if (value < this.minPointSize) {
                        this.minPointSize = value;
                    }
                    this.maxPointSize = value;
                    updatePointSizePreferences();
                    repaint();
                }

            } catch (NumberFormatException numberFormatException) {
                JOptionPane.showMessageDialog(IGV.getInstance().getMainFrame(),
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
        element.setAttribute("useChrColors", String.valueOf(useChrColors));
        element.setAttribute("singleColor", String.valueOf(singleColor));
        element.setAttribute("showAxis", String.valueOf(showAxis));
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
        useChrColors = Boolean.parseBoolean(element.getAttribute("useChrColors"));
        singleColor = Boolean.parseBoolean(element.getAttribute("singleColor"));
        showAxis = Boolean.parseBoolean(element.getAttribute("showAxis")) || Boolean.parseBoolean(element.getAttribute("drawYAxis"));
        alternatingColors = Boolean.parseBoolean(element.getAttribute("alternatingColors"));
        if (element.hasAttribute("primaryColor")) {
            primaryColor = ColorUtilities.stringToColor(element.getAttribute("primaryColor"));
        }
        if (element.hasAttribute("secondaryColor ")) {
            secondaryColor = ColorUtilities.stringToColor(element.getAttribute("secondaryColor"));
        }
    }
}
