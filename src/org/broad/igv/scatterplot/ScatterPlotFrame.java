package org.broad.igv.scatterplot;

import org.jfree.chart.*;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.block.BlockBorder;
import org.jfree.chart.block.BlockContainer;
import org.jfree.chart.block.BorderArrangement;
import org.jfree.chart.block.LabelBlock;
import org.jfree.chart.labels.StandardXYToolTipGenerator;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.chart.title.LegendTitle;
import org.jfree.chart.title.TextTitle;
import org.jfree.chart.title.Title;
import org.jfree.data.xy.XYDataItem;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.HorizontalAlignment;
import org.jfree.ui.RectangleEdge;

import javax.swing.*;
import javax.swing.border.EtchedBorder;
import javax.swing.border.TitledBorder;
import java.awt.*;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: martind
 * Date: Feb 3, 2010
 * Time: 4:45:45 PM
 * To change this template use File | Settings | File Templates.
 */
/*
*   Scatterplot implementation class - provides JFrame display of IGV data in a
*   scatterplot panel, with border panels for selection of plotting parameters.
*
*   Scatterplot of IGV data appears with default settings.
*
*   Selection comboboxes in the NORTH border panel provide selection of
*   data types for the X and Y axis and an attribute to define series plots.
*
*  Selection combobox(s) in the EAST panel provide attribute selection
*  for series plots of categories for selected attributes.
*
*  Note: Currently WEST border panel is not used.
*  Scrollable sample names are displayed in the WEST border panel, and can be
*  used to indicate selected data points for future operations.
*
*  Note: Currently SOUTH border panel is not used.
*  Status can be displayed in the SOUTH border panel, once enabled.
*
* */
public class ScatterPlotFrame extends JFrame {

    // **************** IGV ScatterPlotFrame data and identifier **************
    private String igvTitle;                 // title for scatterplot
    private ScatterPlotData igvPlotData;   // imported IGV data container class
    private JFreeChart igvScatterPlotChart;   // scatter plot created from selected igvPlotData

    // ***************** IGV ScatterPlotFrame Frame **************************
    private int igvXFrameStart;             // X start position for scatterplot frame
    private int igvYFrameStart;             // Y start position for scatterplot frame
    private int igvXFrameDimension;         // X dimension for scatterplot frame
    private int igvYFrameDimension;         // Y dimension for scatterplot frame
    private Container igvFrameContent;      // content associated with scatterplot JFrame
    private BorderLayout igvBorderLayout;   // borderlayout object
    private JPanel igvSampleListPanel;      // scrollable sample names panel - WEST
    private ChartPanel igvPlotPanel;        // scatterplot data panel - CENTER
    private JPanel igvAttributeLegendPanel; // JPanel containing symbol selection - EAST

    //**************** Scatterplot Chart settings *********************
    private JFreeChart igvScatterChart;     // scatterplot as a JFreeChart
    //private JPanel igvUpdatePanel;          // JPanel containing "Update" button - embedded in EAST border
    private JButton igvUpdateButton;        // enables update of scatterplot for new picks (currently disabled)

    private JComboBox igvXAxisComboBox;     // X Axis combobox - shows axis labels
    private int igvXAxisIndex;              // X Axis index for chosen label
    private String igvXAxisName;            // X Axis name selection

    private JComboBox igvYAxisComboBox;     // Y Axis combobox - shows axis labels
    private int igvYAxisIndex;              // Y Axis index for chosen label
    private String igvYAxisName;            // Y axis name selection

    private JComboBox igvSymbolComboBox;    // symbol attributes combobox - shows sample attributes
    private int igvSymbolIndex;             // symbol index for chosen attribute
    private String igvAttributeName;        // attribute name selection

    // *************** scatterplot data and datasets ********************
    // Note: every time there is a selction from igvSymbolComboBox, new shape panels
    // and attribute filtering legend needs to be rebuilt
    private String[] igvDataNames;       // names for axis assigned data types
    private String[] igvSymbolNames;     // names for plot series data attributes

    // ****************** series data, plot shapes and colors *****************
    private XYSeriesCollection igvXYSeriesDataset; // collection of series for plot selections

    // SymbolSettings: category, shape, color, isVisible
    private SymbolSettings[] igvSymbolSettings;
    private SymbolFactory igvSymbolFactory;      // probably don't need to save this

    /**
     * ScatterPlotFrame - ScatterPlotFrame constructor provides plotting choices for IGV
     * scatterplot data display panel.
     * <p/>
     * Parameters:
     * title - String containing JFreeChart title to be displayed
     * scatterPlotData -  ScatterPlotData container for IGV data to be plotted
     */
    public ScatterPlotFrame(String title, ScatterPlotData scatterPlotData) {
        // construct super with version title
        super("IGV Scatter Plot Version 1.0");

        // accept the igv scatterplot data
        igvTitle = title;
        igvPlotData = scatterPlotData;

        // assign default axis and selections
        igvXAxisIndex = 0;      // default choice
        igvYAxisIndex = 1;      // deafult choice
        igvXAxisName = igvPlotData.getDataKeyName(igvXAxisIndex);
        igvYAxisName = igvPlotData.getDataKeyName(igvYAxisIndex);

        // get data names for axis selctions
        ArrayList<String> dataNameList = igvPlotData.getDataNames();
        igvDataNames = dataNameList.toArray(new String[0]);

        // assign symbol selection and symbol filters
        ArrayList<String> symbolNameList = igvPlotData.getSymbolNames();

        igvSymbolNames = symbolNameList.toArray(new String[0]);
        igvSymbolIndex = 0;
        igvAttributeName = igvPlotData.getSymbolKeyName(igvSymbolIndex);

        // create the application frame for IGV ScatterPlotFrame
        Toolkit theKit = this.getToolkit();       // Get the window toolkit
        Dimension wndSize = theKit.getScreenSize();  // Get screen size

        // Set the JFrame position to screen center & size to half screen size
        igvXFrameStart = wndSize.width / 4;
        igvYFrameStart = wndSize.height / 4;
        igvXFrameDimension = wndSize.width / 2;
        igvYFrameDimension = wndSize.height / 2;
        super.setBounds(igvXFrameStart, igvYFrameStart,   // Position
                igvXFrameDimension, igvYFrameDimension);  // Size
        super.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        igvBorderLayout = new BorderLayout();           // Create a layout manager
        igvFrameContent = super.getContentPane();        // Get the igvFrameContent pane
        igvFrameContent.setLayout(igvBorderLayout);      // Set the container layout mgr
        EtchedBorder edge = new EtchedBorder(EtchedBorder.RAISED);  // Button border

        // create axis selection panels - drop in NORTH border panel
        JPanel axisSelectPanel = createDataSelectionPanel(igvDataNames,
                igvXAxisIndex, igvYAxisIndex, igvSymbolIndex, igvSymbolNames);
        igvFrameContent.add(axisSelectPanel, BorderLayout.NORTH);

        // Symbol selection filter panel - drop in EAST panel
        igvAttributeLegendPanel = createAttributeSymbolPanel(igvAttributeName);
        igvFrameContent.add(igvAttributeLegendPanel, BorderLayout.EAST);

        // create JFreeChart scatterplot panel - drop in CENTER panel
        igvPlotPanel = createScatterPlotPanel(igvTitle, igvPlotData, igvXAxisName, igvYAxisName,
                igvAttributeName, igvSymbolSettings);
        igvFrameContent.add(igvPlotPanel, BorderLayout.CENTER);

        // List of sample names - drop in WEST panel
        //igvSampleListPanel = createSampleListPanel(igvPlotData.getSampleNames());
        //igvFrameContent.add(igvSampleListPanel, BorderLayout.WEST);

        // Status and point selection info display - drop in SOUTH panel
        //igvFrameContent.add(igvSampleInfo = new JLabel("Status"), BorderLayout.SOUTH);

    }

    /*
    *   Creates the selection panels for the choosing scatter plot data,
    *   where comboboxes allow selection of the x axis and y axis data measurement
    *   variables and a symbol selection combobox provides a plotting series
    *   for data attribute selection.
    *
    *   Parameters:
    *       dataNames - string array containing data measurement names
    *       xAxisIndex - index to initialize x axis data measurement selection
    *       yAxisIndex - index to initialize y axis data measurement selection
    *       symbolIndex - index to initialize data attribute symbol selection
    *       symbolNames - names of attributes for symbol selection
    *
    *   Returns:
    *       Axis selection panel
    *
    *   Note: These index selections determine the initial plot variables
    *   for the XYPlot of the scatterplot.
    * */
    private JPanel createDataSelectionPanel(String[] dataNames,
                                            int xAxisIndex, int yAxisIndex, int symbolIndex, String[] symbolNames) {

        JPanel axisSelectPanel = new JPanel();        // axis selection panel
        axisSelectPanel.setLayout(new BoxLayout(axisSelectPanel, BoxLayout.X_AXIS)); // left to right
        axisSelectPanel.add(Box.createHorizontalGlue());

        // data meaurement types selected for axis
        JPanel xAxisPanel = createXAxisSelectionPanel(dataNames, xAxisIndex);
        axisSelectPanel.add(xAxisPanel);
        axisSelectPanel.add(Box.createHorizontalStrut(20));
        JPanel yAxisPanel = createYAxisSelectionPanel(dataNames, yAxisIndex);
        axisSelectPanel.add(yAxisPanel);

        // attribute symbol type selections for series plots
        axisSelectPanel.add(Box.createHorizontalStrut(20));
        JPanel symbolSelectionPanel = createAttributeSelectionPanel(symbolNames, igvSymbolIndex);
        axisSelectPanel.add(symbolSelectionPanel);
        axisSelectPanel.add(Box.createHorizontalGlue());

        return axisSelectPanel;
    }

    /*
    *   Method create X Axis data measurement selection panel
    *
    *   Parameters:
    *       dataNames - string array of selectable data measurement labels
    *       selectedIndex - index to display as selected
    *
    *   Returns:
    *       X Axis selection panel
    * */
    private JPanel createXAxisSelectionPanel(String[] dataNames, int selectedIndex) {
        JPanel xAxisPanel = new JPanel();

        // layout top to bottom
        xAxisPanel.setLayout(new BoxLayout(xAxisPanel, BoxLayout.Y_AXIS));

        xAxisPanel.add(new JLabel("X Axis"));
        igvXAxisComboBox = new JComboBox(dataNames);
        igvXAxisComboBox.setMaximumRowCount(dataNames.length);
        igvXAxisComboBox.setMaximumSize(new Dimension(150, 20));
        igvXAxisComboBox.setMinimumSize(new Dimension(150, 20));
        igvXAxisComboBox.setSelectedIndex(selectedIndex);
        xAxisPanel.add(igvXAxisComboBox);           // add combobox to JFrame

        igvXAxisComboBox.addItemListener(
                new ItemListener() // anonymous inner class
                {
                    // handle JComboBox event
                    public void itemStateChanged(ItemEvent event) {
                        // determine whether checkbox selected
                        if (event.getStateChange() == ItemEvent.SELECTED)
                            igvXAxisIndex = igvXAxisComboBox.getSelectedIndex();

                        // recreate the scatterplot
                        drawScatterplot();
                    } // end method itemStateChanged
                } // end anonymous inner class
        ); // end call to addItemListener.addItemListener(

        return xAxisPanel;
    }

    /*
    *   Method create Y Axis data measuremnt selection panel.
    *
    *   Parameters:
    *       dataNames - string array of selectable data measurement labels
    *       selectedIndex - index to display as selected
    *
     *   Returns:
    *       Y Axis selection panel
    * */
    private JPanel createYAxisSelectionPanel(String[] dataNames, int selectedIndex) {
        JPanel yAxisPanel = new JPanel();

        // layout top to bottom
        yAxisPanel.setLayout(new BoxLayout(yAxisPanel, BoxLayout.Y_AXIS));

        yAxisPanel.add(new JLabel("Y Axis"));
        this.igvYAxisComboBox = new JComboBox(dataNames);
        this.igvYAxisComboBox.setMaximumRowCount(dataNames.length);
        this.igvYAxisComboBox.setMaximumSize(new Dimension(150, 20));
        this.igvYAxisComboBox.setMinimumSize(new Dimension(150, 20));
        this.igvYAxisComboBox.setSelectedIndex(selectedIndex);
        yAxisPanel.add(this.igvYAxisComboBox);           // add combobox to JFrame

        this.igvYAxisComboBox.addItemListener(
                new ItemListener() // anonymous inner class
                {
                    // handle JComboBox event
                    public void itemStateChanged(ItemEvent event) {
                        // determine whether checkbox selected
                        if (event.getStateChange() == ItemEvent.SELECTED)
                            igvYAxisIndex = ScatterPlotFrame.this.igvYAxisComboBox.getSelectedIndex();

                        drawScatterplot();
                    } // end method itemStateChanged
                } // end anonymous inner class
        ); // end call to addItemListener.addItemListener

        return yAxisPanel;
    }

    /*
    *   Method creates attribute name selection panel.
    *
    *   Parameters:
    *       symbolNames - array of selectable attribute names
    *       selectedIndex - index to display as the selected attribute
    *
    *   Returns:
    *       Attribute selection panel
    * */
    private JPanel createAttributeSelectionPanel(String[] symbolNames, int selectedIndex) {
        JPanel symbolPanel = new JPanel();

        // layout top to bottom
        symbolPanel.setLayout(new BoxLayout(symbolPanel, BoxLayout.Y_AXIS));

        symbolPanel.add(new JLabel("Attributes"));
        igvSymbolComboBox = new JComboBox(symbolNames);
        igvSymbolComboBox.setMaximumRowCount(symbolNames.length);
        igvSymbolComboBox.setMaximumSize(new Dimension(150, 20));
        igvSymbolComboBox.setMinimumSize(new Dimension(150, 20));
        igvSymbolComboBox.setSelectedIndex(selectedIndex);
        symbolPanel.add(igvSymbolComboBox);      // add combobox to JPanel

        igvSymbolComboBox.addItemListener(
                new ItemListener() // anonymous inner class
                {
                    // handle JComboBox event
                    public void itemStateChanged(ItemEvent event) {
                        // determine whether checkbox selected
                        if (event.getStateChange() == ItemEvent.SELECTED)
                            igvSymbolIndex = ScatterPlotFrame.this.igvSymbolComboBox.getSelectedIndex();

                        // recreate the scatterplot
                        drawScatterplot();
                    } // end method itemStateChanged
                } // end anonymous inner class
        ); // end call to addItemListener.addItemListener(

        return symbolPanel;
    }

    /*
    *   Method creates a panel for display of attribute category for filter selection.
    *
    *   Parameters:
    *       attribute - selected data attribute; categories extracted for symbol selection panes
    *
    * */
    private JPanel createAttributeSymbolPanel(String attribute) {

        // attribute lgend panel with vertical layout
        JPanel attributeLegendPanel = new JPanel();        // symbol selection panel
        TitledBorder attributeBorder = new TitledBorder(attribute + ":");
        attributeLegendPanel.setBorder(attributeBorder);
        attributeLegendPanel.setLayout(new BoxLayout(attributeLegendPanel, BoxLayout.Y_AXIS));

        // create the plot series symbol settings
        // First get the categories for selected attribute
        String[] categories = igvPlotData.getAttributeCategories(attribute);

        // set up symbol shapes for series plotting
        int shapeWidth = 4;
        int shapeHeight = 3;
        boolean isOutlined = true;
        boolean isFilled = true;

        // Note: Example use for three different SymbolFactory constructors is shown.
        // Null constructor requires caller to use SymbolFactory.addSeriesSymbol
        // method to define symbol settings.

        //  igvSymbolFactory = new SymbolFactory();
        // igvSymbolFactory = new SymbolFactory(categories, SymbolFactory.IGVSymbolShape.Rectangle,
        //     shapeWidth, shapeHeight, isOutlined, isFilled);
        igvSymbolFactory = new SymbolFactory(categories, shapeWidth, shapeHeight, isOutlined, isFilled);


        igvSymbolSettings = igvSymbolFactory.getSeriesSymbolSettings();

        // create symbol panels; one panel per attribute category
        int nSeries = igvSymbolSettings.length;
        for (int series = 0; series < nSeries; ++series) {

            SymbolPanel symbolPanel = new SymbolPanel(igvSymbolSettings[series]);
            attributeLegendPanel.add(symbolPanel);
        }

        return attributeLegendPanel;
    }

    /*
   *   Inner class creates symbol selection panel for attribute categories
   * */
    private class SymbolPanel extends JPanel {

        private SymbolSettings settings;
        private JCheckBox categoryCheckBox; // selected = show; unselected = hide

        public SymbolPanel(SymbolSettings settings) {

            //store settings
            this.settings = settings;

            // horizontal layout
            setLayout(new FlowLayout(FlowLayout.LEFT));
            this.setMinimumSize(new Dimension(150, 50));
            this.setMaximumSize(new Dimension(150, 50));
            this.setPreferredSize(new Dimension(150, 50));

            // Note: panel style    color icon    name    checkbox
            // replaced by          checkbox      name
            // because the need to control color assigned to series icon wouldn't work

            //JPanel icon = new LegendIcon(settings.getPlotShape(), settings.getPlotColor());
            //this.add(icon);

            //JLabel name = new JLabel(settings.getCategory());
            //this.add(name);

            categoryCheckBox = new JCheckBox(settings.getCategory());
            this.add(categoryCheckBox);
            categoryCheckBox.setSelected(settings.isVisible());
            CheckBoxHandler handler = new CheckBoxHandler();
            this.categoryCheckBox.addItemListener(handler);
        }

        // private inner class for ItemListener event handling
        private class CheckBoxHandler implements ItemListener {

            // respond to checkbox events
            public void itemStateChanged(ItemEvent event) {

                // enable/disable display of symbol panel category
                // through the global symbol settings
                boolean visible = categoryCheckBox.isSelected();
                igvSymbolSettings[settings.getSeries()].setVisible(visible);

                updateScatterplot();

            } //end of itemStateChanged
        } // end of CheckBoxHandler

        // inner class for SymbolPanel icon generation
        private class LegendIcon extends JPanel {
            private Shape shape;
            private Color color;

            public LegendIcon(Shape shape, Color color) {
                // note: the size needs to be scaled up; panel resolution different than plot?
                Rectangle2D.Float iconShape = new Rectangle2D.Float();
                Rectangle bounds = shape.getBounds();
                iconShape.setRect(0, 0, 3 * bounds.width, 3 * bounds.height);
                this.shape = iconShape;
                this.color = color;
            }

            // draw shapes with Java 2D API
            public void paint(Graphics g) {
                Graphics2D g2 = (Graphics2D) g;
                g2.setColor(color);
                g2.fill(shape);
            }

            public int getIconWidth() {
                Rectangle bounds = shape.getBounds();
                return bounds.width;
            }

            public int getIconHeight() {
                Rectangle bounds = shape.getBounds();
                return bounds.height;
            }

            public Color getColor() {
                return color;
            }
        }

    }  // end of class SymbolPanel

    /*
    *   createToolTips - method generates tooltips for the imported scatterplot
    *   data set igvPlotData.
    * */
    private void createToolTips(XYItemRenderer renderer) {

        // generate tooltips interface
        //StandardXYToolTipGenerator tooltips =

        renderer.setBaseToolTipGenerator(new StandardXYToolTipGenerator() {
            // callback fo ttoltip display returns a String
            public String generateToolTip(XYDataset data, int series, int item) {
                XYSeries xySeries = ((XYSeriesCollection) data).getSeries(series);

                IGVXYDataItem dataItem = (IGVXYDataItem) xySeries.getDataItem(item);
                String htmlResult = dataItem.getLabel();

                return htmlResult;
            }
        });
    }

    /*
    *   IGVXYDataItem - subclass provides storage of sample description
    *   along with the z and y data values of sample points.
    * */
    private class IGVXYDataItem extends XYDataItem {
        private String label = null;

        public IGVXYDataItem(Number number, Number number1, String label) {
            super(number, number1);

            this.label = label;
        }

        public IGVXYDataItem(double n, double n1, String label) {
            super(n, n1);

            this.label = label;
        }

        public String getLabel() {
            return label;
        }
    }

    /**
     * createScatterPlotPanel - method reates a scatterplot chart panel subclassed
     * from JPanel to allow setChart replotting of selected scatterplot data attributes.
     * <p/>
     * Parameters:
     * title - string contains title of plot
     * igvData - ScatterPlotFrame container class for IGV data samples
     * xAxisName - IGV data measurement name for x Axis coordinate;
     * must be a key defined in igvData dataMap
     * yAxisName - IGV data measurement name for y Axis coordinate;
     * must be a key defined in igvData igvDataMap
     * attribute - sample attribute selected for series symbol display
     * symbolSettings - Symbol display settings include shape, color, fill, visibility
     * <p/>
     * Returns:
     * scatterplot chart panel.
     */
    private ChartPanel createScatterPlotPanel(String title, ScatterPlotData igvData,
                                              String xAxisName, String yAxisName, String attribute, SymbolSettings[] symbolSettings) {

        // create the scatterplot as an XYDataSet plot
        igvScatterChart = createChart(title, igvData, xAxisName, yAxisName, attribute, symbolSettings);
        ChartPanel panel = new ChartPanel(igvScatterChart);
        panel.setMouseWheelEnabled(true);
        return panel;
    }

    /*
    *   createChart - method creates a JFree ChartPanel (derived from JPanel)
    *   which plots the XYDataset assigned form ScatterPlotData selections.
    *
    *   Parameters:
    *       title - string contains title of plot
    *       igvData - ScatterPlotFrame container class for IGV data samples
    *       xAxisName - IGV data measurement name for x Axis coordinate;
    *           must be a key defined in igvData dataMap.
    *       yAxisName - IGV data measurement name for y Axis coordinate;
    *           must be a key defined in igvData dataMap.
    *       attribute - sample attribute selected for series symbol display
    *       symbolSettings - symbol settings for series symbol display.
    *
    *   Returns:
    *       Scatterplot chart
    *
    * */
    private JFreeChart createChart(String title, ScatterPlotData igvData,
                                   String xAxisName, String yAxisName, String attribute, SymbolSettings[] symbolSettings) {

        // Selection of X and Y axis choices define data to be displayed
        double[] xValues = igvData.getDataValues(xAxisName);
        double[] yValues = igvData.getDataValues(yAxisName);

        // check for valid data and axis selection - and if error, return null
        if (yValues == null | xValues == null | yValues.length != xValues.length)
            return null;

        // create a composite x,y data set for all sample points
        int allPoints = xValues.length;
        double[][] xyData = new double[2][allPoints];

        // JfreeChart requires a two dimensional array of plot points
        for (int i = 0; i < allPoints; ++i) {
            xyData[0][i] = xValues[i];
            xyData[1][i] = yValues[i];
        }

        // Note: Currently only one attribute can be selected for scatter plot,
        // therefore using the 1st series (key = index 0) to identify the attribute
        int nSeries = symbolSettings.length;

        // extract sample categories for the selected attribute name
        String[] attributeValues = igvData.getSymbolValues(attribute);

        // create series collection to hold xy series datasets for JFreeChart
        igvXYSeriesDataset = new XYSeriesCollection();

        for (int series = 0; series < nSeries; ++series) {

            // get series category
            Comparable category = symbolSettings[series].getCategory();

            // This should be possible, but
            if (category == null) continue;

            // allocate series x,y dataset
            XYSeries xySeries = new XYSeries(category);

            // create plot series
            for (int dataIndex = 0; dataIndex < allPoints; ++dataIndex) {
                // if attribute value is same as category - assign data point to the series
                final String attributeValue = attributeValues[dataIndex];
                if (category.equals(attributeValue)) {
                    // get tooltips and assign data point to series
                    boolean isHTML = true;
                    String tooltips = igvPlotData.getSampleDescription(dataIndex, isHTML);
                    IGVXYDataItem item = new IGVXYDataItem(xyData[0][dataIndex], xyData[1][dataIndex], tooltips);
                    xySeries.add(item);
                }
            }

            // add series  dataset to series collection
            igvXYSeriesDataset.addSeries(xySeries);
        }

        // load specifications for scatterplot chart production
        igvScatterPlotChart = ChartFactory.createScatterPlot(null, xAxisName, yAxisName,
                igvXYSeriesDataset, PlotOrientation.VERTICAL, false, true, false);

        // Note: this should work according to Developer's Guide - but doesn't
        /*
       TextTitle chartTitle = igvScatterPlotChart.getTitle();
       Font jfcFont = chartTitle.getFont();
       float size = 10.0f;
       Font chartFont = jfcFont.deriveFont(size);
       chartTitle.setFont(chartFont);
       igvScatterPlotChart.setTitle(chartTitle);
        */

        // Alternate method - creat subtiltle and replace title by using null
        // for ChartFactory.createScatterPlot(title, ....
        TextTitle chartTitle = new TextTitle(title, new Font("Tahoma", Font.BOLD, 10));
        igvScatterPlotChart.addSubtitle(chartTitle);

        // get renderer and set series plot shape and color  from seriesSettings
        XYPlot plot = igvScatterPlotChart.getXYPlot();
        //XYItemRenderer renderer = plot.getRenderer();
        //Class renderClass = renderer.getClass();    // returns an XYLineAndShapeRenderer
        XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) plot.getRenderer();

        //  Note seriesCount is JFreeChart series count; nSeries is SymbolSettings items
        int seriesCount = igvXYSeriesDataset.getSeriesCount();

        // set the up series symbol settings for plotting
        for (int series = 0; series < nSeries; ++series) {
            Shape shape = symbolSettings[series].getPlotShape();
            renderer.setSeriesShape(series, shape);
            Color color = symbolSettings[series].getPlotColor();
            renderer.setSeriesPaint(series, color);
            if (symbolSettings[series].isFilled()) {
                //renderer.setSeriesFillPaint(series, color);
                //renderer.setUseFillPaint(true);
            }
            if (symbolSettings[series].isOutlined()) {
                //renderer.setSeriesOutlinePaint(series, Color.black);
                //renderer.setDrawOutlines(true);
            }
            boolean visible = symbolSettings[series].isVisible();
            renderer.setSeriesVisible(series, visible);
        }

        // modify the legend to include the selected attribute mapped to series plots
        LegendTitle legend = new LegendTitle(igvScatterPlotChart.getPlot());

        BlockContainer wrapper = new BlockContainer(new BorderArrangement());
        wrapper.setFrame(new BlockBorder(1.0, 1.0, 1.0, 1.0));
        LabelBlock legendTitle = new LabelBlock(attribute + ":");
        legendTitle.setPadding(5, 5, 5, 5);
        wrapper.add(legendTitle, RectangleEdge.TOP);
        BlockContainer items = legend.getItemContainer();
        items.setPadding(2, 10, 5, 2);
        wrapper.add(items);
        legend.setWrapper(wrapper);

        // query series settings for debug
//        for (int series = 0; series < nSeries; ++series) {
//            int itemCount = igvXYSeriesDataset.getSeries(series).getItemCount();
//            Shape shape = renderer.getSeriesShape(series);
//            Paint color = renderer.getSeriesPaint(series);
//            Paint fillcolor = renderer.getSeriesFillPaint(series);
//        }

        // tooltips generated and loaded for all samples
        createToolTips(renderer);

        // load the modified renderer
        plot.setRenderer(renderer);

        // set remaining chart attributes
        plot.setDomainCrosshairVisible(true);
        plot.setDomainCrosshairLockedOnData(true);
        plot.setRangeCrosshairVisible(true);
        plot.setRangeCrosshairLockedOnData(true);
        plot.setDomainZeroBaselineVisible(true);
        plot.setRangeZeroBaselineVisible(true);
        NumberAxis domainAxis = (NumberAxis) plot.getDomainAxis();
        domainAxis.setAutoRangeIncludesZero(false);

        legend.setPosition(RectangleEdge.RIGHT);
        legend.setHorizontalAlignment(HorizontalAlignment.LEFT);
        igvScatterPlotChart.addSubtitle(legend);

        ChartUtilities.applyCurrentTheme(igvScatterPlotChart);
        return igvScatterPlotChart;
    }

    private class IGVStandardChartTheme extends StandardChartTheme {

        public IGVStandardChartTheme(Title title) {
            super("IGVScatterPlotTheme");
            this.applyToTitle((Title) title);
        }
    }

    /*
    *   Method used to regenerate the scatterplot based on current global selection
    *  settings which include:
    *       igvXAxisIndex - which selects the X axis data sample measurement vriable
    *       igvYAxisIndex - which selects the Y axis data sample measurement vriable
    *       igvSymbolIndex - which selects the symbol attribute to define plotted series
    * */
    private void drawScatterplot() {

        // recreate the scatterplot chart
        String xAxisName = igvPlotData.getDataKeyName(igvXAxisIndex);
        String yAxisName = igvPlotData.getDataKeyName(igvYAxisIndex);
        igvAttributeName = igvPlotData.getSymbolKeyName(igvSymbolIndex);

        igvFrameContent.remove(igvAttributeLegendPanel);
        igvAttributeLegendPanel = createAttributeSymbolPanel(igvAttributeName);
        igvFrameContent.add(igvAttributeLegendPanel, BorderLayout.EAST);
        igvFrameContent.validate();

        igvScatterPlotChart = createChart(igvTitle, igvPlotData, xAxisName, yAxisName,
                igvAttributeName, igvSymbolSettings);

        igvPlotPanel.setChart(igvScatterPlotChart);
    }

    /*
    *   Method used to regenerate the scatterplot based on current global selection
    *  settings for plot series visibility.
    *
    *       igvSymbolSettings - selection settings for attribute categories for plotted series
    * */
    private void updateScatterplot() {

        // get renderer and set series plot shape and color  from seriesSettings
        XYPlot plot = igvScatterPlotChart.getXYPlot();
        XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) plot.getRenderer();

        int nSeries = igvSymbolSettings.length;

        // set the series symbol visibility
        // Note: scatter plot chart is notified of changes
        for (int series = 0; series < nSeries; ++series) {

            boolean visible = igvSymbolSettings[series].isVisible();
            renderer.setSeriesVisible(series, visible);
        }

    }

    /*
    *   createSampleListPanel - method creates a ListPanel for a scrollable list of sample names .
    *
    *   Parameters:
    *       sampleNames - array of sample names taken from IGV
    *
     *   Returns:
    *       sample names list selection panel
    * */
    /*
     private JPanel createSampleListPanel(String[] sampleNames)
     {
      JPanel samplePanel = new JPanel();

      samplePanel.setLayout( new FlowLayout() ); // set panel layout

      JList sampleList = new JList( sampleNames ); // create with sampleNames
      sampleList.setVisibleRowCount( 20 ); // display twenty rows at once

      // do not allow multiple selections
      sampleList.setSelectionMode( ListSelectionModel.SINGLE_SELECTION );

      // add a JScrollPane containing JList to frame
      samplePanel.add( new JScrollPane( sampleList ) );

      sampleList.addListSelectionListener(
         new ListSelectionListener() // anonymous inner class
         {
            // handle list selection events
            public void valueChanged( ListSelectionEvent event )
            {
               // take some action on the selected sample

            } // end method valueChanged
         } // end anonymous inner class
      ); // end call to addListSelectionListener

      return samplePanel;

     } // end ListFrame constructor
     */
}
