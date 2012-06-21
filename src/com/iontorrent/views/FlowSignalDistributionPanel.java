/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package com.iontorrent.views;

import com.iontorrent.data.FlowDistribution;
import com.iontorrent.utils.FileTools;
import com.iontorrent.utils.LocationListener;
import java.awt.*;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.StringSelection;
import java.awt.font.FontRenderContext;
import java.awt.font.GlyphVector;
import java.awt.geom.AffineTransform;
import java.awt.geom.GeneralPath;
import java.io.File;
import java.io.IOException;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.*;
import org.jfree.chart.ChartFactory;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.NumberTickUnit;
import org.jfree.chart.plot.*;
import org.jfree.chart.renderer.category.BarRenderer;
import org.jfree.chart.renderer.category.CategoryItemRenderer;
import org.jfree.chart.renderer.category.StackedBarRenderer;
import org.jfree.chart.renderer.xy.*;
import org.jfree.chart.title.TextTitle;
import org.jfree.data.category.CategoryDataset;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

/**
 *
 * @author Chantal Roth
 */
public class FlowSignalDistributionPanel extends javax.swing.JPanel {

    private static org.apache.log4j.Logger log = org.apache.log4j.Logger.getLogger(FlowSignalDistributionPanel.class);
    /**
     * bar or line chart - only static until we move it into user preferences
     */
    private static int chart_type = -1;
    /**
     * The data: key is the flow signal value, such as 654, and the value is how
     * often it was found
     */
    private FlowDistribution distributions[];
    /**
     * current chromosome location
     */
    private int location;
    /**
     * the size of the bin - the map has bin size 1, which is usually too small
     * for anything to see. This is customizable by the user. Should be stored
     * in some properties... for now we make it static so that it is remembered
     * between calls
     */
    private static int binsize;
    /**
     * where user can store .csv data
     */
    private String filename;
    /**
     * jfreechart
     */
    private JComponent chartpanel;
    private JFreeChart freechart;
    /**
     * for navigation left and right
     */
    private LocationListener listener;
    private Font font = new Font("SansSerif", Font.PLAIN, 10);
    /**
     * Chart colors and symbols
     */
    static final String BASES = "ACGT";
    Shape shapea = getGlyphShape("A");
    Shape shapec = getGlyphShape("C");
    Shape shapeg = getGlyphShape("G");
    Shape shapet = getGlyphShape("T");
    Shape shapex = getGlyphShape("?");
    //Color colA = Color.GREEN.darker();
    Color colA = new Color(0, 150, 30);
    Color colC = Color.BLUE.darker();
    // Color colT = Color.RED.darker();
    Color colT = Color.RED.darker();
    Color colG = new Color(30, 30, 30);
    Color colX = Color.orange;
    Color colors[] = {colA, colC, colG, colT, colX};
    Shape shapes[] = {shapea, shapec, shapeg, shapet, shapex};
    float dash[] = {8.0f};
    BasicStroke dashedStroke = new BasicStroke(2.0f, BasicStroke.CAP_BUTT,
            BasicStroke.JOIN_MITER, 8.0f, dash, 0.0f);

    /** max x val in chart */
    private int maxx;
    public FlowSignalDistributionPanel(FlowDistribution distributions[]) {
        this.distributions = distributions;
        if (chart_type < 0) {
            chart_type = ChartConfigPanel.TYPE_LINE;
        }
        if (binsize < 1) {
            binsize = 25;
        }
        
        initComponents();
        recreateChart();
    }

    private void recreateChart() {
        if (distributions == null || distributions.length < 1) {
            JOptionPane.showMessageDialog(this, "I got no flow signal distribution data");
            return;
        }
        if (chartpanel != null) {
            remove(chartpanel);
        }
        chartpanel = createChart();
        add("Center", chartpanel);
        // with freechart, one sometimes just doesn't get a repaint... 
        this.invalidate();
        chartpanel.invalidate();
        chartpanel.repaint();
        this.repaint();
        this.paintAll(getGraphics());



    }

    private CategoryDataset createCategoryDataset() {
        DefaultCategoryDataset dataset = new DefaultCategoryDataset();
        for (int i = 0; i < distributions.length; i++) {
            int[] data = distributions[i].getBinnedData(binsize);
            String seriename = distributions[i].getName();
            for (int j = 0; j < data.length; j++) {
                String cat = "" + (binsize * j);
                dataset.addValue(data[j], seriename, cat);
            }

        }
        return dataset;
    }

    private XYSeriesCollection createXYDataset() {
        XYSeriesCollection dataset = new XYSeriesCollection();
        for (int i = 0; i < distributions.length; i++) {
            TreeMap<Short, Integer> map = distributions[i].getMap();
            dataset.addSeries(createDataset(distributions[i]));
        }
        return dataset;
    }

    private JComponent createChart() {
        maxx = 0;
        this.location = distributions[0].getLocation();
        String information = distributions[0].getInformation();
        String plotTitle = "Flow Signal Distribution";
        String xaxis = "flow signal value";
        String yaxis = "count";


        NumberAxis xax = new NumberAxis(xaxis);
        NumberAxis yax = new NumberAxis(yaxis);
        freechart = null;
        // could be another chart option
        // XYSplineRenderer renderer = new XYSplineRenderer();

        if (chart_type == ChartConfigPanel.TYPE_BAR || chart_type == ChartConfigPanel.TYPE_STACKED) {
            CategoryDataset dataset = createCategoryDataset();
            freechart = ChartFactory.createBarChart(
                    plotTitle, // chart title
                    xaxis, // domain axis label
                    yaxis, // range axis label
                    dataset, // data
                    PlotOrientation.VERTICAL, // orientation
                    true, // include legend
                    true, // tooltips?
                    false // URLs?
                    );

            CategoryPlot plot = freechart.getCategoryPlot();
            plot.setForegroundAlpha(0.5f);
            CategoryItemRenderer renderer;
            if (chart_type == ChartConfigPanel.TYPE_STACKED) {
                renderer = new StackedBarRenderer();
                plot.setRenderer(renderer);
            } else {
                BarRenderer bar = (BarRenderer) plot.getRenderer();
                bar.setShadowVisible(false);
                bar.setItemMargin(0);
                bar.setDrawBarOutline(false);
                renderer = bar;
            }

            int series = dataset.getRowCount();
            for (int s = 0; s < series; s++) {
                char base = distributions[s].getBase();
                int which = BASES.indexOf(base);
                if (which < 0) {
                    which = 4;
                }
                Color c = colors[which];

                if (distributions[s].isReverse() && !distributions[s].isForward()) {
                    renderer.setSeriesOutlineStroke(s, dashedStroke);
                    renderer.setSeriesStroke(s, dashedStroke);
                    c = c.darker();
                } else if (!distributions[s].isReverse() && distributions[s].isForward()) {
                    c = c.brighter();
                }
                renderer.setSeriesOutlinePaint(s, c);
                renderer.setSeriesItemLabelPaint(s, c);
                renderer.setSeriesPaint(s, c);

            }
        } else {
            XYSeriesCollection dataset = null;
            AbstractXYItemRenderer renderer = null;
            dataset = createXYDataset();
            XYPlot plot = null;
            if (chart_type == ChartConfigPanel.TYPE_LINE) {
                renderer = new XYLineAndShapeRenderer();                
                plot = new XYPlot(dataset, xax, yax, renderer);
                freechart = new JFreeChart(plot);
                plot.setForegroundAlpha(0.75f);

            } else if (chart_type == ChartConfigPanel.TYPE_AREA) {
                renderer = new XYAreaRenderer(XYAreaRenderer.AREA_AND_SHAPES);
                plot = new XYPlot(dataset, xax, yax, renderer);
                plot.setForegroundAlpha(0.5f);
                freechart = new JFreeChart(plot);
            }
            plot.setDomainGridlinesVisible(true);
            int series = dataset.getSeriesCount();
          //  plot.setDomainGridlinePaint(Color.gray.darker());  
            for (int x = 100; x < maxx; x+= 100) {
                final Marker line = new ValueMarker(x);
                line.setPaint(Color.gray.darker());            
                plot.addDomainMarker(line);
            }
            
        
            for (int s = 0; s < series; s++) {
                char base = distributions[s].getBase();
                int which = BASES.indexOf(base);
                if (which < 0) {
                    which = 4;
                }
                Color c = colors[which];

//                if (renderer instanceof XYLineAndShapeRenderer) {
//                    XYLineAndShapeRenderer r = (XYLineAndShapeRenderer) renderer;
//                    Shape shape = shapes[which];
//                    //  p("Using shape " + shape.getClass().getName() + " for base " + base);
//                    r.setBaseShapesVisible(true);
//                    r.setSeriesShape(s, shape);
//                    r.setSeriesShapesVisible(s, true);
//                }
                // if reverse and not forward, use dashed line
                if (distributions[s].isReverse() && !distributions[s].isForward()) {
                    renderer.setSeriesOutlineStroke(s, dashedStroke);
                    renderer.setSeriesStroke(s, dashedStroke);
                    c = c.darker();
                } else if (!distributions[s].isReverse() && distributions[s].isForward()) {
                    c = c.brighter();
                }
                renderer.setSeriesFillPaint(s, c);
                renderer.setSeriesOutlinePaint(s, c);
                renderer.setSeriesItemLabelPaint(s, c);
                renderer.setSeriesPaint(s, c);
            }
        }
        
        xax.setMinorTickCount(1);
        xax.setTickUnit(new NumberTickUnit(50));
        
        xax.setTickMarksVisible(true);
        xax.setMinorTickMarksVisible(true);
        freechart.setTitle(plotTitle);
        String bininfo = "bin size=" + binsize;
        if (information
                == null) {
            freechart.addSubtitle(new TextTitle(bininfo));
        } else {
            freechart.addSubtitle(new TextTitle(information + ", " + bininfo));
        }
        ChartPanel chartPanel = new ChartPanel(freechart);

        chartPanel.setPreferredSize(
                new java.awt.Dimension(800, 600));
        return chartPanel;
    }

    public Shape getGlyphShape(String strGlyphs) {
        AffineTransform transform = new AffineTransform();
//        transform.translate(-6, 6);
        transform.scale(1, 1);
        return getGlyphShapes(font, strGlyphs, transform)[0];
    }

    public Shape[] getGlyphShapes(Font font, String strGlyphs, AffineTransform transform) {

        FontRenderContext frc = new FontRenderContext(null, true, true);
        GlyphVector glyphs = font.createGlyphVector(frc, strGlyphs);

        int count = glyphs.getNumGlyphs();
        Shape[] shapes = new Shape[count];
        for (int i = 0; i < count; i++) {

            // get transformed glyph shape
            GeneralPath path = (GeneralPath) glyphs.getGlyphOutline(i);
            shapes[i] = path.createTransformedShape(transform);
        }
        return shapes;

    }

    private XYSeries createDataset(FlowDistribution dist) {
        int[] bins = dist.getBinnedData(binsize);
        if (maxx < dist.getMaxX()) maxx = dist.getMaxX();
        XYSeries xy = new XYSeries(dist.getName());
        for (int b = 0; b < bins.length; b++) {
            xy.add(b * binsize, bins[b]);
        }
        return xy;
    }

    private String getCsvString() {
        StringBuilder csv = new StringBuilder();
        //csv = csv.append(information).append("\n\n");
        for (int i = 0; i < distributions.length; i++) {
            csv = csv.append(distributions[i].toCsv(binsize));
        }
        return csv.toString();
    }

    private String getJsonString() {
        StringBuilder json = new StringBuilder();
        //json = json.append(information).append("\n\n");
        for (int i = 0; i < distributions.length; i++) {
            json = json.append("\n").append(distributions[i].toJson());
        }
        return json.toString();
    }

    private String getReadString() {
        StringBuilder rinfo = new StringBuilder();
        for (int i = 0; i < distributions.length; i++) {
            rinfo = rinfo.append("\n").append(distributions[i].getReadInfoString());
        }
        return rinfo.toString();
    }

    private void p(String msg) {
        log.info(msg);
    }

    private void err(String msg) {
        log.error(msg);
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        buttonToolBar = new javax.swing.JToolBar();
        btnCopy = new javax.swing.JButton();
        btnCopyJson = new javax.swing.JButton();
        jButton1 = new javax.swing.JButton();
        btnSave = new javax.swing.JButton();
        btnConfigure = new javax.swing.JButton();
        btnLeft = new javax.swing.JButton();
        btnRight = new javax.swing.JButton();

        setLayout(new java.awt.BorderLayout());

        buttonToolBar.setRollover(true);

        btnCopy.setIcon(new javax.swing.ImageIcon(getClass().getResource("/com/iontorrent/views/copy.png"))); // NOI18N
        btnCopy.setToolTipText("Copy the data to the clipboard");
        btnCopy.setFocusable(false);
        btnCopy.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        btnCopy.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        btnCopy.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnCopyActionPerformed(evt);
            }
        });
        buttonToolBar.add(btnCopy);

        btnCopyJson.setIcon(new javax.swing.ImageIcon(getClass().getResource("/com/iontorrent/views/copyj.png"))); // NOI18N
        btnCopyJson.setToolTipText("copy to clip board in Json format");
        btnCopyJson.setFocusable(false);
        btnCopyJson.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        btnCopyJson.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        btnCopyJson.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnCopyJsonActionPerformed(evt);
            }
        });
        buttonToolBar.add(btnCopyJson);

        jButton1.setIcon(new javax.swing.ImageIcon(getClass().getResource("/com/iontorrent/views/copyr.png"))); // NOI18N
        jButton1.setToolTipText("show read info");
        jButton1.setFocusable(false);
        jButton1.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        jButton1.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        jButton1.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jButton1ActionPerformed(evt);
            }
        });
        buttonToolBar.add(jButton1);

        btnSave.setIcon(new javax.swing.ImageIcon(getClass().getResource("/com/iontorrent/views/save.png"))); // NOI18N
        btnSave.setToolTipText("Save data in .csv file (to be used in Excel for instance)");
        btnSave.setFocusable(false);
        btnSave.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        btnSave.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        btnSave.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnSaveActionPerformed(evt);
            }
        });
        buttonToolBar.add(btnSave);

        btnConfigure.setIcon(new javax.swing.ImageIcon(getClass().getResource("/com/iontorrent/views/configure.png"))); // NOI18N
        btnConfigure.setToolTipText("Change the bin size");
        btnConfigure.setFocusable(false);
        btnConfigure.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        btnConfigure.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        btnConfigure.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnConfigureActionPerformed(evt);
            }
        });
        buttonToolBar.add(btnConfigure);

        btnLeft.setIcon(new javax.swing.ImageIcon(getClass().getResource("/com/iontorrent/views/arrow-left.png"))); // NOI18N
        btnLeft.setToolTipText("move to the next bas on the left");
        btnLeft.setFocusable(false);
        btnLeft.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        btnLeft.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        btnLeft.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnLeftActionPerformed(evt);
            }
        });
        buttonToolBar.add(btnLeft);

        btnRight.setIcon(new javax.swing.ImageIcon(getClass().getResource("/com/iontorrent/views/arrow-right.png"))); // NOI18N
        btnRight.setToolTipText("move to the next base on the right");
        btnRight.setFocusable(false);
        btnRight.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        btnRight.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        btnRight.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnRightActionPerformed(evt);
            }
        });
        buttonToolBar.add(btnRight);

        add(buttonToolBar, java.awt.BorderLayout.PAGE_START);
    }// </editor-fold>//GEN-END:initComponents

    private void btnConfigureActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnConfigureActionPerformed
        // let user pick bin size
        ChartConfigPanel config = new ChartConfigPanel();
        config.setChartType(chart_type);
        config.setBinSize(binsize);
        int ans = JOptionPane.showConfirmDialog(this, config,
                "Chart Option", JOptionPane.OK_CANCEL_OPTION);
        if (ans != JOptionPane.OK_OPTION) {
            return;
        }
        binsize = config.getBinSize();
        chart_type = config.getChartType();
        recreateChart();
    }//GEN-LAST:event_btnConfigureActionPerformed

    /**
     * Convert to export to Excel
     */
    private void btnSaveActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnSaveActionPerformed


        String[] options = new String[2];
        String msg = "1) Save image of chart\n";
        msg += "2) Save data in .csv file\n";
        options[0] = "1) Image";
        options[1] = "2) Data";

        int ans = JOptionPane.showOptionDialog(this, msg, "Export",
                JOptionPane.OK_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE, null, options, options[0]);
        if (ans < 0) {
            return;
        }
        if (ans == 0) {
            doSaveImageAction();
        } else if (ans == 1) {
            doSaveDataAction();
        }


    }//GEN-LAST:event_btnSaveActionPerformed
    private void doSaveImageAction() {
        filename = FileTools.getFile("File to store chart image", ".png", filename, true);
        if (filename == null || filename.length() < 1) {
            return;
        }
        try {
            ChartUtilities.saveChartAsJPEG(new File(filename), freechart, 800, 600);
        } catch (IOException ex) {
            Logger.getLogger(FlowSignalDistributionPanel.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    private void doSaveDataAction() {
        // get file name from user
        filename = FileTools.getFile("File to store flow chart distribution", ".csv", filename, true);
        if (filename == null || filename.length() < 1) {
            return;
        }
        String csv = getCsvString();

        File fileToSave = new File(filename);
        FileTools.writeStringToFile(fileToSave, csv, false);
    }

    private void btnCopyActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnCopyActionPerformed
        String csv = getCsvString();
        // copy to clipboard
        StringSelection stringSelection = new StringSelection(csv);
        Clipboard clipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
        clipboard.setContents(stringSelection, null);

        // also show as popup for copy paste
        JTextArea area = new JTextArea(15, 30);
        area.setText(csv);
        JOptionPane.showMessageDialog(this, new JScrollPane(area), "Data to paste to Excel (it is already in the clipboard)", JOptionPane.INFORMATION_MESSAGE);
    }//GEN-LAST:event_btnCopyActionPerformed

    private void btnCopyJsonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnCopyJsonActionPerformed
        String json = getJsonString();
        // copy to clipboard
        StringSelection stringSelection = new StringSelection(json);
        Clipboard clipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
        clipboard.setContents(stringSelection, null);

        // also show as popup for copy paste
        JTextArea area = new JTextArea(15, 30);
        area.setText(json);
        JOptionPane.showMessageDialog(this, new JScrollPane(area), "Json string (it is already in the clipboard)", JOptionPane.INFORMATION_MESSAGE);
    }//GEN-LAST:event_btnCopyJsonActionPerformed

    private void btnLeftActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnLeftActionPerformed
        if (getListener() != null) {
            getListener().locationChanged(location - 1);
        }
    }//GEN-LAST:event_btnLeftActionPerformed

    private void btnRightActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnRightActionPerformed
        if (getListener() != null) {
            getListener().locationChanged(location + 1);
        }
    }//GEN-LAST:event_btnRightActionPerformed

    private void jButton1ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jButton1ActionPerformed
        String readinfo = getReadString();
        // copy to clipboard
        StringSelection stringSelection = new StringSelection(readinfo);
        Clipboard clipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
        clipboard.setContents(stringSelection, null);

        // also show as popup for copy paste
        JTextArea area = new JTextArea(15, 30);
        area.setText(readinfo);
        JOptionPane.showMessageDialog(this, new JScrollPane(area), "Read info string (it is already in the clipboard)", JOptionPane.INFORMATION_MESSAGE);
    }//GEN-LAST:event_jButton1ActionPerformed
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton btnConfigure;
    private javax.swing.JButton btnCopy;
    private javax.swing.JButton btnCopyJson;
    private javax.swing.JButton btnLeft;
    private javax.swing.JButton btnRight;
    private javax.swing.JButton btnSave;
    private javax.swing.JToolBar buttonToolBar;
    private javax.swing.JButton jButton1;
    // End of variables declaration//GEN-END:variables

    public void setDistributions(FlowDistribution[] newdist) {
        this.distributions = newdist;
        recreateChart();
    }

    /**
     * @return the listener
     */
    public LocationListener getListener() {
        return listener;
    }

    /**
     * @param listener the listener to set
     */
    public void setListener(LocationListener listener) {
        this.listener = listener;
    }
}
