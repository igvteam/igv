package org.broad.igv.hic;

import javax.swing.*;
import java.awt.*;
import java.io.Serializable;

/**
 * @author jrobinso
 * @date Aug 2, 2010
 */
public class HeatmapPanel extends JComponent implements Serializable {

    private MainWindow mainWindow;
    private double maxCount = 200.0 / (1 * 1);
    private int binWidth = 2;

    HeatmapRenderer renderer = new HeatmapRenderer();


    public HeatmapPanel() {

    }

    @Override
    protected void paintComponent(Graphics g) {

        if (mainWindow != null && mainWindow.zd != null) {
            Rectangle bounds = this.getVisibleRect();
            renderer.render(mainWindow.xContext.getOrigin(), mainWindow.yContext.getOrigin(), mainWindow.zd,
                    binWidth, maxCount, g, bounds);

        }
    }

    public void setMaxCount(double maxCount) {
        this.maxCount = maxCount;
    }

    public int getBinWidth() {
        return binWidth;
    }

    public MainWindow getMainWindow() {
        return mainWindow;
    }

    public void setMainWindow(MainWindow mainWindow) {
        this.mainWindow = mainWindow;
    }

    public void setBinWidth(int binWidth) {
        this.binWidth = binWidth;
    }

    @Override
    public void setSize(int i, int i1) {
        super.setSize(i, i1);    //To change body of overridden methods use File | Settings | File Templates.
    }

    @Override
    public void setSize(Dimension dimension) {
        super.setSize(dimension);    //To change body of overridden methods use File | Settings | File Templates.
    }

    @Override
    public void setBounds(Rectangle rectangle) {
        super.setBounds(rectangle);    //To change body of overridden methods use File | Settings | File Templates.
    }

    @Override
    public void setBounds(int i, int i1, int i2, int i3) {
        super.setBounds(i, i1, i2, i3);    //To change body of overridden methods use File | Settings | File Templates.
    }
}
