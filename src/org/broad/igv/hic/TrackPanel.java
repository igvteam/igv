package org.broad.igv.hic;

import javax.swing.*;
import java.awt.*;
import java.awt.image.BufferedImage;

/**
 * Created by IntelliJ IDEA.
 * User: neva
 * Date: 4/3/12
 * Time: 4:08 PM
 * To change this template use File | Settings | File Templates.
 */
public class TrackPanel extends JPanel {

    private double[] data;
    private BufferedImage image;

    public void setData(double[] data){
        this.data = data;
        image = new BufferedImage(data.length, getPreferredSize().height, BufferedImage.TYPE_INT_ARGB);
        int h = getPreferredSize().height/2;
        Graphics2D g2d = (Graphics2D)image.getGraphics();
        g2d.setColor(Color.blue.darker());
        double data_max = 0;
        for (double aData : data) {
            if (Math.abs(aData) > data_max) data_max = Math.abs(aData);
        }
        for (int i=0; i < data.length; i++){
            int myh = (int)(data[i]*h/data_max);
            if (data[i] > 0)
                g2d.fillRect(i, h-myh, 1, myh);
            else
                g2d.fillRect(i, h, 1, -myh);
        }
    }

    @Override
    protected void paintComponent(Graphics graphics) {
        Graphics2D g2d = (Graphics2D)graphics;
        super.paintComponent(g2d);
        if (image != null) {
            BufferedImage image2 = new BufferedImage(getWidth(), getHeight(), image.getType());
            Graphics2D g = image2.createGraphics();
            g.drawImage(image, 0, 0, getWidth(), getHeight(), null);
            g.dispose();
            g.setComposite(AlphaComposite.Src);

            g.setRenderingHint(RenderingHints.KEY_INTERPOLATION,
                    RenderingHints.VALUE_INTERPOLATION_BILINEAR);
            g.setRenderingHint(RenderingHints.KEY_RENDERING,
                    RenderingHints.VALUE_RENDER_QUALITY);
            g.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
                    RenderingHints.VALUE_ANTIALIAS_ON);
            
            graphics.drawImage(image2, 0, 0, null);
        }
    }
}
