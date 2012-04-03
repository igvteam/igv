package org.broad.igv.hic;

import javax.swing.*;
import java.awt.*;

/**
 * Created by IntelliJ IDEA.
 * User: neva
 * Date: 4/3/12
 * Time: 4:08 PM
 * To change this template use File | Settings | File Templates.
 */
public class TrackPanel extends JPanel {

    private double[] data;


    public void setData(double[] data){
        this.data = data;
    }

    @Override
    public void paintComponent(Graphics graphics) {
        super.paintComponent(graphics);
        if (data != null) {
            int h = getHeight()/2;

            int w = getWidth()/data.length;
            System.out.println(getWidth() + " " + w);
            graphics.setColor(Color.blue.darker());
            double data_max = 0;
            for (double aData : data) {
                if (Math.abs(aData) > data_max) data_max = Math.abs(aData);
            }

            for (int i=0; i < data.length; i++){
                int myh = (int)(data[i]*h/data_max);
                System.out.println(myh + " " + data[i]);
                if (data[i] > 0)
                    graphics.fillRect(i*w, h, w, myh);
                else
                    graphics.fillRect(i*w, h, w, myh);
            }
        }

    }
}
