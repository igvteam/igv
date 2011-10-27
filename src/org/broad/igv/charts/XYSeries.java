package org.broad.igv.charts;

import org.broad.igv.util.collections.DoubleArrayList;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Jim Robinson
 * @date 10/26/11
 */
public class XYSeries {

    String name;
    DoubleArrayList xData;
    DoubleArrayList yData;
    List<String> description;

    public XYSeries(String name) {
        this.name = name;
        xData = new DoubleArrayList(500);
        yData = new DoubleArrayList(500);
        description = new ArrayList<String>(500);
    }

    public void add(double x, double y, String description) {
        xData.add(x);
        yData.add(y);
        this.description.add(description);
    }

    public int getSize() {
        return xData.size();
    }

    public String getName() {
        return name;
    }

    double [] getX() {
        return xData.toArray();
    }

    double [] getY() {
        return yData.toArray();
    }

}


