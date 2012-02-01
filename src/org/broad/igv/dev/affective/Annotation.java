package org.broad.igv.dev.affective;

import java.util.Map;

/**
 * @author Jim Robinson
 * @date 1/23/12
 */
public class Annotation {

    String date;
    int startTime;
    int endTime;
    String type;
    String description;
    String descriptiveText;

    public Annotation(String date, int startTime, int endTime, String type, String description,
                      String [] attributes, String [] values) {
        this.date = date;
        this.description = description;
        this.endTime = endTime;
        this.startTime = startTime;
        this.type = type;

        StringBuffer buf = new StringBuffer();
        buf.append("<b>" + description + "</b><br>");
        for(int i=1; i<attributes.length; i++) {
            buf.append(attributes[i] + ": " + values[i] + "<br>");
        }
        descriptiveText = buf.toString();
    }
}
