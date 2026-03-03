package org.igv.sample;

import org.igv.track.AttributeManager;

import java.util.Comparator;

/**
 * Sort samples by attribute value
 */
public class SampleAttributeComparator implements Comparator<String> {

    private final String[] attributeNames;
    private final boolean[] ascending;

    public SampleAttributeComparator(String[] attributeNames, boolean[] ascending) {
        assert attributeNames.length == ascending.length;
        this.attributeNames = attributeNames;
        this.ascending = ascending;
    }


    protected String getAttributeValue(String sample, String attName) {
        String value = AttributeManager.getInstance().getAttribute(sample, attName);
        if (value == null) {
            value = "";
        }
        return value.toLowerCase();
    }

    public int compare(String t1, String t2) {

        // Loop through the attributes in order (primary, secondary, tertiary, ...).  The
        // first attribute to yield a non-zero comparison wins

        for (int i = 0; i < attributeNames.length; i++) {
            String attName = attributeNames[i];

            if (attName != null) {
                String value1 = getAttributeValue(t1, attName);
                String value2 = getAttributeValue(t2, attName);

                boolean isNumeric = AttributeManager.getInstance().isNumeric(attName);

                int c = 0;
                if (isNumeric) {
                    double d1;
                    try {
                        d1 = Double.parseDouble(value1);
                    } catch (NumberFormatException e) {
                        d1 = Double.MIN_VALUE;
                    }
                    double d2;
                    try {
                        d2 = Double.parseDouble(value2);
                    } catch (NumberFormatException e) {
                        d2 = Double.MIN_VALUE;
                    }
                    c = Double.compare(d1, d2);
                } else {
                    c = value1.compareTo(value2);
                }

                if (c != 0) {
                    return ascending[i] ? c : -c;
                }
            }
        }

        // All compares are equal
        return 0;
    }

}
