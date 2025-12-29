package org.igv.track;

import java.util.Comparator;

/**
 * Sort tracks by attribute value
 */
public abstract class AttributeComparator<T> implements Comparator<T> {

    private final String[] attributeNames;
    private final boolean[] ascending;

    AttributeComparator(String[] attributeNames, boolean[] ascending) {
        assert attributeNames.length == ascending.length;
        this.attributeNames = attributeNames;
        this.ascending = ascending;
    }

    protected abstract String getAttributeValue(T track, String attName);

    public int compare(T t1, T t2) {
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


    public static class TrackAttributeComparator extends AttributeComparator<Track> {

        public TrackAttributeComparator(String[] attributeNames, boolean[] ascending) {
            super(attributeNames, ascending);
        }

        protected String getAttributeValue(Track track, String attName) {
            String value = track.getAttributeValue(attName);

            if (value == null) {
                value = "";
            }

            return value.toLowerCase();
        }
    }


    /**
     * Sort samples by attribute value
     */
    public static class SampleAttributeComparator extends AttributeComparator<String> {

        public SampleAttributeComparator(String[] attributeNames, boolean[] ascending) {
            super(attributeNames, ascending);
        }

        protected String getAttributeValue(String sample, String attName) {
            String value = AttributeManager.getInstance().getAttribute(sample, attName);

            if (value == null) {
                value = "";
            }

            return value.toLowerCase();
        }

    }
}
