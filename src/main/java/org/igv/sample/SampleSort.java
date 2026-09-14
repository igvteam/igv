package org.igv.sample;

import org.json.JSONArray;
import org.json.JSONObject;

/**
 * A sort of a track's samples, saved in sessions as the track "sort" property.  The json form matches igv.js:
 * <pre>
 *   {"option": "GENOTYPE", "direction": "DESC", "chr": "chr1", "start": 100, "end": 101}
 *   {"option": "ATTRIBUTE", "attribute": "Population", "direction": "ASC"}
 * </pre>
 * Sorts igv.js does not support yet are written in the same form: options SAMPLE_NAME, DEPTH, and QUALITY, and
 * attribute sorts on more than one attribute, which write "attribute" and "direction" as parallel arrays.
 */
public class SampleSort {

    public static final String GENOTYPE = "GENOTYPE";
    public static final String DEPTH = "DEPTH";
    public static final String QUALITY = "QUALITY";
    public static final String VALUE = "VALUE";
    public static final String ATTRIBUTE = "ATTRIBUTE";
    public static final String SAMPLE_NAME = "SAMPLE_NAME";

    private final String option;          // null if not specified -- igv.js treats this as the track's default sort
    private final boolean ascending;
    private final String chr;
    private final int start;
    private final int end;
    private final String[] attributes;
    private final boolean[] attributeAscending;

    private SampleSort(String option, boolean ascending, String chr, int start, int end,
                       String[] attributes, boolean[] attributeAscending) {
        this.option = option;
        this.ascending = ascending;
        this.chr = chr;
        this.start = start;
        this.end = end;
        this.attributes = attributes;
        this.attributeAscending = attributeAscending;
    }

    /**
     * A sort by data over a genomic region -- GENOTYPE, DEPTH, QUALITY, or VALUE.
     */
    public static SampleSort locus(String option, String chr, int start, int end, boolean ascending) {
        return new SampleSort(option, ascending, chr, start, end, null, null);
    }

    public static SampleSort sampleName(boolean ascending) {
        return new SampleSort(SAMPLE_NAME, ascending, null, 0, 0, null, null);
    }

    public static SampleSort attributes(String[] attributes, boolean[] ascending) {
        return new SampleSort(ATTRIBUTE, ascending.length > 0 && ascending[0], null, 0, 0, attributes, ascending);
    }

    public String getOption() {
        return option;
    }

    public boolean isAscending() {
        return ascending;
    }

    public String getChr() {
        return chr;
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }

    public String[] getAttributes() {
        return attributes;
    }

    public boolean[] getAttributeAscending() {
        return attributeAscending;
    }

    public JSONObject toJson() {
        JSONObject json = new JSONObject();
        json.put("option", option);
        if (ATTRIBUTE.equals(option)) {
            if (attributes.length == 1) {
                json.put("attribute", attributes[0]);
                json.put("direction", direction(attributeAscending[0]));
            } else {
                JSONArray names = new JSONArray();
                JSONArray directions = new JSONArray();
                for (int i = 0; i < attributes.length; i++) {
                    names.put(attributes[i]);
                    directions.put(direction(attributeAscending[i]));
                }
                json.put("attribute", names);
                json.put("direction", directions);
            }
        } else {
            json.put("direction", direction(ascending));
            if (chr != null) {
                json.put("chr", chr);
                json.put("start", start);
                json.put("end", end);
            }
        }
        return json;
    }

    public static SampleSort fromJson(JSONObject json) {

        String option = json.has("option") ? json.getString("option").toUpperCase() : null;

        if (ATTRIBUTE.equals(option)) {
            JSONArray names = json.optJSONArray("attribute");
            if (names == null) {
                return attributes(new String[]{json.getString("attribute")},
                        new boolean[]{isAscending(json.optString("direction"))});
            }
            JSONArray directions = json.optJSONArray("direction");
            String[] attributes = new String[names.length()];
            boolean[] ascending = new boolean[names.length()];
            for (int i = 0; i < names.length(); i++) {
                attributes[i] = names.getString(i);
                ascending[i] = directions != null && isAscending(directions.optString(i));
            }
            return attributes(attributes, ascending);
        }

        boolean ascending = isAscending(json.optString("direction"));
        if (SAMPLE_NAME.equals(option)) {
            return sampleName(ascending);
        }

        // igv.js accepts a 1-based "position" in place of start and end
        int start = json.has("start") ? json.getInt("start") : json.getInt("position") - 1;
        int end = json.has("end") ? json.getInt("end") : start + 1;
        return locus(option, json.getString("chr"), start, end, ascending);
    }

    private static String direction(boolean ascending) {
        return ascending ? "ASC" : "DESC";
    }

    private static boolean isAscending(String direction) {
        return "ASC".equalsIgnoreCase(direction);
    }
}
