
package org.igv.sample;

import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.session.SessionAttribute;
import org.igv.track.Track;
import org.igv.ui.IGV;
import org.igv.util.FilterElement;
import org.json.JSONArray;
import org.json.JSONObject;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;


public class SampleFilter {

    private static Logger log = LogManager.getLogger(SampleFilter.class);

    private LinkedHashSet<FilterElement> elements;
    private boolean isEnabled = true;

    boolean matchAll = true; // If true, all elements must match for the track to be visible. If flase, any element match will make the track visible.


    public SampleFilter() {
        this.elements = new LinkedHashSet<>();
    }

    public SampleFilter(boolean matchAll, List<FilterElement> elements) {
        this();
        this.matchAll = matchAll;
        if (elements != null) {
            this.elements.addAll(elements);
        }
    }

    public boolean isMatchAll() {
        return matchAll;
    }

    public void removeAll() {
        elements.clear();
    }

    public void setEnabled(boolean value) {
        isEnabled = value;
    }

    public boolean isEnabled() {
        return isEnabled;
    }

    public boolean isEmpty() {
        return elements.isEmpty();
    }

    public Iterator getFilterElements() {
        return elements.iterator();
    }

    public void add(FilterElement element) {
        elements.add(element);
    }

    public void remove(FilterElement element) {
        elements.remove(element);
    }

    /**
     * Evaluate the FilterElement set.
     *
     * @return
     */
    public void evaluate() {


        boolean filterEnabled = isEnabled();

        for (Track track : IGV.getInstance().getAllTracks()) {

            // If filter is not enabled or has no elements just show all tracks
            if (!filterEnabled || elements.isEmpty()) {
                track.setVisible(true);
                continue;
            }

            if (track.isFilterable()) {
                boolean result = matchAll;
                for (FilterElement element : elements) {
                    boolean elementResult = element.evaluate(track);
                    if (matchAll && !elementResult) {
                        result = false; // If matchAll and any element does not match, the track is not visible
                        break;
                    } else if (!matchAll && elementResult) {
                        result = true; // If matchAny and any element matches, the track is visible
                        break;
                    }
                }
                track.setVisible(result);
            }
        }
    }

    public List<String> evaluateSamples(List<String> sampleNames) {

        List<String> filteredSamples = new java.util.ArrayList<>();

        for (String sampleName : sampleNames) {

            boolean result = matchAll;
            for (FilterElement element : elements) {
                boolean elementResult = element.evaluateSample(sampleName);
                if (matchAll && !elementResult) {
                    result = false; // If matchAll and any element does not match, the track is not visible
                    break;
                } else if (!matchAll && elementResult) {
                    result = true; // If matchAny and any element matches, the track is visible
                    break;
                }
            }
            if (result) {
                filteredSamples.add(sampleName);
            }
        }

        return filteredSamples;
    }


    public JSONObject toJson() {

        JSONObject filterJson = new JSONObject();
        if (!IGV.getInstance().isFilterMatchAll()) {
            filterJson.put(SessionAttribute.FILTER_MATCH, "any");
        } else {    // Defaults to match all
            filterJson.put(SessionAttribute.FILTER_MATCH, "all");
        }

        JSONArray filterElements = new JSONArray();
        Iterator iterator = getFilterElements();
        while (iterator.hasNext()) {

            FilterElement trackFilterElement = (FilterElement) iterator.next();

            JSONObject filterElementJson = new JSONObject();
            filterElementJson.put(SessionAttribute.ITEM, trackFilterElement.getAttributeKey());
            filterElementJson.put(SessionAttribute.OPERATOR, trackFilterElement.getComparisonOperator().getValue());
            filterElementJson.put(SessionAttribute.VALUE, trackFilterElement.getValue());

            filterElements.put(filterElementJson);
        }
        filterJson.put("elements", filterElements);
        return filterJson;
    }

    public static SampleFilter fromJson(JSONObject filterJson) {

        // Determine match mode (all or any)
        String matchMode = filterJson.optString(SessionAttribute.FILTER_MATCH, "all");
        boolean matchAll = "all".equalsIgnoreCase(matchMode);

        // Read filter elements
        List<FilterElement> filterElements = new ArrayList<>();
        if (filterJson.has("elements")) {
            JSONArray elementsArray = filterJson.getJSONArray("elements");
            for (int i = 0; i < elementsArray.length(); i++) {
                JSONObject elementJson = elementsArray.getJSONObject(i);

                String attributeKey = elementJson.getString(SessionAttribute.ITEM);
                String operatorString = elementJson.getString(SessionAttribute.OPERATOR);
                String value = elementJson.getString(SessionAttribute.VALUE);

                // Convert operator string to Operator enum
                FilterElement.Operator operator = getOperatorFromValue(operatorString);
                if (operator != null) {
                    filterElements.add(new FilterElement(attributeKey, operator, value));
                } else {
                    log.warn("Unknown filter operator: " + operatorString);
                }
            }
        }
        return new SampleFilter(matchAll, filterElements);
    }

    /**
     * Convert an operator value string to the corresponding Operator enum.
     */
    private static FilterElement.Operator getOperatorFromValue(String value) {
        for (FilterElement.Operator op : FilterElement.Operator.values()) {
            if (op.getValue().equalsIgnoreCase(value)) {
                return op;
            }
        }
        return null;
    }


    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("Filter[");
        sb.append("enabled=").append(isEnabled);
        sb.append(", matchAll=").append(matchAll);
        sb.append(", elements=[");
        for (FilterElement element : elements) {
            sb.append(element.toString()).append(", ");
        }
        if (!elements.isEmpty()) {
            sb.setLength(sb.length() - 2); // Remove trailing comma and space
        }
        sb.append("]");
        sb.append("]");
        return sb.toString();
    }

}
