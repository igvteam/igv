/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broad.igv.util;

import org.broad.igv.track.Track;
import org.broad.igv.ui.IGV;

import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;

/**
 * @author eflakes
 */
public class TrackFilter {

    private LinkedHashSet<FilterElement> elements;
    private boolean isEnabled = true;

    boolean showAll = false; // If true, all tracks are shown regardless of filter
    boolean matchAll = true; // If true, all elements must match for the track to be visible. If flase, any element match will make the track visible.


    public TrackFilter() {
        this.elements = new LinkedHashSet<>();
    }

    public TrackFilter(boolean showAll, boolean matchAll, List<FilterElement> elements) {
        this();
        this.showAll = showAll;
        this.matchAll = matchAll;
        if (elements != null) {
            this.elements.addAll(elements);
        }
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("Filter[");
        sb.append("enabled=").append(isEnabled);
        sb.append(", showAll=").append(showAll);
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

    public boolean isMatchAll() {
        return matchAll;
    }

    public boolean isShowAll() {
        return showAll;
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
                boolean elementResult =  element.evaluateSample(sampleName);
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

}
