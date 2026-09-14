package org.igv.renderer;

import org.igv.ui.color.ColorUtilities;

import java.awt.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 * A color scale defined by color stops -- values with a color -- blended on a log scale between stops.  For values
 * that span orders of magnitude, such as allele frequencies, where equal ratios should look like equal steps.
 * Values at or below the lowest stop, including zero and negative ("missing") values, take the lowest stop's color;
 * values at or above the highest stop take the highest stop's color.  Stop values must be positive.
 */
public class ColorStopScale extends AbstractColorScale {

    public static final String serializedClassName = "ColorStopScale";

    public record Stop(double value, Color color) {
    }

    private final List<Stop> stops;

    public ColorStopScale(List<Stop> stops) {
        if (stops.isEmpty()) {
            throw new IllegalArgumentException("A color stop scale needs at least one stop");
        }
        for (Stop stop : stops) {
            if (!(stop.value() > 0)) {
                throw new IllegalArgumentException("Color stop values must be positive: " + stop.value());
            }
        }
        List<Stop> sorted = new ArrayList<>(stops);
        sorted.sort(Comparator.comparingDouble(Stop::value));
        this.stops = Collections.unmodifiableList(sorted);
    }

    /**
     * Construct from the string form written by {@link #asString()}: {@code ColorStopScale;value r,g,b;value r,g,b;...}
     */
    public ColorStopScale(String string) {
        this(parse(string));
    }

    private static List<Stop> parse(String string) {
        String[] tokens = string.split(";");
        if (!tokens[0].trim().equals(serializedClassName)) {
            throw new IllegalArgumentException("Illegal ColorStopScale: " + string);
        }
        List<Stop> stops = new ArrayList<>();
        for (int i = 1; i < tokens.length; i++) {
            String[] parts = tokens[i].trim().split("\\s+");
            if (parts.length != 2) {
                throw new IllegalArgumentException("Illegal color stop: " + tokens[i]);
            }
            stops.add(new Stop(Double.parseDouble(parts[0]), ColorUtilities.stringToColor(parts[1])));
        }
        return stops;
    }

    /**
     * @return the stops, in increasing order of value
     */
    public List<Stop> getStops() {
        return stops;
    }

    @Override
    public Color getColor(float value) {
        return getColor((double) value);
    }

    public Color getColor(double value) {

        Stop lowest = stops.get(0);
        if (!(value > lowest.value())) {        // Also NaN
            return lowest.color();
        }
        Stop highest = stops.get(stops.size() - 1);
        if (value >= highest.value()) {
            return highest.color();
        }

        int i = 1;
        while (stops.get(i).value() < value) {
            i++;
        }
        Stop low = stops.get(i - 1);            // low.value < value <= high.value
        Stop high = stops.get(i);
        double t = (Math.log(value) - Math.log(low.value())) / (Math.log(high.value()) - Math.log(low.value()));
        return blend(low.color(), high.color(), t);
    }

    private static Color blend(Color c1, Color c2, double t) {
        return new Color(
                (int) Math.round(c1.getRed() + t * (c2.getRed() - c1.getRed())),
                (int) Math.round(c1.getGreen() + t * (c2.getGreen() - c1.getGreen())),
                (int) Math.round(c1.getBlue() + t * (c2.getBlue() - c1.getBlue())),
                (int) Math.round(c1.getAlpha() + t * (c2.getAlpha() - c1.getAlpha())));
    }

    @Override
    public String asString() {
        StringBuilder buf = new StringBuilder(serializedClassName);
        for (Stop stop : stops) {
            buf.append(';').append(stop.value()).append(' ').append(ColorUtilities.colorToString(stop.color()));
        }
        return buf.toString();
    }

    @Override
    public boolean isDefault() {
        return false;
    }
}
