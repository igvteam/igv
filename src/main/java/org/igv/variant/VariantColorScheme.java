package org.igv.variant;

import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.renderer.AbstractColorScale;
import org.igv.renderer.ContinuousColorScale;
import org.igv.renderer.MonocolorScale;
import org.igv.ui.color.ColorUtilities;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.PrintWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Set;

/**
 * A named set of colors for VCF INFO attribute values, read from a tab delimited file.  The format is the "#colors"
 * section of an IGV sample information file, so the same syntax works in both places.
 *
 * <pre>
 * #name=ClinVar significance
 * #description=Benign (blue) through pathogenic (red)
 * #colors
 * CLNSIG   Pathogenic  202,0,32
 * CLNSIG   *           150,150,150
 * CADD_PHRED   0:40    255,255,200     255,0,0
 * </pre>
 *
 * A row is "INFO key", "value", then one or two colors.  A value of "*" sets the color for values the scheme does
 * not otherwise cover.  A value of the form "min:max" defines a continuous scale over a numeric attribute --
 * one color shades from white to that color, two colors shade from the first to the second.
 */
public class VariantColorScheme {

    private static Logger log = LogManager.getLogger(VariantColorScheme.class);

    private static final String NAME_DIRECTIVE = "#name=";
    private static final String DESCRIPTION_DIRECTIVE = "#description=";
    private static final String SOURCE_DIRECTIVE = "#source=";
    private static final String COLORS_SECTION = "#colors";

    /**
     * The value that sets the color for anything the scheme does not explicitly cover.
     */
    static final String WILDCARD = "*";

    /**
     * Written in place of a value to record that a numeric attribute holds categories, not quantities, so it is
     * colored by value rather than by a scale.  IGV cannot tell the two apart, so the answer is kept here.
     */
    static final String CATEGORICAL = "categorical";

    private String name;
    private String description;
    private String source;

    /**
     * The file this scheme was read from, in the IGV scheme directory.  Null for schemes shipped with IGV.
     */
    private File file;

    /**
     * INFO key (upper case) -> attribute value -> color.  Values keep the case they were written in, so a file
     * IGV saves reads back the way the user wrote it; lookup is case insensitive via {@link #valueIndex}.
     */
    private final Map<String, Map<String, Color>> colors = new LinkedHashMap<>();

    /**
     * INFO key (upper case) -> lower case value -> the value as written, for case insensitive lookup.
     */
    private final Map<String, Map<String, String>> valueIndex = new LinkedHashMap<>();

    /**
     * INFO key (upper case) -> color for values not otherwise covered.
     */
    private final Map<String, Color> defaultColors = new HashMap<>();

    /**
     * INFO key (upper case) -> scale, for numeric attributes given as a "min:max" range.
     */
    private final Map<String, AbstractColorScale> scales = new LinkedHashMap<>();

    /**
     * Numeric INFO keys (upper case) declared to hold categories.
     */
    private final Set<String> categoricalKeys = new LinkedHashSet<>();

    VariantColorScheme(String name) {
        this.name = name;
    }

    /**
     * Parse a color scheme.  Unparseable rows are logged and skipped, a scheme is not rejected outright for one
     * bad line.
     *
     * @param reader      the scheme contents
     * @param defaultName name to use if the file has no "#name=" directive
     */
    public static VariantColorScheme parse(BufferedReader reader, String defaultName) throws IOException {

        VariantColorScheme scheme = new VariantColorScheme(defaultName);
        String line;

        while ((line = reader.readLine()) != null) {

            line = line.trim();
            if (line.isEmpty()) {
                continue;
            }

            if (line.startsWith("#")) {
                String lower = line.toLowerCase();
                if (lower.startsWith(NAME_DIRECTIVE)) {
                    scheme.name = line.substring(NAME_DIRECTIVE.length()).trim();
                } else if (lower.startsWith(DESCRIPTION_DIRECTIVE)) {
                    scheme.description = line.substring(DESCRIPTION_DIRECTIVE.length()).trim();
                } else if (lower.startsWith(SOURCE_DIRECTIVE)) {
                    scheme.source = line.substring(SOURCE_DIRECTIVE.length()).trim();
                }
                // Anything else, including the "#colors" section marker, is a comment
                continue;
            }

            String[] tokens = line.split("\t");

            if (tokens.length == 2 && CATEGORICAL.equalsIgnoreCase(tokens[1].trim())) {
                scheme.categoricalKeys.add(tokens[0].trim().toUpperCase());
                continue;
            }

            if (tokens.length < 3) {
                log.warn("Skipping color scheme line, expected at least 3 tab delimited fields: " + line);
                continue;
            }
            scheme.addRow(tokens);
        }

        return scheme;
    }

    private void addRow(String[] tokens) {

        String key = tokens[0].trim().toUpperCase();
        String value = tokens[1].trim();
        Color color = ColorUtilities.stringToColor(tokens[2].trim(), null);
        if (color == null) {
            log.warn("Skipping color scheme row with unparseable color: " + String.join("\t", tokens));
            return;
        }

        if (value.contains(":")) {
            String[] range = value.split(":");
            try {
                if (range.length > 2) {
                    // "min:mid:max" with three colors -- a gradient through a midpoint
                    double min = Double.parseDouble(range[0].trim());
                    double mid = Double.parseDouble(range[1].trim());
                    double max = Double.parseDouble(range[2].trim());
                    Color midColor = tokens.length > 3 ? ColorUtilities.stringToColor(tokens[3].trim(), null) : null;
                    Color maxColor = tokens.length > 4 ? ColorUtilities.stringToColor(tokens[4].trim(), null) : null;
                    if (midColor == null || maxColor == null) {
                        log.warn("Skipping color scheme row, a min:mid:max range needs three colors: "
                                + String.join("\t", tokens));
                        return;
                    }
                    scales.put(key, new ContinuousColorScale(min, mid, max, color, midColor, maxColor));
                    return;
                }

                float min = Float.parseFloat(range[0].trim());
                float max = Float.parseFloat(range[1].trim());
                if (tokens.length > 3) {
                    Color maxColor = ColorUtilities.stringToColor(tokens[3].trim(), null);
                    if (maxColor != null) {
                        scales.put(key, new ContinuousColorScale(min, max, color, maxColor));
                        return;
                    }
                }
                scales.put(key, new MonocolorScale(min, max, color));
            } catch (NumberFormatException e) {
                log.warn("Skipping color scheme row with unparseable range: " + String.join("\t", tokens));
            }
        } else if (WILDCARD.equals(value)) {
            defaultColors.put(key, color);
        } else {
            putColor(key, value, color);
        }
    }

    private void putColor(String key, String value, Color color) {
        // Reuse the spelling already recorded, if any, so repeated edits do not create near duplicate rows
        String existing = valueIndex.computeIfAbsent(key, k -> new LinkedHashMap<>())
                .putIfAbsent(value.toLowerCase(), value);
        colors.computeIfAbsent(key, k -> new LinkedHashMap<>()).put(existing == null ? value : existing, color);
    }

    /**
     * Return the color this scheme assigns to an attribute value, or null if it does not cover it.  Matching is
     * case insensitive.
     */
    public Color getColor(String infoKey, String value) {
        if (infoKey == null || value == null) {
            return null;
        }
        String key = infoKey.toUpperCase();
        Map<String, String> index = valueIndex.get(key);
        String stored = index == null ? null : index.get(value.toLowerCase());
        Map<String, Color> valueColors = colors.get(key);
        Color color = stored == null || valueColors == null ? null : valueColors.get(stored);
        return color != null ? color : defaultColors.get(key);
    }

    /**
     * Return the continuous scale for a numeric attribute, or null if this scheme does not define one.
     */
    public AbstractColorScale getScale(String infoKey) {
        return infoKey == null ? null : scales.get(infoKey.toUpperCase());
    }

    /**
     * @return the INFO keys this scheme assigns colors to.
     */
    public Set<String> getKeys() {
        Set<String> keys = new LinkedHashSet<>(colors.keySet());
        keys.addAll(defaultColors.keySet());
        keys.addAll(scales.keySet());
        keys.addAll(categoricalKeys);
        return keys;
    }

    /**
     * @return true if this scheme declares a numeric attribute to hold categories rather than quantities.
     */
    public boolean isCategorical(String infoKey) {
        return infoKey != null && categoricalKeys.contains(infoKey.toUpperCase());
    }

    /**
     * @return the value -> color assignments for an INFO key, for display in a legend.  Does not include the
     * wildcard default or continuous scales.
     */
    public Map<String, Color> getColors(String infoKey) {
        Map<String, Color> valueColors = infoKey == null ? null : colors.get(infoKey.toUpperCase());
        return valueColors == null ? Collections.emptyMap() : Collections.unmodifiableMap(valueColors);
    }

    /**
     * @return the color for values this scheme does not list, or null if it sets none.
     */
    public Color getDefaultColor(String infoKey) {
        return infoKey == null ? null : defaultColors.get(infoKey.toUpperCase());
    }

    /**
     * Set the color for one value of an INFO attribute.
     */
    public void setColor(String infoKey, String value, Color color) {
        String key = infoKey.toUpperCase();
        if (WILDCARD.equals(value)) {
            defaultColors.put(key, color);
        } else {
            putColor(key, value, color);
        }
    }

    /**
     * Remove the color for one value of an INFO attribute.  Values it no longer lists take their color from the
     * palette, as they did before the scheme covered them.
     */
    public void removeColor(String infoKey, String value) {

        String key = infoKey.toUpperCase();

        if (WILDCARD.equals(value)) {
            defaultColors.remove(key);
            return;
        }

        Map<String, String> index = valueIndex.get(key);
        String stored = index == null ? null : index.remove(value.toLowerCase());
        Map<String, Color> valueColors = colors.get(key);
        if (stored != null && valueColors != null) {
            valueColors.remove(stored);
        }
    }

    /**
     * Return an independent copy, so an editor can discard its changes.
     */
    public VariantColorScheme copy() {

        VariantColorScheme copy = new VariantColorScheme(name);
        copy.description = description;
        copy.source = source;
        copy.file = file;
        copy.defaultColors.putAll(defaultColors);
        copy.scales.putAll(scales);
        copy.categoricalKeys.addAll(categoricalKeys);
        for (Map.Entry<String, Map<String, Color>> entry : colors.entrySet()) {
            copy.colors.put(entry.getKey(), new LinkedHashMap<>(entry.getValue()));
        }
        for (Map.Entry<String, Map<String, String>> entry : valueIndex.entrySet()) {
            copy.valueIndex.put(entry.getKey(), new LinkedHashMap<>(entry.getValue()));
        }
        return copy;
    }

    /**
     * Set the color scale for a numeric INFO attribute.
     */
    public void setScale(String infoKey, AbstractColorScale scale) {
        scales.put(infoKey.toUpperCase(), scale);
    }

    /**
     * Declare that a numeric INFO attribute holds categories rather than quantities.
     */
    public void setCategorical(String infoKey) {
        categoricalKeys.add(infoKey.toUpperCase());
    }

    /**
     * Write this scheme in the tab delimited form it is read from.  This is the only place the format is
     * written, so hand edited and IGV written files stay interchangeable.
     */
    public void write(PrintWriter writer) {

        writer.println(NAME_DIRECTIVE + name);
        if (description != null) {
            writer.println(DESCRIPTION_DIRECTIVE + description);
        }
        if (source != null) {
            writer.println(SOURCE_DIRECTIVE + source);
        }
        writer.println(COLORS_SECTION);

        for (String key : getKeys()) {
            if (categoricalKeys.contains(key)) {
                writer.println(key + "\t" + CATEGORICAL);
            }
            AbstractColorScale scale = scales.get(key);
            if (scale instanceof ContinuousColorScale) {
                writeScale(writer, key, (ContinuousColorScale) scale);
            }
            for (Map.Entry<String, Color> entry : colors.getOrDefault(key, Collections.emptyMap()).entrySet()) {
                writer.println(key + "\t" + entry.getKey() + "\t" + ColorUtilities.colorToString(entry.getValue()));
            }
            Color defaultColor = defaultColors.get(key);
            if (defaultColor != null) {
                writer.println(key + "\t" + WILDCARD + "\t" + ColorUtilities.colorToString(defaultColor));
            }
        }
    }

    private static void writeScale(PrintWriter writer, String key, ContinuousColorScale scale) {
        if (scale.isUseDoubleGradient()) {
            writer.println(key
                    + "\t" + scale.getMinimum() + ":" + scale.getNegStart() + ":" + scale.getMaximum()
                    + "\t" + ColorUtilities.colorToString(scale.getMinColor())
                    + "\t" + ColorUtilities.colorToString(scale.getMidColor())
                    + "\t" + ColorUtilities.colorToString(scale.getMaxColor()));
        } else {
            writer.println(key
                    + "\t" + scale.getMinimum() + ":" + scale.getMaximum()
                    + "\t" + ColorUtilities.colorToString(scale.getMinColor())
                    + "\t" + ColorUtilities.colorToString(scale.getMaxColor()));
        }
    }

    void setName(String name) {
        this.name = name;
    }

    public String getName() {
        return name;
    }

    public String getDescription() {
        return description;
    }

    public String getSource() {
        return source;
    }

    void setSource(String source) {
        this.source = source;
    }

    /**
     * @return the file this scheme was read from, or null if it is shipped with IGV.
     */
    public File getFile() {
        return file;
    }

    void setFile(File file) {
        this.file = file;
    }

    /**
     * @return true if this scheme ships with IGV, and so cannot be removed or edited in place.
     */
    public boolean isBuiltIn() {
        return file == null;
    }

    @Override
    public String toString() {
        return name;
    }
}
