package org.broad.igv.variant;

import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.IGVPreferences;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.renderer.ColorScale;
import org.broad.igv.renderer.ColorScaleFactory;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.ui.color.PaletteColorTable;

import java.awt.*;
import java.util.*;
import java.util.List;

public class AttributeColorManager {

    private static final Map<Key, ColorScale> attributes = new TreeMap<>();
    public static final Map<String, ColorScale> activeColorScalesMap = Collections.synchronizedMap(new HashMap<>());
    public static final String VARIANT_INFO_PREFERENCE_KEY = "VARIANT_INFO_";
    public static final String VARIANT_FORMAT_PREFERENCE_KEY = "VARIANT_FORMAT_";

    public static PaletteColorTable DEFAULT_BOOLEAN_COLORS = new PaletteColorTable(Color.GRAY);
    static {
        DEFAULT_BOOLEAN_COLORS.put("true", Color.BLACK);
        DEFAULT_BOOLEAN_COLORS.put("false", Color.RED);
    }

    /**
     * get the key name for the preferences file for a given color scale
     * @param type
     * @param id
     * @return
     */
    static String getPreferencesKey(Type type, String id) {
        return switch(type) {
            case INFO -> VARIANT_INFO_PREFERENCE_KEY;
            case FORMAT -> VARIANT_FORMAT_PREFERENCE_KEY;
        } + id;
    }

    public enum Type{
        INFO, FORMAT
    }



    public record Key(Type type, String name) implements Comparable<Key> {
        public Key {
            Objects.requireNonNull(type);
            Objects.requireNonNull(name);
        }
    
        final static Comparator<Key> COMPARATOR = Comparator.comparing(Key::type)
                .thenComparing(Key::name);

        @Override
        public int compareTo(Key o) {
            return COMPARATOR.compare(this, o);
        }
    };

    public static PaletteColorTable getBooleanColorTable(Type type, String id){
        return (PaletteColorTable) attributes.computeIfAbsent(new Key(type, id), k -> DEFAULT_BOOLEAN_COLORS);
    }

    private static String getDefaultScaleString(ColorScaleType type) {
        IGVPreferences preferences = PreferencesManager.getPreferences();
        return switch( type ){
            case CONTINUOUS -> preferences.get(Constants.DEFAULT_CONTINUOUS_COLOR_SCALE);
            case DISCRETE -> preferences.get(Constants.DEFAULT_DISCRETE_COLOR_SCALE);
            case FLAG -> preferences.get(Constants.DEFAULT_BOOLEAN_COLOR_SCALE);
        };
    }

    public enum ColorScaleType {
        CONTINUOUS, DISCRETE, FLAG;
    }
    
    /**
     *  Get the scale for the given 
     *  1. check for a saved scale
     *  2. check for a scale from preferences
     *  3. get the default scale for the value type
     *   
     * @param type
     * @param id
     * @param colorScaleType
     * @return
     */
    public static ColorScale getColorTable(Type type, String id, ColorScaleType colorScaleType){
        Objects.requireNonNull(type);
        Objects.requireNonNull(colorScaleType);
        Objects.requireNonNull(id);

        return attributes.computeIfAbsent(new Key(type, id), key -> {
            String preferenceString = PreferencesManager.getPreferences().get(
                    getPreferencesKey(key.type(), key.name()),
                    getDefaultScaleString(colorScaleType));
            return ColorScaleFactory.getScaleFromString(preferenceString);
        });
    }

    public static void putColorTable(Type type, String id, ColorScale colorScale) {
        attributes.put(new Key(type, id), colorScale);
    }

    public static void putColorTable(Key key, ColorScale colorScale){
        attributes.put(key, colorScale);
    }

    public static void remove(Key key) {
        attributes.remove(key);
    }

    public static Key preferenceIdToKey(String id){
        if(id.startsWith(VARIANT_INFO_PREFERENCE_KEY)){
            String name = id.substring(VARIANT_INFO_PREFERENCE_KEY.length());
            return new Key(Type.INFO, name);
        } else if(id.startsWith(VARIANT_FORMAT_PREFERENCE_KEY)){
            String name = id.substring(VARIANT_FORMAT_PREFERENCE_KEY.length());
            return new Key(Type.FORMAT, name);
        } else {
            throw new IllegalArgumentException("Can't convert attribute id to key: " + id);
        }
    }

    static {
        PaletteColorTable svtypeColors = new PaletteColorTable(ColorUtilities.getPalette("Pastel 1"));
        //from VCF 4.5 spec section 1.4.5
        List.of("DEL",
                "INS",
                "DUP",
                "INV",
                "CNV",
                "CNV:TR",
                "DUP:TANDEM",
                "DEL:ME",
                "INS:ME")
                .forEach(svtypeColors::get);
        attributes.put(new Key(Type.INFO, "SVTYPE"), svtypeColors);

        PaletteColorTable defaultBooleanColors = new PaletteColorTable(Color.GRAY);
        defaultBooleanColors.put("true", Color.BLACK);
        defaultBooleanColors.put("false", Color.RED);
    }


}
