package org.broad.igv.session;

import org.broad.igv.track.BlatTrack;


import java.util.HashMap;
import java.util.Map;

/**
 * Session Element types
 */
public class SessionElement {


    public static final String PANEL = "Panel";
    public static final String PANEL_LAYOUT = "PanelLayout";
    public static final String TRACK = "Track";
    public static final String COLOR_SCALE = "ColorScale";
    public static final String COLOR_SCALES = "ColorScales";
    public static final String DATA_TRACK = "DataTrack";
    public static final String DATA_TRACKS = "DataTracks";
    public static final String FEATURE_TRACKS = "FeatureTracks";
    public static final String DATA_FILE = "DataFile";
    public static final String RESOURCE = "Resource";
    public static final String RESOURCES = "Resources";
    public static final String FILES = "Files";
    public static final String FILTER_ELEMENT = "FilterElement";
    public static final String FILTER = "Filter";
    public static final String SESSION = "Session";
    public static final String GLOBAL = "Global";
    public static final String REGION = "Region";
    public static final String REGIONS = "Regions";
    public static final String DATA_RANGE = "DataRange";
    public static final String PREFERENCES = "Preferences";
    public static final String PROPERTY = "Property";
    public static final String GENE_LIST = "GeneList";
    public static final String HIDDEN_ATTRIBUTES = "HiddenAttributes";
    public static final String VISIBLE_ATTRIBUTES = "VisibleAttributes";
    public static final String ATTRIBUTE = "Attribute";
    public static final String VISIBLE_ATTRIBUTE = "VisibleAttribute";
    public static final String FRAME = "Frame";


    /**
     * Map from an XML clazz element -> Java class
     *
     * @param className
     * @return
     */
    public static Class getClass(String className) throws ClassNotFoundException {

        if (classMap.containsKey(className)) {
            return classMap.get(className);
        } else {
            return Class.forName(className);
        }
    }

    public static String getXMLClassName(Class clazz) {

      //  if (classNameMap.containsKey(clazz)) {
      //      return classNameMap.get(clazz);
      //  } else {
            return clazz.getName();     // <= for backward compatibility
      //  }

    }


    static Map<String, Class> classMap = new HashMap<>();
    static Map<Class, String> classNameMap = new HashMap<>();
    static {
        classMap.put("BlatTrack", org.broad.igv.track .BlatTrack.class);
        classMap.put("FeatureTrack", org.broad.igv.track .FeatureTrack.class);
        classMap.put("SequenceTrack", org.broad.igv.track .SequenceTrack.class);
        classMap.put("DataTrack", org.broad.igv.track .DataTrack.class);
        classMap.put("DataSourceTrack", org.broad.igv.track.DataSourceTrack.class);

        for(Map.Entry<String, Class> entry : classMap.entrySet()) {
            classNameMap.put(entry.getValue(), entry.getKey());
        }
    }

}
