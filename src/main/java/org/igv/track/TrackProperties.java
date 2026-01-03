package org.igv.track;

//~--- JDK imports ------------------------------------------------------------

import org.igv.feature.genome.load.TrackConfig;
import org.igv.logging.*;
import org.igv.ui.color.ColorUtilities;
import org.igv.util.ParsingUtils;

import java.awt.*;
import java.util.LinkedHashMap;
import java.util.Map;

import static org.igv.util.StringUtils.stripQuotes;


/**
 * This class is based on the UCSC track configuration parameters
 *
 * @author jrobinso
 */
public class TrackProperties {

    private static Logger log = LogManager.getLogger(TrackProperties.class);


    public enum BaseCoord {
        ZERO, ONE, UNSPECIFIED
    }

    /**
     * The original track line from a track file or track hub.  This is not always set.  Not sure why we need this.
     */
    private String trackLine;


    private Track.DisplayMode displayMode;

    /**
     * Base coordinate system,  either 0 or 1
     */
    private BaseCoord baseCoord = BaseCoord.UNSPECIFIED;


    private String type;

    /**
     * The track name.  Will be displayed to the left of the track.
     */
    private String name;

    /**
     * The track description.  The description, when supplied, will appear in popup text
     * when hovering over the track.  May optionally be displayed centered over the track
     * on certain chart types (currently scatter plot and barchart).
     */
    private String description;

    /**
     * A url to an web page associated with this track.  This is currently not used.
     */
    private String url;

    /**
     * The track height in pixels
     */
    private Integer height;

    private Integer minHeight;

    private boolean gffTags = false;

    /**
     * The default color for the track.  This can be overridden by individual feature lines in
     * certain formats, notably BED.
     */
    private Color color;

    /**
     * An alternate color. The use of this color depends on the particular graph type (renderer)
     * being displayed.  See documentation for specific renderes for details.
     */
    private Color altColor;

    private Color midColor;

    private String genome;

    private Boolean autoscale;

    private Float minValue;

    private Float maxValue;

    private Float midValue;

    private Float neutralFromValue;

    private Float neutralToValue ;

    private boolean drawYLine = false;

    private Class rendererClass;

    private WindowFunction windowingFunction;

    private int smoothingWindow;

    private Boolean itemRGB = null;

    private Boolean useScore = null;

    private int featureVisibilityWindow = -1;

    private boolean logScale = false;

    private Float yLine;

    private boolean sortable = true;

    private boolean alternateExonColor = false;

    private String dataURL;

    private String indexURL;

    private String coverageURL;

    /**
     * Non-standard track field to indicate file format
     */
    private String format;

    /**
     * Track attributes (meta data)
     */
    private Map<String, String> attributes;


    public TrackProperties() {

    }


    /**
     * Initialize TrackProperties from a TrackConfig object.  TrackConfig is a json-friendly bridge object
     * that is used to load track properties from an igv.js json file or a track hub.
     *
     * @param trackConfig
     */
    public TrackProperties(TrackConfig trackConfig) {
        if(trackConfig.color != null) {
            this.color = parseColor(trackConfig.color);
        }
        if(trackConfig.altColor != null) {
            this.altColor = parseColor(trackConfig.altColor);
        }
        if(trackConfig.displayMode != null) {
            this.displayMode = parseDisplayMode(trackConfig.displayMode);
        }
        setFeatureVisibilityWindow(trackConfig.visibilityWindow != null ? trackConfig.visibilityWindow : -1);
        if(trackConfig.min != null) {
            setMinValue(trackConfig.min);
        }
        if(trackConfig.max != null) {
            setMaxValue(trackConfig.max);
        }
        if(trackConfig.autoscale != null) {
            setAutoScale(trackConfig.autoscale);
        }
        if(trackConfig.height != null) {
            setHeight(trackConfig.height);
        }
        if(trackConfig.minHeight != null) {
            setMinHeight(trackConfig.minHeight);
        }
    }


    private Color parseColor(String color) {
        if (color != null) {
            try {
                return ColorUtilities.stringToColor(color);
            } catch (Exception e) {
                log.error("Error parsing color string: " + color, e);
            }
        }
        return null;
    }

    private Track.DisplayMode parseDisplayMode(String displayMode) {
        if (displayMode != null) {
            try {
                return Track.DisplayMode.valueOf(stripQuotes(displayMode));
            } catch (Exception e) {
                log.error("Error parsing displayMode " + displayMode, e);
            }
        }
        return null;
    }

    public void setType(String type) {
        this.type = type;
    }

    public String getType() {
        return type;
    }

    public String getFormat() {
        return format;
    }

    public void setFormat(String format) {
        this.format = format;
    }

    public void setTrackLine(String trackLine) {
        this.trackLine = trackLine;
    }

    public String getTrackLine() {
        return this.trackLine;
    }

    public boolean isSortable() {
        return sortable;
    }

    public void setSortable(boolean sortable) {
        this.sortable = sortable;
    }

    public boolean isLogScale() {
        return logScale;
    }

    public void setLogScale(boolean logScale) {
        this.logScale = logScale;
    }

    public int getFeatureVisibilityWindow() {
        return featureVisibilityWindow;
    }

    public void setFeatureVisibilityWindow(int featureVisibilityWindow) {
        this.featureVisibilityWindow = featureVisibilityWindow;
    }

    public Boolean isUseScore() {
        return useScore;
    }

    public void setUseScore(Boolean useScore) {
        this.useScore = useScore;
    }

    public Boolean isItemRGB() {
        return itemRGB;
    }

    public void setItemRGB(Boolean itemRGB) {
        this.itemRGB = itemRGB;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }


    public String getDescription() {
        return description;
    }


    public void setDescription(String description) {
        this.description = description;
    }


    public String getUrl() {
        return url;
    }


    public void setUrl(String url) {
        this.url = url;
    }


    public Integer getHeight() {
        return height;
    }


    public void setHeight(Integer height) {
        this.height = height;
    }


    public Color getColor() {
        return color;
    }


    public void setColor(Color color) {
        this.color = color;
    }

    public Color getAltColor() {
        return altColor;
    }


    public void setAltColor(Color altColor) {
        this.altColor = altColor;
    }

    public String getGenome() {
        return genome;
    }

    public void setGenome(String genome) {
        this.genome = genome;
    }

    public Boolean getAutoScale() {
        return this.autoscale;
    }

    public void setAutoScale(Boolean autoScale) {
        this.autoscale = autoScale;
    }
    public Float getMinValue() {
        return minValue;
    }

    public void setMinValue(float minValue) {
        this.minValue = minValue;
    }

    public Float getMaxValue() {
        return maxValue;
    }

    public void setMaxValue(Float maxValue) {
        this.maxValue = maxValue;
    }

    public Float getMidValue() {
        return midValue;
    }

    public void setMidValue(Float midValue) {
        this.midValue = midValue;
    }

    public WindowFunction getWindowingFunction() {
        return windowingFunction;
    }


    public void setWindowingFunction(WindowFunction windowingFunction) {
        this.windowingFunction = windowingFunction;
    }


    public int getSmoothingWindow() {
        return smoothingWindow;
    }


    public void setSmoothingWindow(int smoothingWindow) {
        this.smoothingWindow = smoothingWindow;
    }


    public Class getRendererClass() {
        return rendererClass;
    }

    public void setRendererClass(Class rendererClass) {
        this.rendererClass = rendererClass;
    }


    public Color getMidColor() {
        return midColor;
    }

    public void setMidColor(Color midColor) {
        this.midColor = midColor;
    }

    public boolean isDrawYLine() {
        return drawYLine;
    }

    public void setDrawYLine(boolean drawYLine) {
        this.drawYLine = drawYLine;
    }

    public Integer getMinHeight() {
        return minHeight;
    }

    public void setMinHeight(Integer minHeight) {
        this.minHeight = minHeight;
    }

    public BaseCoord getBaseCoord() {
        return baseCoord;
    }

    public void setBaseCoord(BaseCoord baseCoord) {
        this.baseCoord = baseCoord;
    }

    public Float getNeutralFromValue() {
        return neutralFromValue;
    }
    public void setNeutralFromValue(Float neutralFromValue) {
        this.neutralFromValue = neutralFromValue;
    }

    public Float getNeutralToValue() {
        return neutralToValue;
    }

    public void setNeutralToValue(Float neutralToValue) {
        this.neutralToValue = neutralToValue;
    }

    public Float getyLine() {
        return yLine;
    }

    public void setyLine(Float yLine) {
        this.yLine = yLine;
    }

    public void setGffTags(boolean gffTags) {
        this.gffTags = gffTags;
    }

    public boolean isGffTags() {
        return gffTags;
    }

    public boolean isAlternateExonColor() {
        return alternateExonColor;
    }

    public void setAlternateExonColor(boolean alternateExonColor) {
        this.alternateExonColor = alternateExonColor;
    }

    public void setDisplayMode(Track.DisplayMode displayMode) {
        this.displayMode = displayMode;
    }

    public Track.DisplayMode getDisplayMode() {
        return displayMode;
    }

    public String getDataURL() {
        return dataURL;
    }

    public String getIndexURL() {
        return indexURL;
    }

    public String getCoverageURL() {
        return coverageURL;
    }

    public void setDataURL(String dataURL) {
        this.dataURL = dataURL;
    }

    public void setIndexURL(String indexURL) {
        this.indexURL = indexURL;
    }

    public void setCoverageURL(String coverageURL) {
        this.coverageURL = coverageURL;
    }


    public Map<String, String> getAttributes() {
        return attributes;
    }


    /**
     * example:
     * Library=DNA_Lib 1387;Sample=NH-Osteoblast;Antibody=H3K4me3
     *
     * @param value
     */
    public void setMetaData(String value) {

        if (attributes == null) {
            attributes = new LinkedHashMap<String, String>();   // <= maintain order
        }

        String[] attStrings = ParsingUtils.SEMI_COLON_PATTERN.split(value);
        for (String att : attStrings) {
            String[] kv = ParsingUtils.EQ_PATTERN.split(att, 2);
            if (kv.length == 2) {
                attributes.put(kv[0], kv[1]);
            } else {
                log.warn("Skipping meta value: " + value + ".  Missing '=' token?");
            }
        }

    }


}
