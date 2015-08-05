/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */


package org.broad.igv.track;

//~--- JDK imports ------------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.util.ParsingUtils;

import java.awt.*;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;


/**
 * This class is based on the UCSC track configuration parameters
 *
 * @author jrobinso
 */
public class TrackProperties {

    private static Logger log = Logger.getLogger(TrackProperties.class);

    public enum BaseCoord {
        ZERO, ONE, UNSPECIFIED
    }

    private Track.DisplayMode displayMode;


    /**
     * Base coordinate system,  either 0 or 1
     */
    private BaseCoord baseCoord = BaseCoord.UNSPECIFIED;

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
    private int height;

    private int minHeight;

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

    private int offset;

    private boolean autoScaleFlag = false;

    private float minValue = Float.NaN;

    private float maxValue = Float.NaN;

    private float midValue = Float.NaN;

    private float neutralFromValue = Float.NaN;

    private float neutralToValue = Float.NaN;

    private boolean drawYLine = false;

    private Class rendererClass;

    private WindowFunction windowingFunction;

    private int smoothingWindow;

    private boolean itemRGB = true;

    private boolean useScore = false;

    private int featureVisibilityWindow = -1;

    private boolean logScale;

    private float yLine;

    private boolean sortable = true;

    private boolean alternateExonColor = false;

    private String dataURL;

    /**
     * Track attributes (meta data)
     */
    private Map<String, String> attributes;


    public TrackProperties() {

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

    public boolean isUseScore() {
        return useScore;
    }

    public void setUseScore(boolean useScore) {
        this.useScore = useScore;
    }

    public boolean isItemRGB() {
        return itemRGB;
    }

    public void setItemRGB(boolean itemRGB) {
        this.itemRGB = itemRGB;
    }

    /**
     * Method description
     *
     * @return
     */
    public int getOffset() {
        return offset;
    }

    public void setOffset(int offset) {
        this.offset = offset;
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


    public int getHeight() {
        return height;
    }


    public void setHeight(int height) {
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


    public boolean isAutoScale() {
        return autoScaleFlag || Float.isNaN(minValue) || Float.isNaN(maxValue);
    }

    public boolean getAutoScale()    {
        return this.autoScaleFlag;
    }

    public String getGenome() {
        return genome;
    }


    public void setGenome(String genome) {
        this.genome = genome;
    }


    public float getMinValue() {
        return minValue;
    }


    public void setMinValue(float minValue) {
        this.minValue = minValue;
    }


    public float getMaxValue() {
        return maxValue;
    }


    public void setMaxValue(float maxValue) {
        this.maxValue = maxValue;
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

    public void setAutoScale(boolean autoScale) {
        this.autoScaleFlag = autoScale;
    }

    public float getMidValue() {
        return midValue;
    }

    public void setMidValue(float midValue) {
        this.midValue = midValue;
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

    public int getMinHeight() {
        return minHeight;
    }

    public void setMinHeight(int minHeight) {
        this.minHeight = minHeight;
    }

    public BaseCoord getBaseCoord() {
        return baseCoord;
    }

    public void setBaseCoord(BaseCoord baseCoord) {
        this.baseCoord = baseCoord;
    }

    public float getNeutralFromValue() {
        return neutralFromValue;
    }

    public void setNeutralFromValue(float neutralFromValue) {
        this.neutralFromValue = neutralFromValue;
    }

    public float getNeutralToValue() {
        return neutralToValue;
    }

    public void setNeutralToValue(float neutralToValue) {
        this.neutralToValue = neutralToValue;
    }

    public float getyLine() {
        return yLine;
    }

    public void setyLine(float yLine) {
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

    public void setDataURL(String dataURL) {
        this.dataURL = dataURL;
    }


    public Map<String,String> getAttributes() {
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
                log.info("Skipping meta value: " + value + ".  Missing '=' token?");
            }
        }

    }


}
