/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.plugin.mongocollab;

import com.mongodb.ReflectionDBObject;
import org.broad.igv.feature.AbstractFeature;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.renderer.IGVFeatureRenderer;
import org.broad.igv.session.SubtlyImportant;
import org.broad.igv.ui.color.ColorUtilities;
import htsjdk.tribble.Feature;

import java.awt.*;

/**
 * Object mapping to Mongo database
 * ReflectionDBObject works with getters/setters, and
 * doesn't use the Java Beans case convention.
 * So (get/set)Chr maps to a field named "Chr", not "chr"
 * as we might prefer
 *
 * TODO Use existing feature interfaces/classes, which are long past
 * TODO overdue for refactoring
 */
public class DBFeature extends ReflectionDBObject implements Feature {

    static Color DEFAULT_COLOR = IGVFeatureRenderer.DULL_BLUE;

    private String chr;
    private int start;
    private int end;
    private String description;

    private Color color = DEFAULT_COLOR;
    private String name;

    @SubtlyImportant
    public DBFeature(){}

    DBFeature(String chr, int start, int end, String name, String description, Color color){
        this.chr = chr;
        this.start = start;
        this.end = end;
        this.name = name;
        this.description = description;
        this.color = color;
    }

    static DBFeature create(Feature feature){
        if(feature instanceof AbstractFeature){
            return create((AbstractFeature) feature);
        }
        return new DBFeature(feature.getChr(), feature.getStart(), feature.getEnd(), null, null, DEFAULT_COLOR);
    }

    static DBFeature create(AbstractFeature feature){
        Color color = feature.getColor() != null ? feature.getColor() : DEFAULT_COLOR;
        return new DBFeature(feature.getChr(), feature.getStart(), feature.getEnd(), feature.getName(), feature.getDescription(), color);
    }

    public String getChr() {
        return chr;
    }

    public void setChr(String chr) {
        this.chr = chr;
    }

    public String getDescription() {
        return description;
    }

    public void setDescription(String description) {
        this.description = description;
    }

    public int getEnd() {
        return end;
    }

    public void setEnd(int end) {
        this.end = end;
    }

    public int getStart() {
        return start;
    }

    public void setStart(int start) {
        this.start = start;
    }

    public void setName(String name){
        this.name = name;
    }

    /**
     * Get the name of this feature, in upper-case, or null if name is null
     * For case-insensitive searches
     * @return
     */
    public String getUpperName(){
        return this.name != null ? this.name.toUpperCase() : null;
    }

    /**
     * no-op, just here because ReflectionDB needs a setter
     * Does nothing
     * @param upperName
     */
    public void setUpperName(String upperName){
        //pass
    }

    public String getName(){
        return this.name;
    }

    /**
     * Setter for {@code color}, can take a {@link Color}
     * or a String. A String will be converted using
     * {@link ColorUtilities#stringToColor(String)}. We do this weird business
     * because MongoDB doesn't support {@link Color}s natively
     * @param color
     */
    public void setColor(Object color){
        if(color instanceof Color){
            this.color = (Color) color;
        }else if (color instanceof String){
            this.color = ColorUtilities.stringToColor((String) color);
        }
    }

    public Color getColor() {
        return color;
    }


    public IGVFeat createIGVFeature(){
        return new IGVFeat(this);
    }


    /**
     * This is the feature to be returned to {@code FeatureTrack}s.
     * Would really like to inherit both ReflectionDBObject and BasicFeature, but
     * that's not possible.
     */
    public static class IGVFeat extends BasicFeature{

        private DBFeature dbFeat;

        IGVFeat(DBFeature dbFeat){
            super(dbFeat.getChr(), dbFeat.getStart(), dbFeat.getEnd());
            this.dbFeat = dbFeat;
            setName(dbFeat.getName());
            setDescription(dbFeat.getDescription());
            setColor(dbFeat.getColor());
        }

        DBFeature getDBFeature(){
            return this.dbFeat;
        }
    }
}
