/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.session;

import org.broad.igv.renderer.ColorScale;
import org.broad.igv.renderer.ColorScaleFactory;
import org.broad.igv.ui.color.ColorUtilities;

import javax.xml.bind.annotation.adapters.XmlAdapter;

/**
 * User: jacob
 * Date: 2013-Jan-02
 */
public class SessionXmlAdapters{

    public static class Color extends XmlAdapter<String, java.awt.Color>{
        @Override
        public String marshal(java.awt.Color value) throws Exception {
            if (value != null) {
                return ColorUtilities.colorToString(value);
            }
            return null;
        }

        @Override
        public java.awt.Color unmarshal(String colorString) throws Exception {
            if (colorString != null) {
                ColorUtilities.stringToColor(colorString);
            }
            return null;
        }
    }

    public static class Renderer extends XmlAdapter<String, org.broad.igv.renderer.Renderer>{

        @Override
        public String marshal(org.broad.igv.renderer.Renderer renderer) throws Exception {
            if (renderer != null) {
                RendererFactory.RendererType type = RendererFactory.getRenderType(renderer);
                if (type != null) {
                    return type.name();
                }
            }
            return null;
        }

        @Override
        public org.broad.igv.renderer.Renderer unmarshal(String rendererType) throws Exception {
            if (rendererType != null) {
                Class rendererClass = RendererFactory.getRendererClass(rendererType);
                if (rendererClass != null) {
                    return (org.broad.igv.renderer.Renderer) rendererClass.newInstance();
                }
            }
            return null;
        }
    }

    public static class ContinuousColorScale extends XmlAdapter<String, org.broad.igv.renderer.ContinuousColorScale>{

        @Override
        public String marshal(org.broad.igv.renderer.ContinuousColorScale colorScale) throws Exception {
            if (colorScale != null && !colorScale.isDefault()) {
                return colorScale.asString();
            }
            return null;
        }

        @Override
        public org.broad.igv.renderer.ContinuousColorScale unmarshal(String colorScaleString) throws Exception {
            if (colorScaleString != null) {
                ColorScale cs = ColorScaleFactory.getScaleFromString(colorScaleString);
                if(cs instanceof org.broad.igv.renderer.ContinuousColorScale){
                    return (org.broad.igv.renderer.ContinuousColorScale) cs;
                }
            }
            return null;
        }
    }

    public static class Height extends XmlAdapter<String, Integer>{

        @Override
        public String marshal(Integer height) throws Exception {
            if(height > 0){
                return "" + height;
            }
            return null;
        }

        @Override
        public Integer unmarshal(String heightString) throws Exception {
            if(heightString != null){
                return Integer.parseInt(heightString);
            }
            return null;
        }
    }
}
