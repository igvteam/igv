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

package org.broad.igv.session;

import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.renderer.ColorScale;
import org.broad.igv.renderer.ColorScaleFactory;
import org.broad.igv.track.DataSourceTrack;
import org.broad.igv.track.DataTrack;
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
                return ColorUtilities.stringToColor(colorString);
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

    public static class Genome extends XmlAdapter<String, org.broad.igv.feature.genome.Genome>{

        @Override
        public String marshal(org.broad.igv.feature.genome.Genome v) throws Exception {
            return v.getId();
        }

        @Override
        public org.broad.igv.feature.genome.Genome unmarshal(String v) throws Exception {
            org.broad.igv.feature.genome.Genome genome = GenomeManager.getInstance().getCurrentGenome();
            if(genome != null && !genome.getId().equals(v)){
                throw new IllegalStateException("Must load the proper genome before unmarshalling");
            }
            return genome;
        }

    }

    public static class DataTrackIDAdapter extends XmlAdapter<String, org.broad.igv.track.DataTrack> {

        @Override
        public String marshal(DataTrack dataTrack) throws Exception {
            return dataTrack.getId();
        }

        @Override
        public DataTrack unmarshal(String trackId) throws Exception {
            DataTrack dataTrack = (DataTrack) IGVSessionReader.getMatchingTrack(trackId, null);
            if(dataTrack == null){
                dataTrack = new DataSourceTrack(null, trackId, null, null);
            }
            return dataTrack;
        }
    }

}
