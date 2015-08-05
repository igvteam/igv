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

package org.broad.igv.util.encode;

import java.io.File;
import java.util.Collection;
import java.util.Map;

/**
 * @author jrobinso
 *         Date: 10/31/13
 *         Time: 10:11 PM
 */
public class EncodeFileRecord {

    boolean selected = false;
    String path;
    Map<String, String> attributes;
    String trackName;

    public EncodeFileRecord(String path, Map<String, String> attributes) {
        this.path = path;
        this.attributes = attributes;
    }

    public String getPath() {
        return path;
    }

    public String getFileType() {
        //String off trailing gz, if present
        String filetype = path;
        if (filetype.endsWith(".gz")) {
            filetype = filetype.substring(0, filetype.length() - 3);
        }
        int idx = filetype.lastIndexOf(".");
        return filetype.substring(idx + 1);
    }

    public String getAttributeValue(String name) {
        String value = attributes.get(name);
        if (name.equals("type") && value == null) value = getFileType();
        return value;
    }

    public Collection<String> getAttributeNames() {
        return attributes.keySet();
    }

    public boolean containsText(String filter) {
        for (String value : attributes.values()) {
            if (value.contains(filter)) return true;
        }
        return false;
    }

    boolean isSelected() {
        return selected;
    }

    void setSelected(boolean selected) {
        this.selected = selected;
    }

    /**
     * Return a friendly name for the track.  Unfortunately it is neccessary to hardcode certain attributes.
     *
     * @return
     */
    public String getTrackName() {

        if (trackName == null) {
            StringBuffer sb = new StringBuffer();
            if(attributes.containsKey("cell")) sb.append(attributes.get("cell") + " ");
            if(attributes.containsKey("antibody")) sb.append(attributes.get("antibody") + " ");
            if(attributes.containsKey("dataType")) sb.append(attributes.get("dataType") + " ");
            if(attributes.containsKey("view")) sb.append(attributes.get("view") + " ");
            if(attributes.containsKey("replicate")) sb.append("rep " + attributes.get("replicate"));

            trackName = sb.toString().trim();
            if(sb.length() == 0) trackName = (new File(path)).getName();
        }

        return trackName;

    }

    /**
     * Test if record has a eough of meta-data to be interpretable
     *
     * @return
     */
    public boolean hasMetaData() {

        return  (attributes.containsKey("cell")) || (attributes.containsKey("antibody"));

    }
}
