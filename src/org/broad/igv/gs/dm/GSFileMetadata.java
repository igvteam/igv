/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package org.broad.igv.gs.dm;

import org.json.JSONException;
import org.json.JSONObject;

/**
 * Represents a file or directory in GS storage.
 *
 * @author Jim Robinson
 * @date Aug 2, 2011
 */
public class GSFileMetadata {
    private boolean isDirectory;
    private String name;
    private String path;
    private String url;
    private String format;
    private String size;

    public GSFileMetadata(String name, String path, String url, String format, String size, boolean isDirectory) {
        this.isDirectory = isDirectory;
        this.name = name;
        this.path = path;
        this.url = url;
        this.format = format;
        this.size = size;
    }

    public GSFileMetadata(JSONObject o) throws JSONException {
        name = (String) o.get("name");
        path = (String) o.get("path");
        url = (String) o.get("url");
        isDirectory = (Boolean) o.get("isDirectory");
        if (o.has("dataFormat")) {
            JSONObject dataFormat = o.has("dataFormat") ? (JSONObject) o.get("dataFormat") : null;
            format = dataFormat == null ? "" : dataFormat.getString("name");
            size =  o.get("size").toString();
        }

    }

    public String toString() {
        return getName();
    }

    public boolean isDirectory() {
        return isDirectory;
    }

    public String getName() {
        return name;
    }

    public String getPath() {
        return path;
    }


    public String getUrl() {
        return url;
    }

    public String getFormat() {
        return format;
    }

    public String getSize() {
        return size;
    }
}
