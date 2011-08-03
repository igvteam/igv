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

import org.broad.igv.util.IGVHttpClientUtils;
import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;
import org.json.JSONTokener;

import java.io.BufferedInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

/**
 * Utility class for accessing the GenomeSpace data manager web service
 * 
 * @author Jim Robinson
 * @date Aug 2, 2011
 */
public class GSDMUtils {

    /**
     * Fetch the contents of the GenomeSpace directory.
     *
     * @param directoryURL
     * @return
     * @throws IOException
     * @throws JSONException
     */
    public static GSDirectoryListing getDirectoryListing(URL directoryURL) throws IOException, JSONException {
        StringBuffer buf = new StringBuffer();
        InputStream is = null;
        try {
            is = IGVHttpClientUtils.openConnectionStream(directoryURL);
            BufferedInputStream bis = new BufferedInputStream(is);
            int b;
            while ((b = bis.read()) >= 0) {
                buf.append((char) b);
            }

            JSONTokener tk = new JSONTokener(buf.toString());
            JSONObject obj = new JSONObject(tk);

            JSONObject directory = (JSONObject) obj.get("directory");
            String dirUrlString = directory.get("url").toString();

            LinkedList<GSFileMetadata> elements = new LinkedList();
            if (obj.has("contents")) {
                Object c = obj.get("contents");
                List<JSONObject> contents = new ArrayList();
                if (c instanceof JSONObject) {
                    contents.add((JSONObject) c);
                } else {
                    JSONArray tmp = (JSONArray) c;
                    int l = tmp.length();
                    for (int i = 0; i < l; i++) {
                        contents.add((JSONObject) tmp.get(i));
                    }
                }

                ArrayList<GSFileMetadata> dirElements = new ArrayList();
                ArrayList<GSFileMetadata> fileElements = new ArrayList();
                int contentsLength = contents.size();
                for (int i = 0; i < contentsLength; i++) {
                    JSONObject o = contents.get(i);
                    String name = (String) o.get("name");
                    String objurl = (String) o.get("url");
                    if (o.get("directory").equals("true")) {
                        dirElements.add(new GSFileMetadata(name, objurl, "", "", true));
                    } else {
                        JSONObject dataFormat = o.has("dataFormat") ? (JSONObject) o.get("dataFormat") : null;
                        String format = dataFormat == null ? "" : dataFormat.getString("name");
                        String size = (String) o.get("size");
                        fileElements.add(new GSFileMetadata(name, objurl, format, size, false));
                    }
                }

                elements.addAll(dirElements);
                elements.addAll(fileElements);
            }

            return new GSDirectoryListing(dirUrlString, elements);

        } finally {
            is.close();
        }
    }

}
