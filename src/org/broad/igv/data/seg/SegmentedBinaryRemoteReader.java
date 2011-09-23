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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.data.seg;

import org.apache.log4j.Logger;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.ResourceLocator;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author jrobinso
 */
public class SegmentedBinaryRemoteReader implements SegmentedBinaryReader {

    static Logger log = Logger.getLogger(SegmentedBinaryRemoteReader.class);

    String serverURL;

    String filePath;

    Map<String, String> attributes;

    public SegmentedBinaryRemoteReader(ResourceLocator locator) {
        this.serverURL = locator.getServerURL();
        this.filePath = locator.getPath();
    }

    public SegmentedChromosomeData getChromosomeData(String chr) {

        InputStream is = null;
        try {
            URL url = new URL(serverURL + "?method=getChromosomeData&file=" + filePath + "&chr=" + chr);
            is = HttpUtils.getInstance().openConnectionStream(url);
            SegmentedChromosomeData cd = new SegmentedChromosomeData();
            cd.deserialize(is);
            return cd;
        }
        catch (IOException ex) {
            log.error("Error opening file", ex);
            throw new RuntimeException(ex);
        } finally {
            if (is != null) {
                try {
                    is.close();

                } catch (IOException iOException) {
                    log.error("Error closing URL stream", iOException);
                }
            }

        }
    }

    public List<String> getSampleNames() {
        InputStream urlStream = null;

        try {
            List<String> childNames = new ArrayList(100);
            URL url = new URL(serverURL + "?method=getSampleNames&file=" + filePath);
            urlStream = HttpUtils.getInstance().openConnectionStream(url);
            BufferedReader reader = new BufferedReader(new InputStreamReader(urlStream));
            String nextLine = "";
            while ((nextLine = reader.readLine()) != null) {
                childNames.add(nextLine.trim());
            }
            reader.close();
            return childNames;

        } catch (IOException ex) {
            //log.error("Error in getChildNames", ex);   
            throw new RuntimeException(ex);
        } finally {
            if (urlStream != null) {
                try {
                    urlStream.close();

                } catch (IOException iOException) {
                    log.error("Error closing URL stream", iOException);
                }
            }

        }
    }
    /*
     *             
    List<String> lines = readAsStrings("data/attributes.txt");
    if (lines != null) {
    for (String kv : lines) {
    String[] tokens = kv.split("\t");
    attrs.put(tokens[0], tokens[1]);
    }
    }
     */

    public String getStringAttribute(String key) {
        if (attributes == null) {
            attributes = new HashMap();
            InputStream urlStream = null;

            try {
                URL url = new URL(serverURL + "?method=getAttributes&file=" + filePath);
                urlStream = HttpUtils.getInstance().openConnectionStream(url);
                BufferedReader reader = new BufferedReader(new InputStreamReader(urlStream));
                String nextLine = "";
                while ((nextLine = reader.readLine()) != null) {
                    String[] tokens = nextLine.split("=");
                    if (tokens.length > 1) {
                        attributes.put(tokens[0], tokens[1]);
                    }
                }

                reader.close();

            } catch (IOException ex) {
                log.error("Error in getChildNames", ex);

            } finally {
                if (urlStream != null) {
                    try {
                        urlStream.close();

                    } catch (IOException iOException) {
                        log.error("Error closing URL stream", iOException);
                    }
                }

            }
        }
        return attributes.get(key);

    }
}
