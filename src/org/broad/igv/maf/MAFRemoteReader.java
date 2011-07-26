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
package org.broad.igv.maf;

import org.apache.log4j.Logger;
import org.broad.igv.util.IGVHttpClientUtils;
import org.broad.igv.util.ResourceLocator;

import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.net.URLConnection;
import java.util.List;
import java.util.zip.GZIPInputStream;

/**
 * @author jrobinso
 */
public class MAFRemoteReader implements MAFReader {

    static Logger log = Logger.getLogger(MAFRemoteReader.class);
    String serverURL = "http://www.broadinstitute.org/webservices/igv";
    static MAFTileCodec codec = new MAFTileCodec();

    public MAFRemoteReader(ResourceLocator locator) {
        // TODO -- set server URL and path to MAF directory
    }

    public MAFTile loadTile(String chr, int start, int end,
                            List<String> species) {
        InputStream is = null;
        try {
            URL url = new URL(serverURL + "?method=maf&chr=" + chr + "&start=" + start + "&end=" + end);
            is = IGVHttpClientUtils.openConnectionStream(url);
            DataInputStream dis = new DataInputStream(new GZIPInputStream(new BufferedInputStream(is)));
            MAFTile tile = codec.decode(dis);
            return tile;
        } catch (IOException ex) {
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
}
