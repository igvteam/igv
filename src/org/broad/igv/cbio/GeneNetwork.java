/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTIES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package org.broad.igv.cbio;

import org.broad.igv.util.HttpUtils;
import org.broad.tribble.readers.AsciiLineReader;
import org.broad.tribble.readers.LineReader;
import org.jgrapht.EdgeFactory;

import java.io.IOException;
import java.io.InputStream;
import java.io.UnsupportedEncodingException;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLEncoder;

/**
 * Fetch a gene network from cBio portal. WORK IN PROGRESS.
 * User: jacob
 * Date: 2012/01/31
 */
public class GeneNetwork {

    static final String BASE_URL = "http://www.cbioportal.org/public-portal/webservice.do";
    static final String CMD = "getNetwork";
    static final String GENE_LIST_KEY = "gene_list";
    static final String BASIC_URL = BASE_URL + "?cmd=" + CMD + "&" + GENE_LIST_KEY + "=";

    LineReader getCBioData(String[] gene_list) throws IOException {

        try {
            String string_gl = URLEncoder.encode(gene_list[0], "UTF-8");
            for (int gi = 1; gi < gene_list.length; gi++) {
                string_gl += "," + URLEncoder.encode(gene_list[gi], "UTF-8");
            }

            String url = GeneNetwork.BASIC_URL + string_gl;
            System.out.println(url);
            InputStream is = HttpUtils.getInstance().openConnectionStream(new URL(url));
            LineReader reader = new AsciiLineReader(is);
            return reader;
        } catch (UnsupportedEncodingException e) {
            throw new IllegalArgumentException("Bad argument in genelist: " + e.getMessage());
        } catch (MalformedURLException e) {
            //It's not a malformed URL. There's essentially no way it could be,
            //unless the encoding malfunctions but throws no exception
            return null;
        }
    }

    public class testEdgeFactory implements EdgeFactory {

        public Object createEdge(Object o, Object o1) {
            return null; //TODO
        }
    }


}
