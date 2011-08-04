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

package org.broad.igv.gs.atm;

import org.broad.igv.util.IGVHttpClientUtils;
import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;
import org.json.JSONTokener;

import java.io.BufferedInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;

/**
 * Utility class for accessing the GenomeSpace ATM webservice.
 * <p/>
 * https://atmtest.genomespace.org:8443/atm/webtools
 *
 * @author Jim Robinson
 * @date Aug 2, 2011
 */
public class GSATMUtils {

    static final String BASE_URL = "https://atmtest.genomespace.org:8443/atm/webtools";

    public static void main(String[] args) throws IOException, JSONException {

        String url = getWebtoolLaunchURL("Cytoscape");
        System.out.println(url);

        List<WebToolDescriptor> webTools = getWebTools(new URL(BASE_URL));
        for (WebToolDescriptor wt : webTools) {
            wt.print();
        }
        System.exit(-1);
    }

    /**
     * Parse the contents of the URL as a JSON string encoding a list of WebToolDescriptor objects.
     *
     * @param url
     * @return
     * @throws IOException
     * @throws JSONException
     */
    public static List<WebToolDescriptor> getWebTools(URL url) throws IOException, JSONException {

        StringBuffer buf = new StringBuffer();
        InputStream is = null;
        try {
            is = IGVHttpClientUtils.openConnectionStream(url);
            BufferedInputStream bis = new BufferedInputStream(is);
            int b;
            while ((b = bis.read()) >= 0) {
                buf.append((char) b);
            }
            JSONTokener tk = new JSONTokener(buf.toString());
            JSONArray array = new JSONArray(tk);
            JSONArray webDescArray = (JSONArray) array.get(1);
            List<WebToolDescriptor> webTools = parseWebtools(webDescArray);
            return webTools;
        } finally {
            is.close();
        }
    }

    /**
     * Parse a WebToolDescriptor JSONArray
     *
     * @param webDescArray
     * @return
     * @throws JSONException
     */
    private static List<WebToolDescriptor> parseWebtools(JSONArray webDescArray) throws JSONException {
        int count = webDescArray.length();
        List<WebToolDescriptor> webTools = new ArrayList();
        for (int i = 0; i < count; i++) {
            JSONArray ar = webDescArray.getJSONArray(i);
            JSONObject obj = ar.getJSONObject(1);

            String name = obj.get("name").toString();
            String id = obj.get("id").toString();
            String version = obj.get("version").toString();
            String author = obj.get("author").toString();
            String description = obj.get("description").toString();
            String help = obj.get("help").toString();
            String baseUrl = obj.get("baseUrl").toString();

            JSONArray fileParamArray = obj.getJSONArray("fileParameters").getJSONArray(1);
            List<FileParameter> fileParams = parseFileParameters(fileParamArray);

            JSONArray subToolsArray = obj.getJSONArray("subTools").getJSONArray(1);
            List<SubToolDescriptor> subTools = parseSubtools(subToolsArray);

            webTools.add(new WebToolDescriptor(name, id, version, author, description, help, baseUrl,
                    fileParams, subTools));
        }
        return webTools;
    }

    /**
     * Parse a SubToolDescriptor JSONArray
     *
     * @param subToolsArray
     * @return
     * @throws JSONException
     */
    private static List<SubToolDescriptor> parseSubtools(JSONArray subToolsArray) throws JSONException {

        List<SubToolDescriptor> subtoolDescriptors = new ArrayList();
        if (subToolsArray.length() > 0) {
            int nTools = subToolsArray.length();
            for (int n = 0; n < nTools; n++) {
                JSONObject stObj = subToolsArray.getJSONArray(n).getJSONObject(1);
                String stName = stObj.get("name").toString();
                String stId = stObj.get("id").toString();
                String stVersion = stObj.get("version").toString();
                String stAuthor = stObj.get("author").toString();
                String stDescription = stObj.get("description").toString();
                String stHelp = stObj.get("help").toString();
                String stUrlModifier = stObj.get("urlModifier").toString();
                JSONArray stFileParams = stObj.getJSONArray("fileParameters").getJSONArray(1);
                List<FileParameter> fileParams = parseFileParameters(stFileParams);

                subtoolDescriptors.add(new SubToolDescriptor(stName, stId, stVersion, stAuthor,
                        stDescription, stHelp, stUrlModifier, fileParams));
            }
        }
        return subtoolDescriptors;
    }

    /**
     * Parse a FileParameter JSONArray
     *
     * @param fileParamsArray
     * @return
     * @throws JSONException
     */
    private static List<FileParameter> parseFileParameters(JSONArray fileParamsArray) throws JSONException {

        List<FileParameter> fileParameters = new ArrayList();

        if (fileParamsArray.length() > 0) {
            int nFileParams = fileParamsArray.length();
            for (int n = 0; n < nFileParams; n++) {
                JSONObject fObj = fileParamsArray.getJSONArray(n).getJSONObject(1);
                String fpName = fObj.get("name").toString();
                String fpDescription = fObj.get("description").toString();
                String fpRequired = fObj.get("required").toString();
                String fpCompositeFilename = fObj.get("compositeFilename").toString();
                String fpNameDelimiters = fObj.get("nameDelimiters").toString();
                JSONArray formats = fObj.getJSONArray("formats");

                List<GSDataFormat> dataFormats = new ArrayList();
                if (formats.length() > 0) {
                    final JSONArray formatsArray = formats.getJSONArray(1);
                    int nFormats = formatsArray.length();
                    for (int f = 0; f < nFormats; f++) {
                        JSONObject format = formatsArray.getJSONArray(f).getJSONObject(1);
                        String fName = format.get("name").toString();
                        String fVersion = format.get("version").toString();
                        String fUrl = format.get("url").toString();
                        dataFormats.add(new GSDataFormat(fName, fVersion, fUrl));
                    }
                }

                fileParameters.add(new FileParameter(fpName, fpDescription, fpRequired,
                        fpCompositeFilename, fpNameDelimiters, dataFormats));
            }
        }

        return fileParameters;
    }

    public static String getWebtoolLaunchURL(String webtoolname) throws IOException, JSONException {
        URL url = new URL(BASE_URL + "/" + webtoolname);

        StringBuffer buf = new StringBuffer();
        InputStream is = null;
        try {
            is = IGVHttpClientUtils.openConnectionStream(url);
            BufferedInputStream bis = new BufferedInputStream(is);
            int b;
            while ((b = bis.read()) >= 0) {
                buf.append((char) b);
            }
            return buf.toString();
        } finally {
            is.close();
        }
    }

}
