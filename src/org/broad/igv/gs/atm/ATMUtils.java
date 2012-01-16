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

import org.broad.igv.PreferenceManager;
import org.broad.igv.util.HttpUtils;
import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;
import org.json.JSONTokener;

import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;

/**
 * Utility class for accessing the GenomeSpace ATM webservice.
 * <p/>
 *
 * @author Jim Robinson
 * @date Aug 2, 2011
 */
public class ATMUtils {


    /**
     * Parse the contents of the URL as a JSON string encoding a list of WebToolDescriptor objects.
     *
     * @return
     * @throws IOException
     * @throws JSONException
     */
    public static List<WebToolDescriptor> getWebTools() throws IOException, JSONException {
        URL url = new URL(PreferenceManager.getInstance().get(PreferenceManager.GENOME_SPACE_ATM_SERVER) + "webtool/descriptor");
        String contents = HttpUtils.getInstance().getContentsAsString(url);
        JSONTokener tk = new JSONTokener(contents);
        JSONArray array = new JSONArray(tk);
       // JSONObject webDescArray = (JSONObject) array.get(1);
        List<WebToolDescriptor> webTools = parseWebtools(array);
        return webTools;
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
            JSONObject obj = webDescArray.getJSONObject(i);
           //JSONObject obj = ar.getJSONObject(1);

            String name = obj.getAsString("name");
            String id = obj.getAsString("internalId");
            String author = obj.getAsString("author");
            String description = obj.getAsString("description");
            String baseUrl = obj.getAsString("baseUrl");

            JSONArray fileParamArray = obj.getJSONArray("fileParameters");
            List<FileParameter> fileParams = parseFileParameters(fileParamArray);


            webTools.add(new WebToolDescriptor(name, id, author, description, baseUrl, fileParams));
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
                String stName = stObj.getAsString("name");
                String stId = stObj.getAsString("id");
                String stVersion = stObj.getAsString("version");
                String stAuthor = stObj.getAsString("author");
                String stDescription = stObj.getAsString("description");
                String stHelp = stObj.getAsString("help");
                String stUrlModifier = stObj.getAsString("urlModifier");
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
                JSONObject fObj = fileParamsArray.getJSONObject(n);
                String fpName = fObj.getAsString("name");
                String parentInternalId = fObj.getAsString("parentInternalId");
                String fpRequired = fObj.getAsString("required");
                String fpCompositeFilename = fObj.getAsString("compositeFilename");
                String fpNameDelimiters = fObj.getAsString("nameDelimiters");
                JSONArray formats = fObj.getJSONArray("formats");

                List<GSDataFormat> dataFormats = new ArrayList();
                if (formats.length() > 0) {
                    int nFormats = formats.length();
                    for (int f = 0; f < nFormats; f++) {
                        JSONObject format = formats.getJSONObject(f);
                        String fName = format.getAsString("name");
                        String fDescription = format.getAsString("description");
                        String fUrl = format.getAsString("url");
                        String fExt = format.getAsString("fileExtension");
                        dataFormats.add(new GSDataFormat(fName, fExt, fUrl));
                    }
                }

                fileParameters.add(new FileParameter(fpName, fpRequired,
                        fpCompositeFilename, fpNameDelimiters, dataFormats));
            }
        }

        return fileParameters;
    }

    public static String getWebtoolLaunchURL(String webtoolname) throws IOException, JSONException {
        URL url = new URL(PreferenceManager.getInstance().get(PreferenceManager.GENOME_SPACE_ATM_SERVER)+ "webtools/" + webtoolname);
        return HttpUtils.getInstance().getContentsAsString(url);
    }


    public static String getSubtoolLaunchURL(String webtoolname, String subtoolname) throws IOException, JSONException {
        URL url = new URL(PreferenceManager.getInstance().get(PreferenceManager.GENOME_SPACE_ATM_SERVER + "webtools/"  + webtoolname + "/" + subtoolname));
        return HttpUtils.getInstance().getContentsAsString(url);
    }

}
