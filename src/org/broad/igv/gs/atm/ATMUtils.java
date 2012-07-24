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

package org.broad.igv.gs.atm;

import org.broad.igv.PreferenceManager;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.StringUtils;
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
        String contents = HttpUtils.getInstance().getContentsAsJSON(url);
        JSONTokener tk = new JSONTokener(contents);
        JSONArray array = new JSONArray(tk);
        // JSONObject webDescArray = (JSONObject) array.get(1);
        List<WebToolDescriptor> webTools = parseWebtools(array);
        return webTools;
    }


    public static WebToolDescriptor getWebTool(String name) throws IOException, JSONException {

        name = name.replace(" ", "%20");

        URL url = new URL(PreferenceManager.getInstance().get(PreferenceManager.GENOME_SPACE_ATM_SERVER) + "webtool/" +
                name + "/descriptor");
        String contents = HttpUtils.getInstance().getContentsAsJSON(url);
        JSONTokener tk = new JSONTokener(contents);
        JSONObject obj = new JSONObject(tk);
        return parseWebTool(obj);
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

            final WebToolDescriptor webToolDescriptor = parseWebTool(obj);

            webTools.add(webToolDescriptor);
        }
        return webTools;
    }

    private static WebToolDescriptor parseWebTool(JSONObject obj) throws JSONException {
        String name = obj.getAsString("name");
        String id = obj.getAsString("internalId");
        String author = obj.getAsString("author");
        String description = obj.getAsString("description");
        String baseUrl = obj.getAsString("baseUrl");

        JSONArray fileParamArray = obj.getJSONArray("fileParameters");
        List<FileParameter> fileParams = parseFileParameters(fileParamArray);
        return new WebToolDescriptor(name, id, author, description, baseUrl, fileParams);
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


    /**
     * Return a launch URL with no parameters
     */
    public static String getWebtoolLaunchURL(WebToolDescriptor descriptor) throws IOException {

        return getWebtoolLaunchURL(descriptor, null);
    }

    /**
     * Notes:  Currently only a single file parameter is supported, if/when there is a GS client that accepts
     * multiple named file parameters we will deal with it then.
     *
     * @param descriptor
     * @param file
     * @return
     */
    public static String getWebtoolLaunchURL(WebToolDescriptor descriptor, String file) throws IOException {

        String name = descriptor.getName().replace(" ", "%20");
        String url = PreferenceManager.getInstance().get(PreferenceManager.GENOME_SPACE_ATM_SERVER) +
                "webtool/" + name + "/launchurl";

        List<FileParameter> fileParameters = descriptor.getFileParameters();
        if (file != null && fileParameters != null && fileParameters.size() > 0) {
            FileParameter param = fileParameters.get(0);
            url += "?" + param.getName() + "=" + StringUtils.encodeURL(file);
        }

        return HttpUtils.getInstance().getContentsAsString(new URL(url));
    }


}
