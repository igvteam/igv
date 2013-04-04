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

import com.google.gson.Gson;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonParser;
import org.broad.igv.PreferenceManager;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.StringUtils;

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

    private static final Gson gson = new Gson();

    /**
     * Parse the contents of the URL as a JSON string encoding a list of WebToolDescriptor objects.
     *
     * @return
     * @throws IOException
     */
    public static List<WebToolDescriptor> getWebTools() throws IOException {
        URL url = new URL(PreferenceManager.getInstance().get(PreferenceManager.GENOME_SPACE_ATM_SERVER) + "webtool/descriptor");
        String contents = HttpUtils.getInstance().getContentsAsJSON(url);

        JsonParser parser = new JsonParser();
        JsonArray array = parser.parse(contents).getAsJsonArray();
        List<WebToolDescriptor> webTools = parseWebtools(array);
        return webTools;
    }


    public static WebToolDescriptor getWebTool(String name) throws IOException{

        name = name.replace(" ", "%20");

        URL url = new URL(PreferenceManager.getInstance().get(PreferenceManager.GENOME_SPACE_ATM_SERVER) + "webtool/" +
                name + "/descriptor");
        String contents = HttpUtils.getInstance().getContentsAsJSON(url);

        JsonElement obj = (new JsonParser()).parse(contents);
        return parseWebTool(obj);
    }

    /**
     * Parse a WebToolDescriptor JSONArray
     *
     * @param webDescArray
     * @return
     */
    private static List<WebToolDescriptor> parseWebtools(JsonArray webDescArray){
        int count = webDescArray.size();
        List<WebToolDescriptor> webTools = new ArrayList();
        for (int i = 0; i < count; i++) {
            JsonElement obj = webDescArray.get(i);
            final WebToolDescriptor webToolDescriptor = parseWebTool(obj);

            webTools.add(webToolDescriptor);
        }
        return webTools;
    }

    private static WebToolDescriptor parseWebTool(JsonElement obj){
        return gson.fromJson(obj, WebToolDescriptor.class);
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
