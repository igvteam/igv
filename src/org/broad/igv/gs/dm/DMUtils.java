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


import biz.source_code.base64Coder.Base64Coder;
import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.util.HttpUtils;
import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;
import org.json.JSONTokener;

import java.io.*;
import java.net.*;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.*;

/**
 * Utility class for accessing the GenomeSpace data manager web service
 *
 * @author Jim Robinson
 * @date Aug 2, 2011
 */
public class DMUtils {

    private static Logger log = Logger.getLogger(DMUtils.class);
    private static final String UPLOAD_SERVICE = "uploadurl";
    public static final String DEFAULT_DIRECTORY = "defaultdirectory";
    public static final String PERSONAL_DIRECTORY = "personaldirectory";


    public static GSDirectoryListing listDefaultDirectory() {
        try {
            URL defaultURL = new URL(PreferenceManager.getInstance().get(PreferenceManager.GENOME_SPACE_DM_SERVER) +
                    DEFAULT_DIRECTORY);
            return getDirectoryListing(defaultURL);
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Fetch the contents of the GenomeSpace directory.
     *
     * @param directoryURL
     * @return
     * @throws IOException
     * @throws JSONException
     */
    public static GSDirectoryListing getDirectoryListing(URL directoryURL) throws IOException, JSONException {

        String str = HttpUtils.getInstance().getContentsAsJSON(directoryURL);
        JSONTokener tk = new JSONTokener(str);
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
                GSFileMetadata metaData = new GSFileMetadata(o);
                if (metaData.isDirectory()) {
                    dirElements.add(metaData);
                } else {
                    fileElements.add(metaData);
                }
            }

            elements.addAll(dirElements);
            elements.addAll(fileElements);
        }

        return new GSDirectoryListing(dirUrlString, elements);

    }


    /**
     * Upload a file to GenomeSpace.
     *
     * @param localFile
     * @param gsPath    the relative path in the users GenomeSpace account
     * @throws IOException
     * @throws URISyntaxException
     */
    public static void uploadFile(File localFile, String gsPath) throws IOException, URISyntaxException {


        byte[] md5 = computeMD5(localFile);
        String base64String = new String(Base64Coder.encode(md5));
        long contentLength = localFile.length();
        String contentType = "application/text"; //"text";

        String tmp = PreferenceManager.getInstance().get(PreferenceManager.GENOME_SPACE_DM_SERVER) + UPLOAD_SERVICE +
                gsPath + "?Content-Length=" + contentLength +
                "&Content-MD5=" + URLEncoder.encode(base64String, "UTF-8") + "&Content-Type=" + contentType;

        String uploadURL = HttpUtils.getInstance().getContentsAsJSON(new URL(tmp));

        Map<String, String> headers = new HashMap<String, String>();
        headers.put("Content-MD5", base64String);
        headers.put("Content-Length", String.valueOf(contentLength));
        headers.put("Content-Type", contentType);


        HttpUtils.getInstance().uploadGenomeSpaceFile(uploadURL, localFile, headers);

    }


    public static GSFileMetadata createDirectory(String putURL) throws IOException, JSONException {

        JSONObject dirMeta = new JSONObject();
        try {
            dirMeta.put("isDirectory", true);
            System.out.println(dirMeta.toString());
        } catch (JSONException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

        String body = "{\"isDirectory\":true}";
        String response = HttpUtils.getInstance().createGenomeSpaceDirectory(new URL(putURL), body);

        JSONTokener tk = new JSONTokener(response);
        JSONObject obj = new JSONObject(tk);
        return new GSFileMetadata(obj);

    }


    /**
     * Compute the MD5 hash for the given file.
     *
     * @param file
     * @return
     * @throws IOException
     */
    public static byte[] computeMD5(File file) throws IOException {
        BufferedInputStream in = null;
        try {
            MessageDigest md = MessageDigest.getInstance("MD5");
            in = new BufferedInputStream(new FileInputStream(file));
            int theByte = 0;
            while ((theByte = in.read()) != -1) {
                md.update((byte) theByte);
            }
            return md.digest();

        } catch (NoSuchAlgorithmException e) {
            // Should be impossible to get here
            log.error("Error creating MD5 digest");
            throw new RuntimeException(e);
        } finally {
            if (in != null) in.close();
        }


    }

    /**
     * Convert an array of bytes to a hex string.
     *
     * @param v
     * @return
     */
    public static String toHexString(byte[] v) {
        final String HEX_DIGITS = "0123456789abcdef";
        StringBuffer sb = new StringBuffer(v.length * 2);
        for (int i = 0; i < v.length; i++) {
            int b = v[i] & 0xFF;
            sb.append(HEX_DIGITS.charAt(b >>> 4)).append(HEX_DIGITS.charAt(b & 0xF));
        }
        return sb.toString();
    }


    public static void main(String[] args) throws IOException, URISyntaxException {
        final String testFile = "/Users/jrobinso/projects/igv/test/data/bed/Unigene.sample.bed";
        final File localFile = new File(testFile);
        uploadFile(localFile, "/users/test/Unigene.sample.bed");
        System.exit(-1);

        /*
        File file = new File(testFile);

        byte[] md5 = computeMD5(file);

        String hexString = toHexString(md5);
        System.out.println(hexString);

        String base64String = (new BASE64Encoder()).encode(md5);
        System.out.println(base64String);    */
    }

}
