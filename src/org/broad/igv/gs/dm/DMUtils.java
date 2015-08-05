/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broad.igv.gs.dm;


import biz.source_code.base64Coder.Base64Coder;
import com.google.gson.*;
import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.util.HttpUtils;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.net.URLEncoder;
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

    /**
     * Fetch the contents of the GenomeSpace directory.
     *
     * @param directoryURL
     * @return
     * @throws IOException
     */
    public static GSDirectoryListing getDirectoryListing(URL directoryURL) throws IOException{

        String str = HttpUtils.getInstance().getContentsAsJSON(directoryURL);
        JsonParser parser = new JsonParser();
        JsonObject obj = parser.parse(str).getAsJsonObject();

        JsonObject directory = obj.get("directory").getAsJsonObject();
        String dirUrlString = directory.get("url").getAsString();

        LinkedList<GSFileMetadata> elements = new LinkedList();
        if (obj.has("contents")) {
            Object c = obj.get("contents");
            List<JsonObject> contents = new ArrayList();
            if (c instanceof JsonObject) {
                contents.add((JsonObject) c);
            } else {
                JsonArray tmp = (JsonArray) c;
                int l = tmp.size();
                for (int i = 0; i < l; i++) {
                    contents.add((JsonObject) tmp.get(i));
                }
            }

            ArrayList<GSFileMetadata> dirElements = new ArrayList();
            ArrayList<GSFileMetadata> fileElements = new ArrayList();
            int contentsLength = contents.size();
            for (int i = 0; i < contentsLength; i++) {
                JsonObject o = contents.get(i);
                GSFileMetadata metaData = GSFileMetadata.deserializeElement(o);
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
        String contentType = "application/text";

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


    public static GSFileMetadata createDirectory(String putURL) throws IOException{

        JsonObject dirMeta = new JsonObject();

        dirMeta.add("isDirectory", new JsonPrimitive(true));
        System.out.println(dirMeta.toString());

        String body = "{\"isDirectory\":true}";
        String response = HttpUtils.getInstance().createGenomeSpaceDirectory(new URL(putURL), body);

        JsonParser parser = new JsonParser();
        JsonElement obj = parser.parse(response);
        return GSFileMetadata.deserializeElement(obj);

    }

    static void deleteFileOrDirectory(String delURL) throws IOException{
        HttpUtils.getInstance().delete(new URL(delURL));
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
