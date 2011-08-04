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

import com.sun.tools.corba.se.idl.StringGen;
import org.apache.commons.codec.binary.Base64;
import org.broad.igv.gs.GSUtils;
import org.broad.igv.util.IGVHttpClientUtils;
import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;
import org.json.JSONTokener;
import sun.java2d.pipe.BufferedTextPipe;
import sun.misc.BASE64Encoder;

import java.io.*;
import java.net.URI;
import java.net.URISyntaxException;
import java.net.URL;
import java.net.URLEncoder;
import java.security.DigestException;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.*;

/**
 * Utility class for accessing the GenomeSpace data manager web service
 *
 * @author Jim Robinson
 * @date Aug 2, 2011
 */
public class GSDMUtils {

    public static String baseUrl = "https://dmtest.genomespace.org:8444/datamanager/";

    /**
     * Fetch the contents of the GenomeSpace directory.
     *
     * @param directoryURL
     * @return
     * @throws IOException
     * @throws JSONException
     */
    public static GSDirectoryListing getDirectoryListing(URL directoryURL) throws IOException, JSONException {

        String str = IGVHttpClientUtils.getContentsAsString(directoryURL);
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

    }


    public static void uploadFile(File localFile, String gsPath) throws IOException, URISyntaxException {


        byte[] md5 = computeMD5(localFile);
        String base64String = (new BASE64Encoder()).encode(md5);
        String hexString = toHexString(md5);
        long contentLength = localFile.length();
        String contentType = "text";

        // TODO -- is this the correct way to get the user?
        String user = GSUtils.getCachedUsernameForSSO();

        String tmp = baseUrl + "uploadurls/users/" + user + "/" + gsPath + "?Content-Length=" + contentLength +
                "&Content-MD5=" + URLEncoder.encode(base64String, "UTF-8") + "&Content-Type=" + contentType;

        String uploadURL = IGVHttpClientUtils.getContentsAsString(new URL(tmp));

        URI uri = new URI(uploadURL);

        Map<String, String> headers = new HashMap<String, String>();
        headers.put("Content-MD5", base64String);
        headers.put("x-amz-meta-md5-hash", hexString);

        IGVHttpClientUtils.uploadFile(uri, localFile, headers);

    }

    private static byte[] computeMD5(File file) throws IOException {
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
            // Basically impossible to get here
            e.printStackTrace();
            return null;
        } finally {
            if (in != null) in.close();
        }


    }

    private static String toHexString(byte[] v) {
        final String HEX_DIGITS = "0123456789abcdef";
        StringBuffer sb = new StringBuffer(v.length * 2);
        for (int i = 0; i < v.length; i++) {
            int b = v[i] & 0xFF;
            sb.append(HEX_DIGITS.charAt(b >>> 4)).append(HEX_DIGITS.charAt(b & 0xF));
        }
        return sb.toString();
    }


    public static void main(String[] args) throws IOException, URISyntaxException {
        final String testFile = "/Users/jrobinso/projects/igv/259.wgs.muc1.hist.txt";
        final File localFile = new File(testFile);
        uploadFile(localFile, "259.wgs.muc1.hist.txt");
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
