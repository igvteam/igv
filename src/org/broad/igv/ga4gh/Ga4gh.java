package org.broad.igv.ga4gh;

import com.google.gson.Gson;
import org.apache.commons.io.IOUtils;
import org.broad.igv.PreferenceManager;

import javax.net.ssl.HttpsURLConnection;
import java.io.BufferedReader;
import java.io.DataOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.HttpURLConnection;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.HashMap;

/**
 * Created by jrobinso on 6/16/14.
 * <p/>
 * Class to full around with the Google Global Alliance API
 * <p/>
 * TEST CASE
 * HG00096  (GBR low coverage)
 * 1:70,231,000-70,232,000
 * <p/>
 * readset id:  CJDmkYn8ChCcnc7i4KaWqmQ
 */
public class Ga4gh {


    // readset id = CJDmkYn8ChCcnc7i4KaWqmQ


    public static void main(String[] args) throws IOException {

        datasets();
        //readset();
        //readsetSearch();
        //reads();
    }

    public static void datasets() throws IOException {

        String authKey = PreferenceManager.getInstance().get(PreferenceManager.GOOGLE_API_KEY);

        String baseURL = "https://www.googleapis.com/genomics/v1beta/datasets/376902546192?key=" + authKey;


        // Create a URLConnection
        java.net.URLConnection connection = new java.net.URL(baseURL).openConnection();

        // Read the response
        java.io.BufferedReader br = new java.io.BufferedReader(new java.io.InputStreamReader(connection.getInputStream()));
        java.lang.StringBuffer sb = new java.lang.StringBuffer();
        java.lang.String str = br.readLine();
        while (str != null) {
            sb.append(str);
            str = br.readLine();
        }
        br.close();
        java.lang.String responseString = sb.toString();

        System.out.println(responseString);

    }

    public static void readset() throws IOException {

        String authKey = PreferenceManager.getInstance().get(PreferenceManager.GOOGLE_API_KEY);
        String baseURL = "https://www.googleapis.com/genomics/v1beta/readsets/CJDmkYn8ChCcnc7i4KaWqmQ?key=" + authKey;


        // Create a URLConnection
        java.net.URLConnection connection = new java.net.URL(baseURL).openConnection();

        // Read the response
        java.io.BufferedReader br = new java.io.BufferedReader(new java.io.InputStreamReader(connection.getInputStream()));
        java.lang.StringBuffer sb = new java.lang.StringBuffer();
        java.lang.String str = br.readLine();
        while (str != null) {
            sb.append(str);
            str = br.readLine();
        }
        br.close();
        java.lang.String responseString = sb.toString();

        System.out.println(responseString);

    }


    public static void reads() throws IOException {

        String authKey = PreferenceManager.getInstance().get(PreferenceManager.GOOGLE_API_KEY);
        URL baseURL = new URL("https://www.googleapis.com/genomics/v1beta/reads/search?key=" + authKey);
        String datasetId = "376902546192";

        // JSON content to post
        java.lang.String contentToPost = //"key=" + key + " " +
                "{readsetIds: [\"CJDmkYn8ChCcnc7i4KaWqmQ\"], " +
                        "sequenceName: 1, " +
                        "sequenceStart: 70231000, " +
                        "sequenceEnd: 70232000}";
        System.out.println(contentToPost);
        byte[] bytes = contentToPost.getBytes();

        // Create a URLConnection
        HttpURLConnection connection = (HttpURLConnection) baseURL.openConnection();
        connection.setUseCaches(false);
        connection.setDoInput(true);
        connection.setDoOutput(true);
        connection.setRequestMethod("POST");
        //connection.setRequestProperty("Content-Length", "" + bytes.length);
        connection.setRequestProperty("Content-Type", "application/json");
        connection.setRequestProperty("Cache-Control", "no-cache");

        // Post  content
        java.io.OutputStream stream = connection.getOutputStream();
        stream.write(bytes);
        stream.close();

        // Read the response
        java.io.BufferedReader br = new java.io.BufferedReader(new java.io.InputStreamReader(connection.getInputStream()));
        java.lang.String str;
        Gson gson = new Gson();
        while ((str = br.readLine()) != null) {
           // Object obj = gson.fromJson(str, HashMap.class);
            System.out.println(str);
        }
        br.close();


    }

    public static void readsetSearch() throws IOException {

        String authKey = PreferenceManager.getInstance().get(PreferenceManager.GOOGLE_API_KEY);
        URL baseURL = new URL("https://www.googleapis.com/genomics/v1beta/readsets/search?key=" + authKey);
        String datasetId = "376902546192";

        // JSON content to post
        java.lang.String contentToPost = //"key=" + key + " " +
                "{datasetIds: [" + datasetId + "]}";
        System.out.println(contentToPost);
        byte[] bytes = contentToPost.getBytes();

        // Create a URLConnection
        HttpURLConnection connection = (HttpURLConnection) baseURL.openConnection();
        connection.setUseCaches(false);
        connection.setDoInput(true);
        connection.setDoOutput(true);
        connection.setRequestMethod("POST");
        //connection.setRequestProperty("Content-Length", "" + bytes.length);
        connection.setRequestProperty("Content-Type", "application/json");
        connection.setRequestProperty("Cache-Control", "no-cache");

        // Post  content
        java.io.OutputStream stream = connection.getOutputStream();
        stream.write(bytes);
        stream.close();

        // Read the response
        java.io.BufferedReader br = new java.io.BufferedReader(new java.io.InputStreamReader(connection.getInputStream()));
        java.lang.String str;
        while ((str = br.readLine()) != null) {
            System.out.println(str);
        }
        br.close();


    }

}
