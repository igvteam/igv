

package org.broad.igv.sam.mutreview;

import com.google.gson.*;
import org.apache.log4j.Logger;
import org.broad.igv.ga4gh.OAuthUtils;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.MessageUtils;

import javax.imageio.ImageIO;
import javax.swing.*;
import java.awt.image.BufferedImage;
import java.io.*;
import java.net.HttpURLConnection;
import java.net.URL;
import java.util.*;

/**
 * Helper class for  Google Cloud Storage API
 * <p/>
 * <p/>
 * Created by jrobinso on 11/3/17.
 */
public class GoogleCloudStorageHelper {

    private static Logger log = Logger.getLogger(GoogleCloudStorageHelper.class);


    public static void upload(BufferedImage image, VariantReviewMetadata metadata) throws IOException {

        String filename = "test_" + metadata.userId + "_" + metadata.timestamp + "_" + metadata.score;

        String r1 = uploadJpeg(image, filename);
        log.info(r1);

        String json = (new Gson()).toJson(metadata);
        String r2 = uploadJson(json, filename);
        log.info(r2);

    }

    private static String uploadJpeg(BufferedImage image, String filename) throws IOException {

        ByteArrayOutputStream stream = new ByteArrayOutputStream(2000000);
        ImageIO.write(image, "jpg", stream);
        byte [] bytes = stream.toByteArray();
        return doPost(bytes, "image/jpeg", filename + ".jpeg");
    }

    private static String uploadJson(String json, String filename) throws IOException {
        return doPost(json.getBytes(), "text/plain", filename + ".json");
    }


    private static String doPost(byte [] bytes, String contentType, String filename) throws IOException {

        String token = OAuthUtils.getInstance().getAccessToken();

        String fullUrl = "https://www.googleapis.com/upload/storage/v1/b/igv-screenshots/o?uploadType=media&name=" +
                filename;


        URL url = new URL(fullUrl);


        // Create a URLConnection
        HttpURLConnection connection = null;
        BufferedReader br;
        OutputStream outputStream;

        try {
            connection = (HttpURLConnection) url.openConnection();
            connection.setUseCaches(false);
            connection.setDoInput(true);
            connection.setDoOutput(true);
            connection.setRequestMethod("POST");
            connection.setRequestProperty("Content-Length", "" + bytes.length);
            connection.setRequestProperty("Content-Type", contentType);


            if (token != null) {
                connection.setRequestProperty("Authorization", "Bearer " + token);
            }

            // Post  content
            outputStream = connection.getOutputStream();
            outputStream.write(bytes);
            outputStream.close();

            // Read the response
            br = new BufferedReader(new InputStreamReader(connection.getInputStream()));
            StringBuffer sb = new StringBuffer();
            String str = br.readLine();
            while (str != null) {
                sb.append(str);
                str = br.readLine();
            }
            br.close();

            return sb.toString();
        } catch (Exception e) {

                handleHttpException(url, connection, e);


            return null;
        }
    }

    static void handleHttpException(URL url, HttpURLConnection connection, Exception e) throws IOException {
        int rs = connection.getResponseCode();
        String sb = getErrorMessage(connection);
        if (sb != null && sb.length() > 0) {
            MessageUtils.showErrorMessage(sb, e);
        } else if (rs == 404) {
            MessageUtils.showErrorMessage("The requested resource was not found<br>" + url, e);
        } else if (rs == 401 || rs == 403) {
            displayAuthorizationDialog(url.getHost());
        } else {
            MessageUtils.showErrorMessage("Error accessing resource", e);
            log.error("Error accessing resource", e);
        }
    }

    private static String getErrorMessage(HttpURLConnection connection) throws IOException {
        BufferedReader br = new BufferedReader(new InputStreamReader(connection.getErrorStream()));
        StringBuffer sb = new StringBuffer();
        String str = br.readLine();
        while (str != null) {
            sb.append(str);
            str = br.readLine();
        }
        br.close();

        JsonParser parser = new JsonParser();
        JsonObject obj = parser.parse(sb.toString()).getAsJsonObject();

        JsonObject errorObject = obj.getAsJsonObject("error");
        if (errorObject != null) {
            JsonPrimitive msg = errorObject.getAsJsonPrimitive("message");
            if (msg != null) return msg.getAsString();
        }


        return sb.toString();
    }

    static void displayAuthorizationDialog(String host) {

        String message = "The requested resource at '" + host + "' requires authorization.";
        Icon icon = null;
        int option = JOptionPane.showOptionDialog(IGV.getMainFrame(),
                message,
                "Error",
                JOptionPane.OK_CANCEL_OPTION,
                JOptionPane.QUESTION_MESSAGE,
                icon,
                new String[]{"Cancel", "Authorize"},
                JOptionPane.YES_OPTION
        );

        if (option == 1) {
            try {
                OAuthUtils.getInstance().openAuthorizationPage();
            } catch (Exception e) {
                MessageUtils.showErrorMessage("Error fetching oAuth token", e);
                log.error("Error fetching oAuth tokens", e);
            }
        }

    }


    static Map<String, String> genomeIdMap = new HashMap<String, String>();  // A hack until readset meta data is available

    static {
        genomeIdMap = new HashMap<String, String>();
        genomeIdMap.put("Google 10473108253681171589", "hg19");
        genomeIdMap.put("Google 383928317087", "hg19");
        genomeIdMap.put("Google 461916304629", "hg19");
        genomeIdMap.put("Google 337315832689", "hg19");

    }


}
