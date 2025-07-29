package org.broad.igv.util;

import org.broad.igv.logging.*;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.ui.util.MessageUtils;
import org.json.JSONObject;

import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.net.URL;
import java.net.URLEncoder;
import java.util.HashMap;
import java.util.Map;

import static org.broad.igv.prefs.Constants.GOOGLE_PROJECT;
import static org.broad.igv.prefs.Constants.SAVE_GOOGLE_CREDENTIALS;

/**
 * Created by jrobinson on 7/15/16.
 */
public class GoogleUtils {

    private static Logger log = LogManager.getLogger(GoogleUtils.class);

    /**
     * Google project ID for "user pays" requests.
     */
    private static String ProjectID;
    public static final String GOOGLE_API_HOST = "www.googleapis.com";
    public static final String GOOGLE_DRIVE_HOST = "drive.google.com";

    public static boolean isGoogleURL(String url) {
        return url != null && (isGoogleDrive(url) || isGoogleStorageURL(url));
    }


    public static boolean isGoogleDrive(String url) {
        return url != null && (url.contains(GOOGLE_DRIVE_HOST) || url.contains("www.googleapis.com/drive"));
    }

    public static boolean isGoogleStorageURL(String url) {
        return url != null &&
                (url.startsWith("gs://") ||
                        url.startsWith("https://www.googleapis.com/storage") ||
                        url.startsWith("https://storage.cloud.google.com") ||
                        url.startsWith("https://storage.googleapis.com"));
    }


    /**
     * gs://igv-bam-test/NA12878.bam
     * https://www.googleapis.com/storage/v1/b/igv-bam-test/o/NA12878.bam
     * <p>
     * https://storage.googleapis.com/download/storage/v1/b/BUCKET_NAME/o/OBJECT_NAME?alt=media
     *
     * @param gsUrl
     * @return
     */
    public static String translateGoogleCloudURL(String gsUrl) {

        int i = gsUrl.indexOf('/', 5);
        if (i < 0) {
            log.error("Invalid gs url: " + gsUrl);
            return gsUrl;
        }

        String bucket = gsUrl.substring(5, i);
        String object = googleObjectEncode(gsUrl.substring(i + 1));

        return "https://storage.googleapis.com/storage/v1/b/" + bucket + "/o/" + object + "?alt=media";

    }


    public static void enterGoogleProjectID() {
        String projectID = MessageUtils.showInputDialog("Enter Google project ID (for \"Requestor Pays\")",
                GoogleUtils.getProjectID());
        if (projectID != null) {
            GoogleUtils.setProjectID(projectID);

        }
    }


    public static String getProjectID() {
        if (ProjectID == null && PreferencesManager.getPreferences().getAsBoolean(SAVE_GOOGLE_CREDENTIALS)) {
            ProjectID = PreferencesManager.getPreferences().get(GOOGLE_PROJECT);
        }
        return ProjectID;
    }

    public static void setProjectID(String projectID) {
        ProjectID = projectID;
        if (ProjectID != null && PreferencesManager.getPreferences().getAsBoolean(SAVE_GOOGLE_CREDENTIALS)) {
            PreferencesManager.getPreferences().put(GOOGLE_PROJECT, projectID);
        }
    }

    public static String driveDownloadURL(String link) {

        // Return a google drive download url for the sharable link
        //https://drive.google.com/open?id=0B-lleX9c2pZFbDJ4VVRxakJzVGM
        //https://drive.google.com/file/d/1_FC4kCeO8E3V4dJ1yIW7A0sn1yURKIX-/view?usp=sharing

        String id = getGoogleDriveFileID(link);

        return id == null ? link :
                "https://www.googleapis.com/drive/v3/files/" + id + "?alt=media&supportsTeamDrives=true";
    }

    public static String getGoogleDriveFileID(String link) {

        //https://drive.google.com/file/d/1_FC4kCeO8E3V4dJ1yIW7A0sn1yURKIX-/view?usp=sharing
        if (link.contains("/open?id=")) {
            int i1 = link.indexOf("/open?id=") + 9;
            int i2 = link.indexOf("&");
            if (i1 > 0 && i2 > i1) {
                return link.substring(i1, i2);
            } else if (i1 > 0) {
                return link.substring(i1);
            }
        } else if (link.contains("/file/d/")) {
            int i1 = link.indexOf("/file/d/") + 8;
            int i2 = link.lastIndexOf("/");
            return link.substring(i1, i2);
        }
        return null;

    }


    public static JSONObject getDriveFileInfo(String googleDriveURL) {

        try {
            String id = getGoogleDriveFileID(googleDriveURL);
            String endPoint = "https://www.googleapis.com/drive/v3/files/" + id + "?supportsTeamDrives=true";

            String json = HttpUtils.getInstance().getContentsAsJSON(new URL(endPoint));
            org.json.JSONObject obj = new org.json.JSONObject(json);
            return obj;

        } catch (IOException e) {

            log.error("Error fetching google drive info", e);
            return null;
        }
    }

    /**
     * Percent a GCS object name. See https://cloud.google.com/storage/docs/request-endpoints
     * Specific characters to encode: !, #, $, &, ', (, ), *, +, ,, /, :, ;, =, ?, @, [, ], and space characters.
     * @param objectName
     * @return
     */

    static String googleObjectEncode(String objectName) {

        String result = "";
        for (int i = 0; i < objectName.length(); i++) {
            Character c = objectName.charAt(i);
            if (encodings.containsKey(c)) {
                result += encodings.get(c);
            } else {
                result += c;
            }
        }
        return result;
    }

    //	%23	%24	%25	%26	%27	%28	%29	%2A	%2B	%2C	%2F	%3A	%3B	%3D	%3F	%40	%5B	%5D
    static Map<Character, String> encodings = new HashMap<>();

    static {
        encodings.put(' ', "%20");
        encodings.put('!', "%21");
        encodings.put('#', "%23");
        encodings.put('$', "%24");
        encodings.put('%', "%25");
        encodings.put('&', "%26");
        encodings.put('\'', "%27");
        encodings.put('(', "%28");
        encodings.put(')', "%29");
        encodings.put('*', "%2A");
        encodings.put('+', "%2B");
        encodings.put(',', "%2C");
        encodings.put('/', "%2F");
        encodings.put(':', "%3A");
        encodings.put(';', "%3B");
        encodings.put('=', "%3D");
        encodings.put('?', "%3F");
        encodings.put('@', "%40");
        encodings.put('[', "%5B");
        encodings.put(']', "%5D");

    }

}
