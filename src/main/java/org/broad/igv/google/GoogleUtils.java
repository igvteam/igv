package org.broad.igv.google;

import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import org.apache.log4j.Logger;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.HttpUtils;

import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.net.URISyntaxException;
import java.net.URL;
import java.net.URLEncoder;

import static org.broad.igv.prefs.Constants.GOOGLE_PROJECT;
import static org.broad.igv.prefs.Constants.SAVE_GOOGLE_CREDENTIALS;

/**
 * Created by jrobinson on 7/15/16.
 */
public class GoogleUtils {

    private static Logger log = Logger.getLogger(GoogleUtils.class);

    private static String ProjectID;
    public static final String GOOGLE_API_HOST = "www.googleapis.com";
    public static final String GOOGLE_DRIVE_HOST = "drive.google.com";

    public static boolean isGoogleURL(String url) {
        return url != null && (isGoogleCloud(url) || isGoogleDrive(url) || isGoogleStorageURL(url));
    }

    public static boolean isGoogleCloud(String url) {
        return url != null && (url.startsWith("gs://") || url.contains(GOOGLE_API_HOST));
    }

    public static boolean isGoogleDrive(String url) {
        return url != null && (url.contains(GOOGLE_DRIVE_HOST) || url.contains("www.googleapis.com/drive"));
    }

    public static boolean isGoogleStorageURL(String url) {
        return url != null &&
                (url.startsWith("https://www.googleapis.com/storage") ||
                        url.startsWith("https://storage.cloud.google.com")  ||
                        url.startsWith("https://storage.googleapis.com"));
    }

    public static void checkLogin() {
        if (!OAuthUtils.getInstance().getProvider().isLoggedIn()) {
            OAuthUtils.getInstance().getProvider().doSecureLogin();
        }
    }


    /**
     * gs://igv-bam-test/NA12878.bam
     * https://www.googleapis.com/storage/v1/b/igv-bam-test/o/NA12878.bam
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
        String object = gsUrl.substring(i + 1);
        try {
            object = URLEncoder.encode(object, "UTF8");
        } catch (UnsupportedEncodingException e) {
            // This isn't going to happen
            log.error(e);
        }

        return "https://www.googleapis.com/storage/v1/b/" + bucket + "/o/" + object + "?alt=media";

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


    public static JsonObject getDriveFileInfo(String googleDriveURL) {

        try {
            String id = getGoogleDriveFileID(googleDriveURL);
            String endPoint = "https://www.googleapis.com/drive/v3/files/" + id + "?supportsTeamDrives=true";

            String json = HttpUtils.getInstance().getContentsAsJSON(new URL(endPoint));
            JsonParser parser = new JsonParser();
            JsonObject obj = parser.parse(json).getAsJsonObject();
            return obj;

        } catch (IOException e) {

            log.error("Error fetching google drive info", e);
            return null;
        }
    }


    public static void main(String[] args) throws IOException, URISyntaxException {

        String fileUrl = "https://drive.google.com/file/d/1Q4uEV2tv0aIyIvGqvwnOpBdmqVeFAPIw/view?usp=sharing";

        //OAuthUtils.getInstance().getProvider().openAuthorizationPage();
        JsonObject obj = GoogleUtils.getDriveFileInfo(fileUrl);

        System.out.println(obj);

    }
}
