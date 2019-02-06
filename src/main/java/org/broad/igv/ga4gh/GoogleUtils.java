package org.broad.igv.ga4gh;

import org.apache.log4j.Logger;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.ui.util.MessageUtils;

import java.io.UnsupportedEncodingException;
import java.net.URLEncoder;

import static org.broad.igv.prefs.Constants.GOOGLE_PROJECT;
import static org.broad.igv.prefs.Constants.SAVE_GOOGLE_CREDENTIALS;

/**
 * Created by jrobinson on 7/15/16.
 */
public class GoogleUtils {

    private static Logger log = Logger.getLogger(GoogleUtils.class);

    private static String ProjectID;
    public static String GOOGLE_API_HOST = "www.googleapis.com";

    public static boolean isGoogleCloud(String url) {
        return url.startsWith("gs://") || url.contains(GOOGLE_API_HOST);
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
}
