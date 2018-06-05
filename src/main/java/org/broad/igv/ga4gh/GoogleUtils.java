package org.broad.igv.ga4gh;

import org.apache.log4j.Logger;
import org.broad.igv.exceptions.HttpResponseException;
import org.broad.igv.ui.action.LoadFromURLMenuAction;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.HttpUtils;

import java.io.IOException;
import java.io.InputStream;
import java.io.UnsupportedEncodingException;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLEncoder;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by jrobinson on 7/15/16.
 */
public class GoogleUtils {

    private static Logger log = Logger.getLogger(GoogleUtils.class);

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





}
