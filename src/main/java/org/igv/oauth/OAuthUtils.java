package org.igv.oauth;

import org.igv.logging.*;
import org.igv.DirectoryManager;
import org.igv.prefs.Constants;
import org.igv.prefs.PreferencesManager;
import org.igv.ui.IGVMenuBar;
import org.igv.ui.util.MessageUtils;
import org.igv.util.AmazonUtils;
import org.igv.util.HttpMappings;
import org.igv.util.HttpUtils;
import org.json.JSONArray;
import org.json.JSONObject;


import java.io.*;
import java.net.URL;
import java.net.URLConnection;
import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.zip.GZIPInputStream;

/**
 * Created by jrobinso on 11/19/14.
 */
public class OAuthUtils {

    private static Logger log = LogManager.getLogger(OAuthUtils.class);

    private static OAuthUtils theInstance;

    static OAuthProvider awsProvider;

    static OAuthProvider googleProvider;

    static Map<String, OAuthProvider> providerCache;

    public static synchronized OAuthUtils getInstance() {
        if (theInstance == null) {
            theInstance = new OAuthUtils();
        }
        return theInstance;
    }

    private OAuthUtils() {
        try {
            providerCache = new LinkedHashMap<>();   // Ordered (linked) map is important
            fetchOauthConfigs();
        } catch (Exception e) {
            log.error("Error fetching oAuth properties", e);
        }
    }

    /**
     * Fetch user-configured oAuth configurations, if any*
     * @throws IOException
     */
    private void fetchOauthConfigs() throws IOException {

        // Load a provider config specified in preferences
        String provisioningURL = PreferencesManager.getPreferences().getProvisioningURL();
        log.debug("The provisioning URL from prefs.properties is: " + provisioningURL);
        if (provisioningURL != null && provisioningURL.length() > 0) {
            String json = loadAsString(provisioningURL);
            parseProviderJson(json);
        }

        // Local config takes precendence, overriding URL provisioned and Broad's default
        String oauthConfig = DirectoryManager.getIgvDirectory() + "/oauth-config.json";
        if (!(new File(oauthConfig)).exists()) {
            oauthConfig = DirectoryManager.getIgvDirectory() + "/oauth-config-custom.json";
        }
        if ((new File(oauthConfig)).exists()) {
            try {
                log.debug("Loading Oauth properties from: " + oauthConfig);
                String json = loadAsString(oauthConfig);
                parseProviderJson(json);
            } catch (IOException e) {
                log.error(e);
            }
        }
    }

    /**
     * Called by AWS code only
     */
    public OAuthProvider getAWSProvider() {
        if (awsProvider == null) {
            throw new RuntimeException("AWS Oauth is not configured");
        }
        return awsProvider;
    }

    public OAuthProvider getGoogleProvider() {
        if (googleProvider == null) {
            try {
                log.info("Loading Google oAuth properties");
                googleProvider = loadDefaultOauthProperties();
                if (googleProvider != null && IGVMenuBar.getInstance() != null) {
                    IGVMenuBar.getInstance().enableGoogleMenu(true);
                }
            } catch (IOException e) {
                log.error("Error loading Google oAuth properties", e);
                MessageUtils.showErrorMessage("Error loading Google oAuth properties", e);
            }
        }
        return googleProvider;
    }


    /**
     * Load the default (Google) oAuth properties
     *
     * @throws IOException
     */
    private OAuthProvider loadDefaultOauthProperties() throws IOException {
        String propertiesURL = PreferencesManager.getPreferences().get(Constants.PROVISIONING_URL_DEFAULT);
        if (propertiesURL != null && propertiesURL.length() > 0) {
            String json = loadAsString(propertiesURL);
            JSONObject obj = new JSONObject(json);
            return parseProviderObject(obj);
        } else {
            return null;
        }
    }

    /**
     * Called from preferences editor
     *
     * @param provisioningURL
     * @throws IOException
     */
    public void updateOauthProvider(String provisioningURL) throws IOException {
        if (provisioningURL != null && provisioningURL.trim().length() > 0) {
            String json = loadAsString(provisioningURL);
            parseProviderJson(json);
        }
    }


    private String loadAsString(String urlOrPath) throws IOException {
        InputStream is = null;
        try {
            is = openInputStream(urlOrPath);
            byte[] bytes = is.readAllBytes();
            if ((bytes[0] == (byte) (GZIPInputStream.GZIP_MAGIC)) && (bytes[1] == (byte) (GZIPInputStream.GZIP_MAGIC >> 8))) {
                bytes = new GZIPInputStream(new ByteArrayInputStream(bytes)).readAllBytes();
            }

            return new String(bytes, "UTF-8");
        } finally {
            if (is != null) is.close();
        }
    }

    /**
     * Parse Oauth provider configuration and update state
     * <p>
     *
     * @param json
     * @throws IOException
     */
    private void parseProviderJson(String json) throws IOException {

        try {
            JSONArray jsonArray = new org.json.JSONArray(json);
            for (int i = 0; i < jsonArray.length(); i++) {
                parseProviderObject(jsonArray.getJSONObject(i));
            }
        } catch (org.json.JSONException e) {
            JSONObject jsonObject = new org.json.JSONObject(json);
            parseProviderObject(jsonObject);
        }
    }

    private OAuthProvider parseProviderObject(JSONObject obj) throws IOException {
        OAuthProvider p = new OAuthProvider(obj);
        providerCache.put(p.getState(), p);
        if ((obj.has("auth_provider")
                && obj.getString("auth_provider").equals("Amazon")) ||
                obj.has("aws_region")) {
            awsProvider = p;
            AmazonUtils.setCognitoConfig(obj);
        } else if (p.isGoogle()) {
            googleProvider = p;
            googleProvider.setAuthProvider("Google");
        }
        return p;
    }

    /**
     * Called during authorization flow from CommandListener
     *
     * @param state
     * @return
     * @throws IOException
     */
    public OAuthProvider getProviderForState(String state) throws IOException {
        if (providerCache.containsKey(state)) {
            return providerCache.get(state);
        } else {
            // This should never happen, perhaps an error should be thrown.
            log.warn("No oAuth provider found for callback");
            return null;
        }
    }

    public OAuthProvider getProviderForURL(URL url) throws IOException {
        for (OAuthProvider provider : providerCache.values()) {
            if (provider.appliesToUrl(url)) {
                if (provider.isGoogle()) {
                    IGVMenuBar.getInstance().enableGoogleMenu(true);
                }
                return provider;
            }
        }
        return null;
    }

    public Collection<OAuthProvider> getAllProviders() {
        return providerCache.values();
    }

    /**
     * Open an input stream for reading a local or remote file.
     *
     * @param urlOrPath -- either an http URL or path to a local file.  Can be gzipped
     * @return
     * @throws IOException
     */
    public static InputStream openInputStream(String urlOrPath) throws IOException {

        InputStream inputStream = null;
        if (urlOrPath.startsWith("http://") || urlOrPath.startsWith("https://")) {
            URL url = new URL(urlOrPath);
            URLConnection conn = HttpUtils.getInstance().openProxiedConnection(url);
            inputStream = conn.getInputStream();
        } else {
            if (urlOrPath.startsWith("file://")) {
                urlOrPath = urlOrPath.substring(7);
            }
            File file = new File(urlOrPath);
            inputStream = new FileInputStream(file);
        }

        return inputStream;

    }


}
