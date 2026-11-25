package org.igv.util;

import org.igv.AbstractHeadlessTest;
import org.igv.prefs.PreferencesManager;
import org.junit.After;
import org.junit.Test;

import java.io.FileNotFoundException;
import java.net.HttpURLConnection;
import java.net.URL;
import java.util.Collections;
import java.util.Map;

import static org.junit.Assert.assertEquals;

/**
 * Tests substitution of the UCSC backup host when the primary host cannot be reached.  The failure is simulated by
 * setting the connect timeout to 1 ms, which no real connection can meet.
 */
public class UCSCBackupHostTest extends AbstractHeadlessTest {

    private static final String UCSC_HOST = "hgdownload.soe.ucsc.edu";
    private static final String BACKUP_HOST = "genome-browser.s3.us-east-1.amazonaws.com";
    private static final String UCSC_URL = "https://" + UCSC_HOST + "/goldenPath/hg38/bigZips/hg38.2bit";

    // Request a single byte, we're interested in the connection not the content
    private static final Map<String, String> RANGE_HEADER = Collections.singletonMap("Range", "bytes=0-0");

    @After
    public void resetBackupHost() {
        HttpUtils.getInstance().resetUCSCBackupHost();
        // Restore the default timeout, so a failed test does not leave the 1 ms timeout in place for the next
        PreferencesManager.getPreferences().remove("UCSC_CONNECT_TIMEOUT");
    }

    /**
     * A connection failure switches to the backup host, and the substitution is remembered for subsequent requests.
     */
    @Test
    public void testConnectTimeoutSwitchesToBackupHost() throws Exception {

        PreferencesManager.getPreferences().put("UCSC_CONNECT_TIMEOUT", "1");

        HttpURLConnection conn = HttpUtils.getInstance().openConnection(new URL(UCSC_URL), RANGE_HEADER);
        assertEquals(BACKUP_HOST, conn.getURL().getHost());
        conn.disconnect();

        // Restore a workable timeout.  The primary host would now succeed, but the cached substitution should
        // send us to the backup host without attempting it.
        PreferencesManager.getPreferences().put("UCSC_CONNECT_TIMEOUT", "10000");

        conn = HttpUtils.getInstance().openConnection(new URL(UCSC_URL), RANGE_HEADER);
        assertEquals(BACKUP_HOST, conn.getURL().getHost());
        conn.disconnect();

        // Clearing the substitution, as happens when it expires, returns us to the primary host
        HttpUtils.getInstance().resetUCSCBackupHost();

        conn = HttpUtils.getInstance().openConnection(new URL(UCSC_URL), RANGE_HEADER);
        assertEquals(UCSC_HOST, conn.getURL().getHost());
        conn.disconnect();
    }

    /**
     * An error response means the host answered -- it is not a reason to switch hosts.
     */
    @Test
    public void testFileNotFoundDoesNotSwitchHosts() throws Exception {

        try {
            HttpUtils.getInstance().openConnection(
                    new URL("https://" + UCSC_HOST + "/goldenPath/no_such_file_igv_test"), RANGE_HEADER);
            throw new AssertionError("Expected FileNotFoundException");
        } catch (FileNotFoundException expected) {
            // Expected
        }

        HttpURLConnection conn = HttpUtils.getInstance().openConnection(new URL(UCSC_URL), RANGE_HEADER);
        assertEquals(UCSC_HOST, conn.getURL().getHost());
        conn.disconnect();
    }
}
