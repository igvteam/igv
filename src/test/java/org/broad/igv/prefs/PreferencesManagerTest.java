package org.broad.igv.prefs;

import org.broad.igv.util.TestUtils;
import org.junit.BeforeClass;
import org.junit.Ignore;
import org.junit.Test;

import static org.broad.igv.prefs.Constants.TRACK_HEIGHT_KEY;
import static org.junit.Assert.*;

/**
 * Created by jrobinso on 1/22/17.
 */
public class PreferencesManagerTest {

    @BeforeClass
    public static void setupClass() {
        PreferencesManager.setPrefsFile(TestUtils.DATA_DIR + "prefs/testUserPrefs.properties");
    }

    @Test @Ignore("Fails unless tests are run in separate JVMs")
    public void getPreferences() throws Exception {
        IGVPreferences preferences = PreferencesManager.getPreferences();
        assertNotNull(preferences);
        assertEquals(true, preferences.getAsBoolean(Constants.SAM_FLAG_UNMAPPED_PAIR));
        assertEquals(25, preferences.getAsInt(Constants.SAM_MAX_VISIBLE_RANGE));
    }

    @Test
    public void getThirdGenPreferences() throws Exception {
        IGVPreferences preferences = PreferencesManager.getPreferences(Constants.THIRD_GEN);
        assertNotNull(preferences);

        assertEquals(1000, preferences.getAsInt(Constants.SAM_MAX_VISIBLE_RANGE));
    }

    @Test
    public void getRNAPreferences() throws Exception {
        IGVPreferences preferences = PreferencesManager.getPreferences(Constants.RNA);
        assertNotNull(preferences);
        assertEquals(300, preferences.getAsInt(Constants.SAM_MAX_VISIBLE_RANGE));

        assertEquals(15, preferences.getAsInt(TRACK_HEIGHT_KEY));   // Test picking up default property from parent
    }


    @Test
    public void getThrirdGenDefaults() throws Exception {
        PreferencesManager.setPrefsFile(TestUtils.DATA_DIR + "prefs/testUserPrefs2.properties");
        IGVPreferences preferences = PreferencesManager.getPreferences(Constants.THIRD_GEN);
        assertNotNull(preferences);
        assertEquals(1, preferences.getAsInt(Constants.SAM_LARGE_INDELS_THRESHOLD));
        assertEquals(true, preferences.getAsBoolean(Constants.SAM_FLAG_LARGE_INDELS));
        assertEquals(1000, preferences.getAsInt(Constants.SAM_MAX_VISIBLE_RANGE));

        assertEquals(15, preferences.getAsInt(TRACK_HEIGHT_KEY));   // Test picking up property from parent
    }

}
