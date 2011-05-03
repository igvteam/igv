/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.ui;

//~--- JDK imports ------------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.track.TrackType;

import java.awt.*;

/**
 * @author jrobinso
 */
public class UIConstants {

    private static Logger log = Logger.getLogger(UIConstants.class);
    /**
     * Field description
     */
    final public static int groupGap = 10;
    /**
     * Field description
     */
    final public static String APPLICATION_NAME = "IGV";
    /**
     * Field description
     */
    final public static String APPLICATION_LONG_NAME = "Integrative Genomics Viewer";
    /**
     * Field description
     */
    final public static Dimension preferredSize = new Dimension(1000, 750);
    // To support mutation track overlay.  Generalize later.  Cancer specific option.
    public static TrackType overlayTrackType = TrackType.MUTATION;
    // Default user folder


    private static int doubleClickInterval = -1;
    // General
    final static public String CLICK_ITEM_TO_EDIT_TOOLTIP = "Click this item bring up its editor";// Toolbar Menu Item Tooltips
    final static public String JUMP_TO_WHOLE_GENOME_VIEW_TOOLTIP = "Jump to whole genome view";
    final static public String JUMP_TO_LOCUS_TOOLTIP = "Jump to Gene or Locus";
    final static public String SELECT_CHROMOSOME_TOOLTIP = "Select a chromosome to view";
    final static public String ZOOM_TOOL_TOOLTIP = "Click + to zoom in - Click - to zoom out";// File Menu Item Tooltips
    final static public String LOAD_TRACKS_TOOLTIP = "Load data, features or sample information";
    final static public String LOAD_SERVER_DATA_TOOLTIP = "Load data, features or sample information from a server";
    final static public String LOAD_ATTRIBUTES_TOOLTIP = "Load track attributes";
    final static public String SAVE_IMAGE_TOOLTIP = "Capture and save an image";
    final static public String NEW_SESSION_TOOLTIP = "Create a new session";
    final static public String SAVE_SESSION_TOOLTIP = "Save the current session";
    final static public String SAVE_SESSION_AS_TOOLTIP = "Save the current session to the specified file";
    final static public String RESTORE_SESSION_TOOLTIP = "Reload the named session";
    final static public String EXIT_TOOLTIP = "Exit the application";
    final static public String EXPORT_REGION_TOOLTIP = "Allows currently selected regions to be exported to a file";
    final static public String IMPORT_REGION_TOOLTIP = "Allows previously exported regions to be reloaded";
    final static public String CLEAR_REGION_TOOLTIP = "Clear all regions of interest";// Edit Menu Item Tooltips
    final static public String CHANGE_GENOME_TOOLTIP = "Switch the current genome";
    final static public String IMPORT_GENOME_TOOLTIP =
            "Create a user-defined genome and makes it available for use in the application";
    final static public String LOAD_GENOME_TOOLTIP =
            "Make an externally defined genome available for use in the application";
    final static public String REMOVE_USER_DEFINE_GENOME_TOOLTIP =
            "Removes user-defined genomes from the drop-down list";
    final static public String CLEAR_GENOME_CACHE_TOOLTIP =
            "Clears locally cached versions of IGV hosted genomes.";// IGVPanel Menu Item Tooltips
    final static public String PREFERENCE_TOOLTIP =
            "Set user specific preferences";
    final static public String SHOW_ATTRIBUTE_DISPLAY_TOOLTIP =
            "Show or hide the attribute display";
    final static public String ENABLE_REGIONS_OF_INTEREST_TOOLTIP =
            "Enable the \"Region of Interest\" tool";
    final static public String REFRESH_TOOLTIP =
            "Refresh the application's display";
    final static public String SHOW_HEATMAP_LEGEND_TOOLTIP =
            "IGVPanel or edit color legends and scales";
    final static public String SELECT_DISPLAYABLE_ATTRIBUTES_TOOLTIP =
            "Customize attribute display to show only checked attributes";
    final static public String USE_IMAGE_CACHING_TOOLTIP =
            "Use image caching to improve display performance";
    final static public String DIRECT_DRAW_DISABLED_TOOLTIP =
            "Checked this item to prevent the use of direct draw" +
                    " when rendering images - Uncheck it to use the default system behavior";// Tracks Menu Item Tooltips
    final static public String SORT_TRACKS_TOOLTIP =
            "Sort tracks by attribute value";
    final static public String GROUP_TRACKS_TOOLTIP =
            "Group tracks";
    final static public String FILTER_TRACKS_TOOLTIP =
            "Filter tracks by attribute value";
    final static public String RESET_DEFAULT_TRACK_HEIGHT_TOOLTIP =
            "Reset all track to the default track height";
    final static public String SET_DEFAULT_TRACK_HEIGHT_TOOLTIP =
            "Set the height for all tracks";
    final static public String FIT_DATA_TO_WINDOW_TOOLTIP =
            "Resizes all track heights in order to make all tracks visible in " +
                    "their display with no vertical scrolling";// Help Menu Item Tooltips
    final static public String HELP_TOOLTIP =
            "Open web help page";
    final static public String TUTORIAL_TOOLTIP =
            "Open tutorial web page";
    final static public String ABOUT_TOOLTIP =
            "Display application information";// DataPanelTool Menu
    final static public String MACRO_SNAPSHOTS =
            "Macro Snapshots";
    final static public String RESET_FACTORY_TOOLTIP =
            "Restores all user preferences to their default settings.";
    public static final String NAVIGATE_REGION_TOOLTIP = "Navigate regions";

    public static Font boldFont = FontManager.getScalableFont(Font.BOLD, 12);

    public static int getDoubleClickInterval() {

        if (doubleClickInterval < 0) {
            try {
                Number obj = (Number) Toolkit.getDefaultToolkit().getDesktopProperty("awt.multiClickInterval");
                doubleClickInterval = obj.intValue();
            } catch (Exception e) {
                log.info("Error retrieving doubleClickInterval", e);
                doubleClickInterval = 500;
            }

        }
        return doubleClickInterval;
    }




    /**
     * Field description
     */
    final public static String OVERWRITE_SESSION_MESSAGE =
            "<html>Opening a session will unload all current data. " + "<br>Are you sure you wish to continue?";
    /**
     * Field description
     */
    final public static String NEW_SESSION_MESSAGE =
            "<html>Creating a new session will unload all current data. " + "<br>Are you sure you wish to continue?";
    final public static String DEFAULT_SERVER_GENOME_ARCHIVE_LIST = "http://igv.broadinstitute.org/genomes/genomes.txt";
    final public static String CANNOT_ACCESS_SERVER_GENOME_LIST = "The Genome server is currently inaccessible.";
    final public static String INVALID_SERVER_GENOME_LIST_HEADER = "Genomes cannot be retrieved from the server. " + "The server-side genome list is invalid!";
    // Session Folder
    /**
     * Field description
     */
    final public static String SESSION_FOLDER = "/igv_sessions";
    /**
     * Field description
     */
    final public static int NUMBER_OF_RECENT_SESSIONS_TO_LIST = 3;
    /**
     * Field description
     */
    final public static String DEFAULT_SESSION_FILE = "igv_session" + Globals.SESSION_FILE_EXTENSION;
    // URLs
    /**
     * Field description
     */
    final static public String SERVER_BASE_URL = "http://www.broadinstitute.org/";
    /**
     * Field description
     */
    final public static String IGV_LOG_SERVER_URL = SERVER_BASE_URL + "igv/LogServlet";
    // Colors
    // float[] hsb = Color.RGBtoHSB(255, 255, 210, null);

    final static public Color LIGHT_YELLOW = new Color(255, 244, 201);

    final public static Color ZOOMED_OUT_COLOR = new Color(238, 239, 240);

    final public static Color TRACK_BORDER_GRAY = new Color(200, 200, 210);

    public static Color NO_DATA_COLOR = new Color(200, 200, 200, 150);
    // GENOME
    /**
     * Field description
     */
    final static public String IMPORT_GENOME_LIST_MENU_ITEM = "Import Genome...";
    /**
     * Field description
     */
    final static public String LOAD_GENOME_LIST_MENU_ITEM = "Load Genome...";
    /**
     * Field description
     */
    final static public String REMOVE_GENOME_LIST_MENU_ITEM = "Remove Imported Genomes...";
    /**
     * Field description
     */
    final static public String GENOME_LIST_SEPARATOR = "--SEPARATOR--";

}
