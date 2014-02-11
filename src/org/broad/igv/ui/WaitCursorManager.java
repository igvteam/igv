/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */


package org.broad.igv.ui;

//~--- JDK imports ------------------------------------------------------------

import java.awt.*;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;


/**
 * Utility class for managing IGV cursors.  The main purpose of the class is to centrally manage
 * a global wait cursor.  When in "wait" mode component set cursor events are ignored, or rather
 * saved in a cached until the wait cursor is removed.
 *
 * @author jrobinso
 */
public class WaitCursorManager {


    /**
     * A set of tokens, one for each call to "showWaitCursor".  These are removed in the
     * "removeWaitCursor" method.  The presence of a token in this list indicates that IGV is
     * in the wait state.
     */
    static Set<CursorToken> tokens = Collections.synchronizedSet(new HashSet());

    /**
     * The wait cursor, defined statically for convenience.
     */
    static Cursor waitCursor = Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR);

    /**
     * Show the wait cursor on all components.  Add a token to represent this invocation of
     * showWaitCursor
     *
     * @return token representing this invocation.  This token should be used by clients to remove
     *         the wait cursor.  This should be done in a finally block to insure removal.
     */
    public static CursorToken showWaitCursor() {
        IGV.getRootPane().getGlassPane().setVisible(true);
        CursorToken token = new CursorToken();
        tokens.add(token);
        // Return a token representing this wait cursor set.  The token is used to release the
        // wait cursor.
        return token;
    }

    /**
     * Remove the token for a showWaitCursor() invocation.  This indicates that the client has completed
     * its task and removed the wait cursor request.  If the last token has been removed reset
     * the cursors on the components to their last requested cursor, or the default cursor if
     * there are no outstanding requests.
     *
     * @param token
     */
    public static void removeWaitCursor(CursorToken token) {
        tokens.remove(token);
        if (tokens.isEmpty()) {
            IGV.getRootPane().getGlassPane().setVisible(false);
            IGV.getInstance().getContentPane().getStatusBar().deactivateCancelButton();
        }
    }


    /**
     * A class to represent a token.
     */
    public static class CursorToken {
    }

}
