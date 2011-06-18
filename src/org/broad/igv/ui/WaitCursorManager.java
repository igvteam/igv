/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
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


package org.broad.igv.ui;

//~--- JDK imports ------------------------------------------------------------

import java.awt.*;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;


/**
 * Utility class for managing IGV cursors.  The batch purpose of the class is to centrally manage
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
     *         the wait cursr.  This should be done in a finally block to insure removal.
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
        }
    }


    /**
     * A class to represent a token.
     */
    public static class CursorToken {
    }

}
