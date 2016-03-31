/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */


package org.broad.igv.ui;

//~--- JDK imports ------------------------------------------------------------

import org.broad.igv.ui.util.UIUtilities;

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
     * Show the wait cursor on all components.  Add a token to represent this invocation of
     * showWaitCursor
     *
     * @return token representing this invocation.  This token should be used by clients to remove
     * the wait cursor.  This should be done in a finally block to insure removal.
     */
    public static CursorToken showWaitCursor() {

        IGV.getRootPane().getGlassPane().setVisible(true);
        CursorToken token = new CursorToken();
        tokens.add(token);

        //if(tokens.size() == 1) UIUtilities.activateMainFrame();   // neccessary for busy cursor
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
