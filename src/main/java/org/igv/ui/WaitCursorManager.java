package org.igv.ui;

//~--- JDK imports ------------------------------------------------------------

import org.igv.logging.*;
import org.igv.ui.util.UIUtilities;

import javax.swing.*;
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
    private static Logger log = LogManager.getLogger(WaitCursorManager.class);

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

        UIUtilities.invokeOnEventThread(() -> IGV.getInstance().getRootPane().getGlassPane().setVisible(true));
        CursorToken token = new CursorToken();
        tokens.add(token);
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
            UIUtilities.invokeOnEventThread(() -> IGV.getInstance().getRootPane().getGlassPane().setVisible(false));
        }
    }


    /**
     * A class to represent a token.
     */
    public static class CursorToken {
    }

}
