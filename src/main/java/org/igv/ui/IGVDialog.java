package org.igv.ui;

import javax.swing.*;
import java.awt.*;

/**
 * Extension of JDialog that disables IGV hot key dispatching while dialogs are open.
 */
public class IGVDialog extends JDialog {

    public IGVDialog() {
    }

    public IGVDialog(Dialog owner) {
        super(owner);
    }

    public IGVDialog(Dialog owner, boolean modal) {
        super(owner, modal);
    }

    public IGVDialog(Dialog owner, String title) {
        super(owner, title);
    }

    public IGVDialog(Frame owner) {
        super(owner);
    }

    public IGVDialog(Frame owner, boolean modal) {
        super(owner, modal);
    }

    public IGVDialog(Frame owner, String title) {
        super(owner, title);
    }

    public IGVDialog(Frame owner, String title, boolean modal) {
        super(owner, title, modal);
    }

}
