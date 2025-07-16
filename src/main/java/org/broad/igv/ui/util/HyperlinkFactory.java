package org.broad.igv.ui.util;

import org.broad.igv.Globals;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.net.URI;
import java.net.URISyntaxException;

/**
 * Creates JLabels that behave like hyperlinks.
 */
public class HyperlinkFactory {

    private static Logger log = LogManager.getLogger(HyperlinkFactory.class);

    public static JLabel createLink(String label, String link) throws HeadlessException {

        JLabel hyperLink = new JLabel(label);
        hyperLink.setForeground(Globals.isDarkMode() ? Color.CYAN : Color.BLUE.darker());
        hyperLink.setCursor(new Cursor(Cursor.HAND_CURSOR));
        hyperLink.setToolTipText(link);

        try {
            final URI uri = new URI(link);

            hyperLink.addMouseListener(new MouseAdapter() {
                @Override
                public void mouseClicked(MouseEvent e) {
                    try {
                        Desktop.getDesktop().browse(uri);
                    } catch (Exception ex) {
                        log.error("Error following hyperlink: " + link, ex);
                    }
                }
            });

            return hyperLink;
        } catch (URISyntaxException e) {
            // Not a url
            return new JLabel(link);
        }
    }
}
