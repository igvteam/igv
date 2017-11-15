package org.broad.igv.ui.javafx;

import javafx.scene.text.Font;
import javafx.scene.text.Text;

/**
 * A substitute for the old Swing FontMetrics class since JavaFX doesn't have a (public) version.
 * This only implements the method(s) we need, so it's not general purpose.
 *
 * @author eby
 */
public class FontMetrics {

    // From https://stackoverflow.com/questions/13015698/how-to-calculate-the-pixel-width-of-a-string-in-javafx

    /**
     * Get the width of a text String with the given Font applied.
     * This is taken from https://stackoverflow.com/questions/13015698/how-to-calculate-the-pixel-width-of-a-string-in-javafx
     * Note: for now, this code takes the simple form of working with the font ONLY.  It may need to be more complex for
     * dealing with CSS.
     *
     * @param textString
     * @param font
     * @return
     */
    public static double getTextWidthInFont(String textString, Font font) {
        // We might consider caching these Text objects if it's too expensive to create them. 
        Text text = new Text(textString);
        text.setFont(font);
        return text.getBoundsInLocal().getWidth();
    }

    /**
     * Get the width of a text String with the Font of the given Text instance.  This is an optimization for where it's
     * necessary to calculate widths for many Strings.
     *
     * @param textString
     * @param text
     * @return
     */
    public static double getTextWidthInFont(String textString, Text text) {
        text.setText(textString);
        return text.getBoundsInLocal().getWidth();
    }

}
