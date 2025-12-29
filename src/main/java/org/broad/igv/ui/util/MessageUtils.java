/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.ui.util;

import org.broad.igv.Globals;
import org.broad.igv.logging.Level;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.ui.IGV;

import javax.swing.*;
import java.awt.*;

/**
 * Provides thread-safe, Swing-safe, utilities for interacting with JOptionPane.  Accounts for
 * (1) Swing is not thread safe => synchronize access
 * (2) JOptionPane methods must be invoked on event dispatch thread
 *
 * @author jrobinso
 */
public class MessageUtils {

    private static Logger log = LogManager.getLogger(MessageUtils.class);

    public enum InputType {INT, FLOAT, STRING}

    // Somewhat silly class, needed to pass values between threads
    static class ValueHolder {
        Object value;
    }

    /**
     * A simple data holder for returning both a string value and a checkbox state from dialogs or UI components.
     * <p>
     * This class is typically used when a dialog needs to return both a user-entered value and the state of a checkbox
     * (for example, "Do not show this message again").
     * </p>
     * <ul>
     *   <li>{@code value}: The string value entered or selected by the user. May be {@code null} if no value was provided.</li>
     *   <li>{@code isChecked}: {@code true} if the checkbox was selected, {@code false} otherwise.</li>
     * </ul>
     */
    public static class ValueCheckboxHolder {
        /**
         * The value entered or selected by the user. May be {@code null} if no value was provided.
         */
        public String value;
        /**
         * {@code true} if the checkbox was selected, {@code false} otherwise.
         */
        public boolean isChecked;
    }

    /**
     * Log the exception and show {@code message} to the user
     *
     * @param e
     * @param message
     */
    public static void showErrorMessage(String message, Exception e) {
        log.error(message, e);
        showMessage(Level.ERROR, message);
    }

    public static void showMessage(String message) {
        showMessage(null, message);
    }

    public static void showMessage(Level level, String message) {

        if (level != null) log.log(level, message);
        boolean showDialog = !(Globals.isHeadless() || Globals.isSuppressMessages() || Globals.isTesting() || Globals.isBatch());
        if (showDialog) {
            UIUtilities.invokeAndWaitOnEventThread(() -> {
                // Always use HTML for message displays, but first remove any embedded <html> tags.
                String dlgMessage = "<html>" + message.replaceAll("<html>", "");
                Frame parent = IGV.hasInstance() ? IGV.getInstance().getMainFrame() : null;

                //JEditorPane So users can select text
                JEditorPane content = new JEditorPane();
                content.setContentType("text/html");
                content.setText(dlgMessage);
                content.setEditable(false);
                Component dispMessage = content;

                //Really long messages should be scrollable
                if (dlgMessage.length() > 200) {
                    Dimension size = new Dimension(1000, content.getHeight() + 100);
                    content.setPreferredSize(size);
                    JScrollPane pane = new JScrollPane(content);
                    dispMessage = pane;
                }

                JOptionPane.showMessageDialog(parent, dispMessage);
            });
        }
    }

    public static void setStatusBarMessage(final String message) {
        log.debug("Status bar: " + message);
        if (IGV.hasInstance()) {
            IGV.getInstance().setStatusBarMessage(message);
        }
    }

    public static synchronized boolean confirm(final String message) {
        if (Globals.isHeadless()) {
            log.error("Attempted to confirm while running headless with the following message:\n" + message);
            return true;
        }

        if (Globals.isBatch()) {
            return true;
        }

        final Frame parent = IGV.hasInstance() ? IGV.getInstance().getMainFrame() : null;
        return confirm(parent, message);
    }

    /**
     * Show a yes/no confirmation dialog.
     *
     * @param component
     * @param message
     * @return
     */
    public static synchronized boolean confirm(final Component component, final String message) {


        if (Globals.isHeadless() || Globals.isBatch()) {
            return true;
        }

        if (SwingUtilities.isEventDispatchThread()) {
            int opt = JOptionPane.showConfirmDialog(component, message, "Confirm", JOptionPane.YES_NO_OPTION);
            return opt == JOptionPane.YES_OPTION;
        } else {
            final ValueHolder returnValue = new ValueHolder();
            Runnable runnable = new Runnable() {
                public void run() {
                    int opt = JOptionPane.showConfirmDialog(component, message, "Confirm", JOptionPane.YES_NO_OPTION);
                    returnValue.value = (opt == JOptionPane.YES_OPTION);
                }
            };

            UIUtilities.invokeAndWaitOnEventThread(runnable);

            return (Boolean) (returnValue.value);

        }
    }

    public static String showInputDialog(String message, String defaultValue) {

        final Frame parent = IGV.hasInstance() ? IGV.getInstance().getMainFrame() : null;
        //Pad message with spaces so it's as wide as the defaultValue
        if (defaultValue != null && message.length() < defaultValue.length()) {
            message = String.format("%-" + defaultValue.length() + "s", message);
        }
        final String actMsg = message;

        if (SwingUtilities.isEventDispatchThread()) {
            String val = JOptionPane.showInputDialog(parent, actMsg, defaultValue);
            return val;
        } else {
            final ValueHolder returnValue = new ValueHolder();
            Runnable runnable = () -> {
                String val = JOptionPane.showInputDialog(parent, actMsg, defaultValue);
                returnValue.value = val;
            };

            UIUtilities.invokeAndWaitOnEventThread(runnable);

            return (String) (returnValue.value);
        }
    }

public static ValueCheckboxHolder showInputDialog(String message, String defaultValue, String checkboxMessage) {
    final Frame parent = IGV.hasInstance() ? IGV.getInstance().getMainFrame() : null;

    if (defaultValue != null && message.length() < defaultValue.length()) {
        message = String.format("%-" + defaultValue.length() + "s", message);
    }

    JPanel panel = new JPanel(new GridLayout(0, 1));
    JTextField textField = new JTextField(defaultValue);
    panel.add(new JLabel(message));
    panel.add(textField);

    JCheckBox checkBox = null;
    if (checkboxMessage != null && !checkboxMessage.isEmpty()) {
        checkBox = new JCheckBox(checkboxMessage, true);
        panel.add(checkBox);
    }

    final JCheckBox finalCheckBox = checkBox;

    // Use an array to hold the return value, allowing modification from the lambda
    final ValueCheckboxHolder[] returnValueHolder = new ValueCheckboxHolder[1];

    Runnable dialogLogic = () -> {
        int result = JOptionPane.showConfirmDialog(parent, panel, "Enter content", JOptionPane.OK_CANCEL_OPTION, JOptionPane.PLAIN_MESSAGE);
        if (result == JOptionPane.OK_OPTION) {
            ValueCheckboxHolder holder = new ValueCheckboxHolder();
            holder.value = textField.getText();
            holder.isChecked = finalCheckBox != null && finalCheckBox.isSelected();
            returnValueHolder[0] = holder;
        } else {
            returnValueHolder[0] = null;
        }
    };

    if (SwingUtilities.isEventDispatchThread()) {
        dialogLogic.run();
    } else {
        UIUtilities.invokeAndWaitOnEventThread(dialogLogic);
    }

    return returnValueHolder[0];
}
    public static String showInputDialog(final String message) {

        final Frame parent = IGV.hasInstance() ? IGV.getInstance().getMainFrame() : null;
        if (SwingUtilities.isEventDispatchThread()) {
            String val = JOptionPane.showInputDialog(parent, message);
            return val;
        } else {
            final ValueHolder returnValue = new ValueHolder();
            Runnable runnable = new Runnable() {
                public void run() {
                    String val = JOptionPane.showInputDialog(parent, message);
                    returnValue.value = val;
                }
            };

            UIUtilities.invokeAndWaitOnEventThread(runnable);

            return (String) (returnValue.value);
        }
    }


    /**
     * Test program - call all methods from both main and swing threads
     *
     * @param args
     * @throws Exception
     */

    public static void main(String[] args) throws Exception {

        Runnable runnable = new Runnable() {
            public void run() {
                showMessage("showMessage");

                confirm("confirm");

                confirm(null, "confirm with parent");

                showInputDialog("showInputDialog", "default");

                showInputDialog("showInputDialog");
            }
        };

        // Test on main thread
        runnable.run();


        // Test on swing thread
        UIUtilities.invokeOnEventThread(runnable);

    }

}
