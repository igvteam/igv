package org.broad.igv.ui;

import org.broad.igv.util.StringUtils;

import javax.swing.*;
import java.awt.*;

/**
 * @author jrobinso
 * @date Mar 2, 2011
 */
public class UserPasswordPane {

    JPanel passPanel;
    JTextField userField;
    JPasswordField passwordField;
    String resourceString;

    public UserPasswordPane(String resourceString) {
        this.resourceString = resourceString;
        init();
    }

    private void init() {
        passPanel = new JPanel();
        passPanel.setLayout(new GridLayout(6, 1));

        JLabel message = new JLabel("Please enter your username and password for:");
        JLabel location = new JLabel(StringUtils.checkLength(resourceString, 80));
        JLabel username = new JLabel("User:");
        JLabel password = new JLabel("Pass:");
        userField = new JTextField();
        passwordField = new JPasswordField();
        passPanel.add(message);
        passPanel.add(location);
        passPanel.add(username);
        passPanel.add(userField);
        passPanel.add(password);
        passPanel.add(passwordField);
    }


    public boolean showDialog() {
        int a = JOptionPane.showConfirmDialog(IGV.getInstance().getMainFrame(), passPanel, "Authentication Required", JOptionPane.OK_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE);

        if (a == JOptionPane.CANCEL_OPTION) {
            return false;
        } else {
            return true;
        }
    }

    public String getUser() {
        return userField.getText();
    }

    public String getPassword() {
        return new String(passwordField.getPassword());
    }
}
