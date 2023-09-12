package org.broad.igv.oauth;

import javax.swing.*;
import java.awt.*;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.StringSelection;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

public class OAuthURLForm extends JDialog {

    private static JPanel getPanel(String url) {
        JPanel panel = new JPanel();
        panel.setLayout(new BoxLayout(panel, BoxLayout.Y_AXIS));
        panel.setMinimumSize(new Dimension(400, 300));
        panel.setPreferredSize(new Dimension(400, 300));
        //this.add(panel);

        JPanel headerPanel = new JPanel();
        headerPanel.setLayout(new FlowLayout());
        JLabel label = new JLabel("Copy this authorization URL to your web browser");
        headerPanel.add(label);

        Button copyButton = new Button("Copy");
        copyButton.addActionListener(e -> {
            StringSelection selection = new StringSelection(url);
            Clipboard clipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
            clipboard.setContents(selection, selection);
        });
        headerPanel.add(copyButton);
        panel.add(headerPanel);

        JTextArea textArea = new JTextArea(url);
        textArea.setMargin(new Insets(10, 10, 10, 10));
        textArea.setLineWrap(true);
        panel.add(textArea);

        panel.setAlignmentX(Component.LEFT_ALIGNMENT);

        return panel;
    }

    public static void open(Frame owner, String url) {

        JPanel panel = getPanel(url);
        JOptionPane.showMessageDialog(owner,
                panel,
                "OAuth Authorization URL",
                JOptionPane.PLAIN_MESSAGE);

    }


    public static void main(String[] args) {

       open(null, "https://docs.oracle.com/javase/tutorial/displayCode.html?code=https://docs.oracle.com/javase/tutorial/uiswing/examples/layout/BoxLayoutDemoProject/src/layout/BoxLayoutDemo.java");
    }
}


