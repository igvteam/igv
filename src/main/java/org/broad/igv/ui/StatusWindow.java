/*
 * Created by JFormDesigner on Tue Oct 18 00:46:35 EDT 2011
 */

package org.broad.igv.ui;

import java.awt.*;
import javax.swing.*;

/**
 * @author Stan Diamond
 */
public class StatusWindow extends JFrame {
    public StatusWindow() {
        initComponents();
        editorPane1.setContentType("text/html");
    }

    public void updateText(String text) {
        this.editorPane1.setText(text);
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        scrollPane1 = new JScrollPane();
        editorPane1 = new JEditorPane();

        //======== this ========
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== scrollPane1 ========
        {
            scrollPane1.setViewportView(editorPane1);
        }
        contentPane.add(scrollPane1, BorderLayout.CENTER);
        setSize(405, 355);
        setLocationRelativeTo(getOwner());
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JScrollPane scrollPane1;
    private JEditorPane editorPane1;
    // JFormDesigner - End of variables declaration  //GEN-END:variables
}
