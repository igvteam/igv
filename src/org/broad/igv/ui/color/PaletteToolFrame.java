/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

/*
 * Created by JFormDesigner on Wed Nov 02 23:02:37 EDT 2011
 */

package org.broad.igv.ui.color;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;

/**
 * @author Stan Diamond
 */
public class PaletteToolFrame extends JFrame {
    public PaletteToolFrame() {
        initComponents();
        desaturateCheckbox.setSelected(colorPanel.showGrayScale);
    }

    private void saveItemActionPerformed(ActionEvent e) {
        java.util.List<ColorPanel.Palette> paletteList = colorPanel.paletteList;
        if (paletteList != null) {
            java.awt.FileDialog fd = new FileDialog(this);
            fd.setMode(FileDialog.SAVE);
            fd.setVisible(true);
            String f = fd.getFile();
            if (f != null) {
                try {
                    File file = new File(fd.getDirectory(), f);
                    PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(file)));
                    for (ColorPanel.Palette p : paletteList) {
                        pw.println(p.label);
                        for (ColorPanel.Swatch s : p.swatches) {
                            pw.println(ColorUtilities.colorToString(s.color));
                        }
                    }
                } catch (Exception e1) {
                    e1.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                }
            }
        }
    }

    private void showGrayScaleActionPerformed(ActionEvent e) {
        colorPanel.showGrayScale = desaturateCheckbox.isSelected();
        colorPanel.repaint();
        ;
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        menuBar1 = new JMenuBar();
        menu1 = new JMenu();
        saveItem = new JMenuItem();
        hSpacer1 = new JPanel(null);
        desaturateCheckbox = new JCheckBox();
        colorPanel = new ColorPanel();

        //======== this ========
        setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== menuBar1 ========
        {

            //======== menu1 ========
            {
                menu1.setText("File");

                //---- saveItem ----
                saveItem.setText("Save...");
                saveItem.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        saveItemActionPerformed(e);
                    }
                });
                menu1.add(saveItem);
            }
            menuBar1.add(menu1);
            menuBar1.add(hSpacer1);

            //---- desaturateCheckbox ----
            desaturateCheckbox.setText("De-saturate");
            desaturateCheckbox.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    showGrayScaleActionPerformed(e);
                }
            });
            menuBar1.add(desaturateCheckbox);
        }
        setJMenuBar(menuBar1);
        contentPane.add(colorPanel, BorderLayout.CENTER);
        setSize(890, 570);
        setLocationRelativeTo(getOwner());
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JMenuBar menuBar1;
    private JMenu menu1;
    private JMenuItem saveItem;
    private JPanel hSpacer1;
    private JCheckBox desaturateCheckbox;
    private ColorPanel colorPanel;
    // JFormDesigner - End of variables declaration  //GEN-END:variables

    public static void main(String[] args) {
        (new PaletteToolFrame()).setVisible(true);
    }
}
