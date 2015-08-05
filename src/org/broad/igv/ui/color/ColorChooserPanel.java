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

/*
 * Created by JFormDesigner on Sat Mar 17 09:23:03 EDT 2012
 */

package org.broad.igv.ui.color;

import javax.swing.*;
import javax.swing.border.LineBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.io.Serializable;
import java.util.ArrayList;

/**
 * @author Jim Robinson
 */
public class ColorChooserPanel extends JPanel implements Serializable {

    Color selectedColor;

    java.util.List<ActionListener> listeners = new ArrayList<ActionListener>();


    public ColorChooserPanel() {
        this(Color.white);
    }

    public ColorChooserPanel(Color selectedColor) {
        this.selectedColor = selectedColor;
        initComponents();
        setSelectedColor(selectedColor);
    }

    public void setSelectedColor(Color selectedColor) {
        this.selectedColor = selectedColor;
        colorPanel.setBackground(selectedColor);
    }


    public Color getSelectedColor() {
        return selectedColor;
    }

    public void addActionListener(ActionListener listener) {
        listeners.add(listener);
    }


    private void changeButtonActionPerformed(ActionEvent e) {
        changeColorAction();

    }

    private void colorPanelMouseClicked(MouseEvent e) {
        changeColorAction();
    }

    private void changeColorAction() {
        JColorChooser colorChooser = new JColorChooser(selectedColor);
        JDialog dialog = JColorChooser.createDialog(this, "Dialog Title", true,
                colorChooser, null, null);
        dialog.setVisible(true);
        Color c = colorChooser.getColor();
        if (c != null) {
            setSelectedColor(c);
        }

        ActionEvent evt = new ActionEvent(this, 0, "");
        for (ActionListener l : listeners) {
            l.actionPerformed(evt);
        }
    }

    private void initComponents() {

        colorPanel = new JPanel();
        changeButton = new JButton();
        //Workaround of Metal LAF bug, see http://bugs.openjdk.java.net/show_bug.cgi?id=100280
        if (UIManager.getLookAndFeel().getName().toLowerCase().contains("metal")) {
            changeButton = new JMenuItem();
        }

        //======== this ========
        setLayout(new GridBagLayout());
        ((GridBagLayout) getLayout()).columnWidths = new int[]{30, 0, 0};
        ((GridBagLayout) getLayout()).rowHeights = new int[]{0, 0};
        ((GridBagLayout) getLayout()).columnWeights = new double[]{1.0, 1.0, 1.0E-4};
        ((GridBagLayout) getLayout()).rowWeights = new double[]{1.0, 1.0E-4};

        //======== colorPanel ========
        {
            colorPanel.setPreferredSize(new Dimension(20, 15));
            colorPanel.setMinimumSize(new Dimension(20, 5));
            colorPanel.setMaximumSize(new Dimension(20, 15));
            colorPanel.setFocusable(false);
            colorPanel.setBorder(LineBorder.createBlackLineBorder());
            colorPanel.setBackground(new Color(204, 204, 255));
            colorPanel.setLayout(null);
            colorPanel.addMouseListener(new MouseAdapter() {
                @Override
                public void mouseClicked(MouseEvent mouseEvent) {
                    colorPanelMouseClicked(mouseEvent);
                }
            });

        }
        add(colorPanel, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0,
                GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                new Insets(0, 0, 0, 0), 0, 0));

        //---- changeButton ----
        changeButton.setHorizontalAlignment(SwingConstants.LEFT);
        changeButton.setIcon(UIManager.getIcon("Menu.arrowIcon"));
        changeButton.setToolTipText("Change color");
        changeButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                changeButtonActionPerformed(e);
            }
        });
        add(changeButton, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0,
                GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                new Insets(0, 0, 0, 0), 0, 0));
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel colorPanel;
    private AbstractButton changeButton;
    // JFormDesigner - End of variables declaration  //GEN-END:variables
}
