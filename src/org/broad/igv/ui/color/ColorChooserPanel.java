/*
 * Created by JFormDesigner on Sat Mar 17 09:23:03 EDT 2012
 */

package org.broad.igv.ui.color;

import java.awt.*;
import java.awt.event.*;
import java.io.Serializable;
import java.util.*;
import javax.swing.*;
import javax.swing.border.*;

import com.jidesoft.swing.*;

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
        changeColorAction(e);

    }

    private void changeColorAction(ActionEvent e) {
        JColorChooser colorChooser = new JColorChooser(selectedColor);
        JDialog dialog = JColorChooser.createDialog(this, "Dialog Title", true,
                colorChooser, null, null);
        dialog.setVisible(true);
        Color c = colorChooser.getColor();
        if (c != null) {
            setSelectedColor(c);
        }

        for (ActionListener l : listeners) {
            l.actionPerformed(e);
        }
    }

    private void colorPanelMouseClicked(MouseEvent e) {
        // TODO add your code here
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        colorPanel = new JPanel();
        changeButton = new JButton();

        //======== this ========
        setLayout(new GridBagLayout());
        ((GridBagLayout)getLayout()).columnWidths = new int[] {30, 0, 0};
        ((GridBagLayout)getLayout()).rowHeights = new int[] {0, 0};
        ((GridBagLayout)getLayout()).columnWeights = new double[] {1.0, 1.0, 1.0E-4};
        ((GridBagLayout)getLayout()).rowWeights = new double[] {1.0, 1.0E-4};

        //======== colorPanel ========
        {
            colorPanel.setPreferredSize(new Dimension(20, 10));
            colorPanel.setFocusable(false);
            colorPanel.setBorder(LineBorder.createBlackLineBorder());
            colorPanel.setMinimumSize(new Dimension(20, 10));
            colorPanel.setBackground(new Color(204, 204, 255));
            colorPanel.setLayout(null);

            { // compute preferred size
                Dimension preferredSize = new Dimension();
                for(int i = 0; i < colorPanel.getComponentCount(); i++) {
                    Rectangle bounds = colorPanel.getComponent(i).getBounds();
                    preferredSize.width = Math.max(bounds.x + bounds.width, preferredSize.width);
                    preferredSize.height = Math.max(bounds.y + bounds.height, preferredSize.height);
                }
                Insets insets = colorPanel.getInsets();
                preferredSize.width += insets.right;
                preferredSize.height += insets.bottom;
                colorPanel.setMinimumSize(preferredSize);
                colorPanel.setPreferredSize(preferredSize);
            }
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
    private JButton changeButton;
    // JFormDesigner - End of variables declaration  //GEN-END:variables
}
