/*
 * Created by JFormDesigner on Fri Mar 09 15:54:00 EST 2012
 */

package org.broad.igv.cbio;

import javax.swing.*;
import java.awt.*;

/**
 * @author User #2
 */
public class AttributeFilter {

    public AttributeFilter() {
        initComponents();
        attrName.setModel(new DefaultComboBoxModel(GeneNetwork.attribute_map.keySet().toArray()));
        attrName.addItem(GeneNetwork.PERCENT_ALTERED);
    }

    public Component getComponent() {
        return this.filterRow;
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        panel1 = new JPanel();
        filterRow = new JPanel();
        attrName = new JComboBox();
        label1 = new JLabel();
        minVal = new JTextField();
        label2 = new JLabel();
        maxVal = new JTextField();
        delRow = new JButton();

        //======== panel1 ========
        {
            panel1.setLayout(new BorderLayout());

            //======== filterRow ========
            {
                filterRow.setMaximumSize(new Dimension(100000, 29));
                filterRow.setPreferredSize(new Dimension(600, 29));
                filterRow.setLayout(new BoxLayout(filterRow, BoxLayout.X_AXIS));

                //---- attrName ----
                attrName.setMaximumRowCount(12);
                attrName.setMaximumSize(new Dimension(300, 28));
                attrName.setToolTipText("Attribute by which to filter");
                attrName.setAlignmentX(0.0F);
                filterRow.add(attrName);

                //---- label1 ----
                label1.setText("Min");
                filterRow.add(label1);

                //---- minVal ----
                minVal.setText("0.0");
                minVal.setMaximumSize(new Dimension(100, 28));
                minVal.setPreferredSize(new Dimension(50, 28));
                minVal.setMinimumSize(new Dimension(50, 28));
                filterRow.add(minVal);

                //---- label2 ----
                label2.setText("Max");
                filterRow.add(label2);

                //---- maxVal ----
                maxVal.setText("100.0");
                maxVal.setMaximumSize(new Dimension(100, 28));
                maxVal.setPreferredSize(new Dimension(50, 28));
                maxVal.setMinimumSize(new Dimension(50, 28));
                filterRow.add(maxVal);

                //---- delRow ----
                delRow.setText("-");
                delRow.setMaximumSize(new Dimension(20, 29));
                delRow.setMinimumSize(new Dimension(20, 29));
                delRow.setPreferredSize(new Dimension(20, 29));
                filterRow.add(delRow);
            }
            panel1.add(filterRow, BorderLayout.CENTER);
        }
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel panel1;
    private JPanel filterRow;
    JComboBox attrName;
    private JLabel label1;
    JTextField minVal;
    private JLabel label2;
    JTextField maxVal;
    private JButton delRow;
    // JFormDesigner - End of variables declaration  //GEN-END:variables

    public JButton getDelRow() {
        return delRow;
    }
}
