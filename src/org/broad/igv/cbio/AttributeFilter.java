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
 * Created by JFormDesigner on Fri Mar 09 15:54:00 EST 2012
 */

package org.broad.igv.cbio;

import org.broad.igv.util.StringUtils;

import javax.swing.*;
import java.awt.*;
import java.util.HashMap;
import java.util.Map;

/**
 * @author User #2
 */
public class AttributeFilter {

    private static Map<String, String> machineToHumanMap = new HashMap<String, String>(GeneNetwork.attributeMap.size());

    static {

        machineToHumanMap.put(GeneNetwork.PERCENT_MRNA_WAY_UP, "% mRNA High");
        machineToHumanMap.put(GeneNetwork.PERCENT_MRNA_WAY_DOWN, "% mRNA Low");
    }

    public static String keyToLabel(String key) {
        if (machineToHumanMap.containsKey(key)) {
            return machineToHumanMap.get(key);
        }
        String label = key.replace('_', ' ');
        label = label.replace("PERCENT", "%");
        //Looks kinda funny with 'CNA' term
        label = StringUtils.capWords(label);
        label = label.replace("Cna", "CNA");
        return label;
    }


    AttributeFilter() {
        initComponents();

        attrName.setModel(new DefaultComboBoxModel(GeneNetwork.attributeMap.keySet().toArray()));
        attrName.insertItemAt(GeneNetwork.PERCENT_ALTERED, 0);
        attrName.setSelectedIndex(0);

        attrName.setRenderer(new ListCellRenderer() {
            @Override
            public Component getListCellRendererComponent(JList list, Object value, int index, boolean isSelected, boolean cellHasFocus) {
                String label = AttributeFilter.keyToLabel("" + value);
                JLabel comp = new JLabel(label);
                return comp;
            }
        });
    }

    JPanel getPanel() {
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
        addRow = new JButton();

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
                minVal.setMaximumSize(new Dimension(80, 28));
                minVal.setPreferredSize(new Dimension(50, 28));
                minVal.setMinimumSize(new Dimension(50, 28));
                filterRow.add(minVal);

                //---- label2 ----
                label2.setText("Max");
                filterRow.add(label2);

                //---- maxVal ----
                maxVal.setText("100.0");
                maxVal.setMaximumSize(new Dimension(80, 28));
                maxVal.setPreferredSize(new Dimension(50, 28));
                maxVal.setMinimumSize(new Dimension(50, 28));
                filterRow.add(maxVal);

                //---- delRow ----
                delRow.setText("-");
                delRow.setMaximumSize(new Dimension(20, 29));
                delRow.setMinimumSize(new Dimension(20, 29));
                delRow.setPreferredSize(new Dimension(20, 29));
                delRow.setToolTipText("Delete this filter");
                delRow.setMargin(new Insets(2, 2, 2, 2));
                filterRow.add(delRow);

                //---- addRow ----
                addRow.setText("+");
                addRow.setMaximumSize(new Dimension(20, 29));
                addRow.setMinimumSize(new Dimension(20, 29));
                addRow.setPreferredSize(new Dimension(20, 29));
                addRow.setToolTipText("Add a new filter");
                addRow.setMargin(new Insets(2, 2, 2, 2));
                filterRow.add(addRow);
            }
            panel1.add(filterRow, BorderLayout.CENTER);
        }
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel panel1;
    private JPanel filterRow;
    private JComboBox attrName;
    private JLabel label1;
    JTextField minVal;
    private JLabel label2;
    JTextField maxVal;
    private JButton delRow;
    private JButton addRow;
    // JFormDesigner - End of variables declaration  //GEN-END:variables

    JButton getDelRow() {
        return delRow;
    }

    JButton getAddRow() {
        return addRow;
    }

    void setShowDel(boolean showDel) {
        delRow.setVisible(showDel);
    }

    void setIsLast(boolean isLast) {
        addRow.setVisible(isLast);
    }

    JComboBox getAttrName() {
        return attrName;
    }

}
