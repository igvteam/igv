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
 * Created by JFormDesigner on Sat Feb 08 20:06:23 EST 2014
 */

package org.broad.igv.cursor;

import com.jidesoft.swing.JideBoxLayout;
import com.jidesoft.swing.JideButton;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.border.*;

/**
 *
 */
public class CursorFilterDialog extends JDialog {

    java.util.List<CursorTrack> tracks;
    RegionFilter filter;
    boolean canceled;

    String[] trackNames;

    public static void main(String[] args) {

        (new CursorFilterDialog(null)).setVisible(true);
    }

    // Constructor for testing
    public CursorFilterDialog(Frame owner) {
        this(owner, null, new RegionFilter());
    }

    public CursorFilterDialog(Frame owner, java.util.List<CursorTrack> tracks, RegionFilter filter) {
        super(owner);
        initComponents();

        this.tracks = tracks;
        if (tracks == null || tracks.isEmpty()) {
            this.trackNames = new String[]{};
        } else {
            this.trackNames = new String[tracks.size()];
            for (int i = 0; i < tracks.size(); i++) {
                trackNames[i] = tracks.get(i).getName();
            }
        }

        if (filter == null) {
            filter = new RegionFilter();
        }
        java.util.List<RegionFilter.Clause> clauses = filter.getClauses();
        if (clauses == null || clauses.isEmpty()) {
            addFilter(null);
        } else {
            for (RegionFilter.Clause c : clauses) {
                addFilter(c);
            }
        }
    }

    private void okButtonActionPerformed(ActionEvent e) {

        // TODO -- create RegionFilter and update model
        canceled = false;
        setVisible(false);
    }

    private void cancelButtonActionPerformed(ActionEvent e) {
        canceled = true;
        setVisible(false);
    }

    public boolean isCanceled() {
        return canceled;
    }

    private void addFilter(RegionFilter.Clause clause) {

        final JPanel p = new JPanel();
        p.setMaximumSize(new Dimension(1000, 30));
        p.setLayout(new JideBoxLayout(p, BoxLayout.X_AXIS));

        JComboBox trackCB = new JComboBox(trackNames);
        p.add(trackCB, JideBoxLayout.VARY);

        JComboBox signalCB = new JComboBox(new String[] {"Score", "Signal"});
        p.add(signalCB, JideBoxLayout.FIX);

        JComboBox predCB = new JComboBox(new String[]{"is greater than", "is less than"});
        predCB.setSize(new Dimension(100, 20));
        p.add(predCB, JideBoxLayout.FLEXIBLE);

        JTextField valueField = new JTextField();
        valueField.setSize(new Dimension(50, 20));
        valueField.setPreferredSize(new Dimension(50, 20));
        valueField.setMaximumSize(new Dimension(50, 20));
        p.add(valueField, JideBoxLayout.FIX);

        JideButton plusButton = new JideButton("+");
        plusButton.setSize(new Dimension(20, 20));
        plusButton.setPreferredSize(new Dimension(20, 20));
        plusButton.setMaximumSize(new Dimension(20, 20));
        plusButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                addFilter(null);
            }
        });
        p.add(plusButton, JideBoxLayout.FIX);

        JideButton minusButton = new JideButton("-");
        minusButton.setSize(new Dimension(20, 20));
        minusButton.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent e) {
                System.out.println("Remove " + p);
            }
        });
        p.add(minusButton, JideBoxLayout.FIX);


        this.filterPanel.add(p);
        this.filterPanel.revalidate();

    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        dialogPane = new JPanel();
        contentPanel = new JPanel();
        panel1 = new JPanel();
        matchAllRB = new JRadioButton();
        matchAnyRB = new JRadioButton();
        filterPanel = new JPanel();
        buttonBar = new JPanel();
        okButton = new JButton();
        cancelButton = new JButton();

        //======== this ========
        setModal(true);
        setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
        setTitle("Define Regions");
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== dialogPane ========
        {
            dialogPane.setBorder(new EmptyBorder(12, 12, 12, 12));
            dialogPane.setLayout(new BorderLayout());

            //======== contentPanel ========
            {
                contentPanel.setLayout(new BorderLayout(0, 5));

                //======== panel1 ========
                {
                    panel1.setLayout(new FlowLayout(FlowLayout.LEFT));

                    //---- matchAllRB ----
                    matchAllRB.setText("Match all of the following");
                    matchAllRB.setSelected(true);
                    panel1.add(matchAllRB);

                    //---- matchAnyRB ----
                    matchAnyRB.setText("Match any of the following");
                    panel1.add(matchAnyRB);
                }
                contentPanel.add(panel1, BorderLayout.NORTH);

                //======== filterPanel ========
                {
                    filterPanel.setLayout(new BoxLayout(filterPanel, BoxLayout.Y_AXIS));
                }
                contentPanel.add(filterPanel, BorderLayout.CENTER);
            }
            dialogPane.add(contentPanel, BorderLayout.CENTER);

            //======== buttonBar ========
            {
                buttonBar.setBorder(new EmptyBorder(12, 0, 0, 0));
                buttonBar.setLayout(new GridBagLayout());
                ((GridBagLayout)buttonBar.getLayout()).columnWidths = new int[] {0, 85, 80};
                ((GridBagLayout)buttonBar.getLayout()).columnWeights = new double[] {1.0, 0.0, 0.0};

                //---- okButton ----
                okButton.setText("OK");
                okButton.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        okButtonActionPerformed(e);
                        okButtonActionPerformed(e);
                    }
                });
                buttonBar.add(okButton, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                    new Insets(0, 0, 0, 5), 0, 0));

                //---- cancelButton ----
                cancelButton.setText("Cancel");
                cancelButton.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        cancelButtonActionPerformed(e);
                        cancelButtonActionPerformed(e);
                    }
                });
                buttonBar.add(cancelButton, new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                    new Insets(0, 0, 0, 0), 0, 0));
            }
            dialogPane.add(buttonBar, BorderLayout.SOUTH);
        }
        contentPane.add(dialogPane, BorderLayout.CENTER);
        setSize(765, 570);
        setLocationRelativeTo(getOwner());

        //---- buttonGroup1 ----
        ButtonGroup buttonGroup1 = new ButtonGroup();
        buttonGroup1.add(matchAllRB);
        buttonGroup1.add(matchAnyRB);
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel dialogPane;
    private JPanel contentPanel;
    private JPanel panel1;
    private JRadioButton matchAllRB;
    private JRadioButton matchAnyRB;
    private JPanel filterPanel;
    private JPanel buttonBar;
    private JButton okButton;
    private JButton cancelButton;
    // JFormDesigner - End of variables declaration  //GEN-END:variables
}
