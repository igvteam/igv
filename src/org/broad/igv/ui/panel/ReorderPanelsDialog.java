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
 * Created by JFormDesigner on Mon May 16 00:05:02 EDT 2011
 */

package org.broad.igv.ui.panel;

import java.awt.*;
import java.awt.event.*;
import java.util.ArrayList;
import javax.swing.*;
import javax.swing.border.*;

import org.broad.igv.track.Track;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.*;

/**
 * @author Stan Diamond
 */
public class ReorderPanelsDialog extends JDialog {

    private boolean canceled = false;

    public ReorderPanelsDialog(Frame owner) {
        super(owner);
        initComponents();
        init();
    }


    public void init() {
        java.util.List<TrackPanel> panes = IGV.getInstance().getMainPanel().getTrackPanels();
        java.util.List<Wrapper> wrappers = new ArrayList(panes.size());
        for(TrackPanel pane : panes) {
           wrappers.add(new Wrapper(pane));
        }
        list.setElements(wrappers);
    }

    private void okButtonActionPerformed(ActionEvent e) {
        setVisible(false);
        java.util.List<String> orderedNames = new ArrayList();
        for(Object obj : list.getElements()) {
            Wrapper w = (Wrapper) obj;
            orderedNames.add(w.panelName);
        }
        IGV.getInstance().getMainPanel().reorderPanels(orderedNames);
    }

    private void cancelButtonActionPerformed(ActionEvent e) {
        canceled = true;
        setVisible(false);
    }

    public boolean isCanceled() {
        return canceled;
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        dialogPane = new JPanel();
        contentPanel = new JPanel();
        scrollPane1 = new JScrollPane();
        list = new ReorderableJList();
        label1 = new JLabel();
        buttonBar = new JPanel();
        okButton = new JButton();
        cancelButton = new JButton();

        //======== this ========
        setTitle("Reorder Panels");
        setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== dialogPane ========
        {
            dialogPane.setBorder(new EmptyBorder(12, 12, 12, 12));
            dialogPane.setLayout(new BorderLayout());

            //======== contentPanel ========
            {
                contentPanel.setLayout(new BorderLayout());

                //======== scrollPane1 ========
                {
                    scrollPane1.setViewportView(list);
                }
                contentPanel.add(scrollPane1, BorderLayout.CENTER);

                //---- label1 ----
                label1.setText("To reorder panels drag and drop entries in the list below.");
                contentPanel.add(label1, BorderLayout.NORTH);
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
                    public void actionPerformed(ActionEvent e) {
                        okButtonActionPerformed(e);
                    }
                });
                buttonBar.add(okButton, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                    new Insets(0, 0, 0, 5), 0, 0));

                //---- cancelButton ----
                cancelButton.setText("Cancel");
                cancelButton.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
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
        pack();
        setLocationRelativeTo(getOwner());
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    static class Wrapper {
        String panelName;
        String printName;
        Wrapper(TrackPanel panel) {
            this.panelName = panel.getName();

            StringBuffer buffer = new StringBuffer(50);
            int nChars = 0;
            for(Track track : panel.getTracks()) {
                if(nChars > 0) {
                    buffer.append(", ");
                }
                nChars += track.getName().length();
                buffer.append(track.getName());
            }
            this.printName = buffer.toString();

        }

        public String toString() {
            return printName;
        }
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel dialogPane;
    private JPanel contentPanel;
    private JScrollPane scrollPane1;
    private ReorderableJList list;
    private JLabel label1;
    private JPanel buttonBar;
    private JButton okButton;
    private JButton cancelButton;
    // JFormDesigner - End of variables declaration  //GEN-END:variables
}
