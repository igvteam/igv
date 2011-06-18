/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

/*
 * Created by JFormDesigner on Sat Apr 23 22:58:24 EDT 2011
 */

package org.broad.igv.peaks;

import javax.swing.event.*;
import org.broad.igv.ui.IGV;
import java.awt.*;
import java.awt.event.*;
import java.beans.*;
import javax.swing.*;
import javax.swing.border.*;
/**
 * @author Stan Diamond
 */
public class PeakControlDialog extends JDialog {



    public PeakControlDialog(Frame owner) {
        super(owner);
        initComponents();
        scoreSlider.setValue((int) PeakTrack.getScoreThreshold());
        foldChangeSlider.setValue((int) PeakTrack.getFoldChangeThreshold());
        setAlwaysOnTop(true);
        //setUndecorated(true);
        this.setLocation(owner.getBounds().width - 200, 50);
    }


    private void scoreSliderStateChanged(ChangeEvent e) {
        PeakTrack.setScoreThreshold(scoreSlider.getValue());
        IGV.getInstance().repaintDataPanels();

    }

    private void foldChangeSliderStateChanged(ChangeEvent e) {
        PeakTrack.setFoldChangeThreshold(foldChangeSlider.getValue());
        IGV.getInstance().repaintDataPanels();
    }

    private void slider1PropertyChange(PropertyChangeEvent e) {
        // TODO add your code here
    }

    private void slider1StateChanged(ChangeEvent e) {
        // TODO add your code here
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        dialogPane = new JPanel();
        contentPanel = new JPanel();
        scoreSlider = new JSlider();
        label1 = new JLabel();
        label2 = new JLabel();
        foldChangeSlider = new JSlider();

        //======== this ========
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== dialogPane ========
        {
            dialogPane.setBorder(new EmptyBorder(12, 12, 12, 12));
            dialogPane.setLayout(new BorderLayout());

            //======== contentPanel ========
            {
                contentPanel.setLayout(null);

                //---- scoreSlider ----
                scoreSlider.setPaintTicks(true);
                scoreSlider.setPaintLabels(true);
                scoreSlider.setToolTipText("Adjust score threshold");
                scoreSlider.setMajorTickSpacing(20);
                scoreSlider.addChangeListener(new ChangeListener() {
                    public void stateChanged(ChangeEvent e) {
                        scoreSliderStateChanged(e);
                    }
                });
                contentPanel.add(scoreSlider);
                scoreSlider.setBounds(5, 15, 215, scoreSlider.getPreferredSize().height);

                //---- label1 ----
                label1.setText("Score threshold:");
                contentPanel.add(label1);
                label1.setBounds(new Rectangle(new Point(5, 0), label1.getPreferredSize()));

                //---- label2 ----
                label2.setText("Fold change threshold:");
                contentPanel.add(label2);
                label2.setBounds(new Rectangle(new Point(5, 85), label2.getPreferredSize()));

                //---- foldChangeSlider ----
                foldChangeSlider.setPaintTicks(true);
                foldChangeSlider.setPaintLabels(true);
                foldChangeSlider.setToolTipText("Adjust score threshold");
                foldChangeSlider.setMajorTickSpacing(2);
                foldChangeSlider.setMinorTickSpacing(1);
                foldChangeSlider.setMaximum(10);
                foldChangeSlider.setValue(0);
                foldChangeSlider.addChangeListener(new ChangeListener() {
                    public void stateChanged(ChangeEvent e) {
                        foldChangeSliderStateChanged(e);
                    }
                });
                contentPanel.add(foldChangeSlider);
                foldChangeSlider.setBounds(5, 100, 215, 52);

                { // compute preferred size
                    Dimension preferredSize = new Dimension();
                    for(int i = 0; i < contentPanel.getComponentCount(); i++) {
                        Rectangle bounds = contentPanel.getComponent(i).getBounds();
                        preferredSize.width = Math.max(bounds.x + bounds.width, preferredSize.width);
                        preferredSize.height = Math.max(bounds.y + bounds.height, preferredSize.height);
                    }
                    Insets insets = contentPanel.getInsets();
                    preferredSize.width += insets.right;
                    preferredSize.height += insets.bottom;
                    contentPanel.setMinimumSize(preferredSize);
                    contentPanel.setPreferredSize(preferredSize);
                }
            }
            dialogPane.add(contentPanel, BorderLayout.CENTER);
        }
        contentPane.add(dialogPane, BorderLayout.CENTER);
        pack();
        setLocationRelativeTo(getOwner());
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel dialogPane;
    private JPanel contentPanel;
    private JSlider scoreSlider;
    private JLabel label1;
    private JLabel label2;
    private JSlider foldChangeSlider;
    // JFormDesigner - End of variables declaration  //GEN-END:variables
}
