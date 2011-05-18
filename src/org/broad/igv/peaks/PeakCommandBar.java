/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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
 * Created by JFormDesigner on Tue May 17 18:20:38 EDT 2011
 */

package org.broad.igv.peaks;

import java.awt.event.*;

import org.broad.igv.ui.IGV;

import java.awt.*;
import javax.swing.*;
import javax.swing.event.*;

/**
 * @author Stan Diamond
 */
public class PeakCommandBar extends JPanel {
    public PeakCommandBar() {
        initComponents();
        scoreSlider.setValue((int) PeakTrack.getScoreThreshold());
        foldChangeSlider.setValue((int) PeakTrack.getFoldChangeThreshold());
        colorByFoldChangeCB.setSelected(PeakTrack.getColorOption() == PeakTrack.ColorOption.FOLD_CHANGE);
        peaksCB.setSelected(PeakTrack.isShowPeaks());
        signalsCB.setSelected(PeakTrack.isShowSignals());
    }

    private void scoreSliderStateChanged(ChangeEvent e) {
        PeakTrack.setScoreThreshold(scoreSlider.getValue());
        IGV.getInstance().repaintDataPanels();
    }

    private void foldChangeSliderStateChanged(ChangeEvent e) {
        PeakTrack.setFoldChangeThreshold(foldChangeSlider.getValue());
        IGV.getInstance().repaintDataPanels();

    }

    private void peaksCBActionPerformed(ActionEvent e) {
        PeakTrack.setShowPeaks(peaksCB.isSelected());
        IGV.getInstance().repaint();
    }

    private void signalsCBActionPerformed(ActionEvent e) {
        PeakTrack.setShowSignals(signalsCB.isSelected());
        IGV.getInstance().repaint();
    }

    private void colorByFoldChangeCBActionPerformed(ActionEvent e) {
        if (colorByFoldChangeCB.isSelected()) {
            PeakTrack.setShadeOption(PeakTrack.ColorOption.FOLD_CHANGE);
        } else {
            PeakTrack.setShadeOption(PeakTrack.ColorOption.SCORE);
        }
        IGV.getInstance().repaint();
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        panel3 = new JPanel();
        label1 = new JLabel();
        foldChangeSlider = new JSlider();
        panel2 = new JPanel();
        label2 = new JLabel();
        scoreSlider = new JSlider();
        peaksCB = new JCheckBox();
        signalsCB = new JCheckBox();
        colorByFoldChangeCB = new JCheckBox();

        //======== this ========
        setLayout(new FlowLayout(FlowLayout.LEFT, 5, 0));

        //======== panel3 ========
        {
            panel3.setLayout(new FlowLayout());

            //---- label1 ----
            label1.setText("Fold change");
            panel3.add(label1);

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
            panel3.add(foldChangeSlider);
        }
        add(panel3);

        //======== panel2 ========
        {
            panel2.setLayout(new FlowLayout());

            //---- label2 ----
            label2.setText("Score");
            panel2.add(label2);

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
            panel2.add(scoreSlider);
        }
        add(panel2);

        //---- peaksCB ----
        peaksCB.setText("Show peaks");
        peaksCB.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                peaksCBActionPerformed(e);
            }
        });
        add(peaksCB);

        //---- signalsCB ----
        signalsCB.setText("Show signals");
        signalsCB.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                signalsCBActionPerformed(e);
            }
        });
        add(signalsCB);

        //---- colorByFoldChangeCB ----
        colorByFoldChangeCB.setText("Color by fold change");
        colorByFoldChangeCB.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                colorByFoldChangeCBActionPerformed(e);
            }
        });
        add(colorByFoldChangeCB);
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel panel3;
    private JLabel label1;
    private JSlider foldChangeSlider;
    private JPanel panel2;
    private JLabel label2;
    private JSlider scoreSlider;
    private JCheckBox peaksCB;
    private JCheckBox signalsCB;
    private JCheckBox colorByFoldChangeCB;
    // JFormDesigner - End of variables declaration  //GEN-END:variables
}
