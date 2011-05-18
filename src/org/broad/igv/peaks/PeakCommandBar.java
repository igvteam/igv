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
import javax.swing.border.*;

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
        this.colorByScoresButton.setSelected(PeakTrack.getColorOption() == PeakTrack.ColorOption.SCORE);
        this.colorByChangeButton.setSelected(PeakTrack.getColorOption() == PeakTrack.ColorOption.FOLD_CHANGE);
        this.peaksButton.setSelected(PeakTrack.isShowPeaks() && !PeakTrack.isShowSignals());
        this.signalsButton.setSelected(!PeakTrack.isShowPeaks() && PeakTrack.isShowSignals());
        this.bothButton.setSelected(PeakTrack.isShowPeaks() && PeakTrack.isShowSignals());

        //peaksCB.setSelected(PeakTrack.isShowPeaks());
        //signalsCB.setSelected(PeakTrack.isShowSignals());
    }

    private void scoreSliderStateChanged(ChangeEvent e) {
        PeakTrack.setScoreThreshold(scoreSlider.getValue());
        IGV.getInstance().repaintDataPanels();
    }

    private void foldChangeSliderStateChanged(ChangeEvent e) {
        PeakTrack.setFoldChangeThreshold(foldChangeSlider.getValue());
        IGV.getInstance().repaintDataPanels();

    }

    private void radioButtonActionPerformed(ActionEvent e) {
        if (bothButton.isSelected()) {
            PeakTrack.setShowPeaks(true);
            PeakTrack.setShowSignals(true);
        } else if (peaksButton.isSelected()) {
            PeakTrack.setShowPeaks(true);
            PeakTrack.setShowSignals(false);

        } else {
            PeakTrack.setShowPeaks(false);
            PeakTrack.setShowSignals(true);
        }
        IGV.getInstance().repaint();
    }

    /* private void peaksCBActionPerformed(ActionEvent e) {
        PeakTrack.setShowPeaks(peaksCB.isSelected());
        IGV.getInstance().repaint();
    }

    private void signalsCBActionPerformed(ActionEvent e) {
        PeakTrack.setShowSignals(signalsCB.isSelected());
        IGV.getInstance().repaint();
    }
    */

    private void colorByFoldChangeCBActionPerformed(ActionEvent e) {
    }

    private void colorByActionPeformed(ActionEvent e) {
        if (colorByChangeButton.isSelected()) {
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
        hSpacer1 = new JPanel(null);
        panel1 = new JPanel();
        label3 = new JLabel();
        peaksButton = new JRadioButton();
        signalsButton = new JRadioButton();
        bothButton = new JRadioButton();
        hSpacer2 = new JPanel(null);
        panel4 = new JPanel();
        label4 = new JLabel();
        colorByScoresButton = new JRadioButton();
        colorByChangeButton = new JRadioButton();

        //======== this ========
        setLayout(new FlowLayout(FlowLayout.LEFT, 5, 0));

        //======== panel3 ========
        {
            panel3.setLayout(new FlowLayout());

            //---- label1 ----
            label1.setText("<html>Fold change<br>threshold");
            panel3.add(label1);

            //---- foldChangeSlider ----
            foldChangeSlider.setPaintTicks(true);
            foldChangeSlider.setToolTipText("Adjust score threshold");
            foldChangeSlider.setMajorTickSpacing(2);
            foldChangeSlider.setMinorTickSpacing(1);
            foldChangeSlider.setMaximum(10);
            foldChangeSlider.setValue(0);
            foldChangeSlider.setPaintLabels(true);
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
            label2.setText("<html>Score<br>threshold");
            panel2.add(label2);

            //---- scoreSlider ----
            scoreSlider.setPaintTicks(true);
            scoreSlider.setToolTipText("Adjust score threshold");
            scoreSlider.setMajorTickSpacing(20);
            scoreSlider.setPaintLabels(true);
            scoreSlider.addChangeListener(new ChangeListener() {
                public void stateChanged(ChangeEvent e) {
                    scoreSliderStateChanged(e);
                }
            });
            panel2.add(scoreSlider);
            panel2.add(hSpacer1);
        }
        add(panel2);

        //======== panel1 ========
        {
            panel1.setBorder(new EtchedBorder(EtchedBorder.RAISED));
            panel1.setLayout(new FlowLayout());

            //---- label3 ----
            label3.setText("Show:");
            panel1.add(label3);

            //---- peaksButton ----
            peaksButton.setText("Peaks");
            peaksButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    radioButtonActionPerformed(e);
                }
            });
            panel1.add(peaksButton);

            //---- signalsButton ----
            signalsButton.setText("Signals");
            signalsButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    radioButtonActionPerformed(e);
                }
            });
            panel1.add(signalsButton);

            //---- bothButton ----
            bothButton.setText("Both");
            bothButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    radioButtonActionPerformed(e);
                }
            });
            panel1.add(bothButton);
        }
        add(panel1);
        add(hSpacer2);

        //======== panel4 ========
        {
            panel4.setBorder(new EtchedBorder(EtchedBorder.RAISED));
            panel4.setLayout(new FlowLayout(FlowLayout.LEFT));

            //---- label4 ----
            label4.setText("Color by:");
            panel4.add(label4);

            //---- colorByScoresButton ----
            colorByScoresButton.setText("Scores");
            colorByScoresButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    colorByActionPeformed(e);
                }
            });
            panel4.add(colorByScoresButton);

            //---- colorByChangeButton ----
            colorByChangeButton.setText("Change");
            colorByChangeButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    colorByActionPeformed(e);
                }
            });
            panel4.add(colorByChangeButton);
        }
        add(panel4);

        //---- buttonGroup1 ----
        ButtonGroup buttonGroup1 = new ButtonGroup();
        buttonGroup1.add(peaksButton);
        buttonGroup1.add(signalsButton);
        buttonGroup1.add(bothButton);

        //---- buttonGroup2 ----
        ButtonGroup buttonGroup2 = new ButtonGroup();
        buttonGroup2.add(colorByScoresButton);
        buttonGroup2.add(colorByChangeButton);
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
    private JPanel hSpacer1;
    private JPanel panel1;
    private JLabel label3;
    private JRadioButton peaksButton;
    private JRadioButton signalsButton;
    private JRadioButton bothButton;
    private JPanel hSpacer2;
    private JPanel panel4;
    private JLabel label4;
    private JRadioButton colorByScoresButton;
    private JRadioButton colorByChangeButton;
    // JFormDesigner - End of variables declaration  //GEN-END:variables
}
