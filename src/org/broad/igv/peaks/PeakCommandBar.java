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
 * Created by JFormDesigner on Tue May 17 18:20:38 EDT 2011
 */

package org.broad.igv.peaks;

import java.awt.event.*;
import javax.swing.border.*;

import org.broad.igv.feature.RegionOfInterest;
import org.broad.igv.track.RegionScoreType;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.util.MessageUtils;

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

    private void colorByActionPeformed(ActionEvent e) {
        if (colorByChangeButton.isSelected()) {
            PeakTrack.setShadeOption(PeakTrack.ColorOption.FOLD_CHANGE);
        } else {
            PeakTrack.setShadeOption(PeakTrack.ColorOption.SCORE);
        }
        IGV.getInstance().repaint();

    }

    private void sortButtonActionPerformed(ActionEvent e) {
        // Create a region based on current view

        // TODO -- do right thing if in geneListView
        if (FrameManager.isGeneListMode()) {
            MessageUtils.showMessage("To sort in gene list mode right-click the header of a gene panel.");
        } else {
            ReferenceFrame frame = FrameManager.getDefaultFrame();
            RegionOfInterest roi = new RegionOfInterest(frame.getChrName(), (int) frame.getOrigin(), (int) frame.getEnd(), "");
            IGV.getInstance().sortByRegionScore(roi, RegionScoreType.SCORE, frame);
        }
    }

    private void signalsButtonActionPerformed(ActionEvent e) {
        colorByScoresButton.setEnabled(false);
        colorByChangeButton.setEnabled(false);
    }

    private void peaksButtonActionPerformed(ActionEvent e) {
        colorByScoresButton.setEnabled(true);
        colorByChangeButton.setEnabled(true);
    }

    private void bothButtonActionPerformed(ActionEvent e) {
        colorByScoresButton.setEnabled(true);
        colorByChangeButton.setEnabled(true);
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
        colorByPanel = new JPanel();
        label4 = new JLabel();
        colorByScoresButton = new JRadioButton();
        colorByChangeButton = new JRadioButton();
        sortButton = new JButton();

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
            foldChangeSlider.setPreferredSize(new Dimension(150, 52));
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
            scoreSlider.setPreferredSize(new Dimension(150, 52));
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
            panel1.setMinimumSize(new Dimension(250, 37));
            panel1.setLayout(new FlowLayout());

            //---- label3 ----
            label3.setText("Show:");
            panel1.add(label3);

            //---- peaksButton ----
            peaksButton.setText("Peaks");
            peaksButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    radioButtonActionPerformed(e);
                    peaksButtonActionPerformed(e);
                }
            });
            panel1.add(peaksButton);

            //---- signalsButton ----
            signalsButton.setText("Signals");
            signalsButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    radioButtonActionPerformed(e);
                    signalsButtonActionPerformed(e);
                }
            });
            panel1.add(signalsButton);

            //---- bothButton ----
            bothButton.setText("Both");
            bothButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    radioButtonActionPerformed(e);
                    bothButtonActionPerformed(e);
                }
            });
            panel1.add(bothButton);
        }
        add(panel1);
        add(hSpacer2);

        //======== colorByPanel ========
        {
            colorByPanel.setBorder(new EtchedBorder(EtchedBorder.RAISED));
            colorByPanel.setLayout(new FlowLayout(FlowLayout.LEFT));

            //---- label4 ----
            label4.setText("Color peaks by:");
            colorByPanel.add(label4);

            //---- colorByScoresButton ----
            colorByScoresButton.setText("Factor");
            colorByScoresButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    colorByActionPeformed(e);
                }
            });
            colorByPanel.add(colorByScoresButton);

            //---- colorByChangeButton ----
            colorByChangeButton.setText("Change");
            colorByChangeButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    colorByActionPeformed(e);
                }
            });
            colorByPanel.add(colorByChangeButton);
        }
        add(colorByPanel);

        //---- sortButton ----
        sortButton.setText("Sort");
        sortButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                sortButtonActionPerformed(e);
            }
        });
        add(sortButton);

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
    private JPanel colorByPanel;
    private JLabel label4;
    private JRadioButton colorByScoresButton;
    private JRadioButton colorByChangeButton;
    private JButton sortButton;
    // JFormDesigner - End of variables declaration  //GEN-END:variables
}
