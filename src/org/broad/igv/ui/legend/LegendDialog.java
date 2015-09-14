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
 * Created by JFormDesigner on Thu Jun 16 10:44:21 EDT 2011
 */

package org.broad.igv.ui.legend;

import java.awt.*;
import java.awt.Component;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.border.*;
import javax.swing.event.*;

import org.broad.igv.track.TrackType;
import org.broad.igv.ui.util.UIUtilities;
import org.jdesktop.layout.GroupLayout;
import org.jdesktop.layout.LayoutStyle;

/**
 * @author Stan Diamond
 */
public class LegendDialog extends JDialog {

    public LegendDialog(Frame owner) {
        super(owner);
        initComponents();
    }


    public LegendDialog(Frame owner, boolean b) {
        super(owner, b);
    }

    private void okButtonActionPerformed(ActionEvent e) {
        setVisible(false);
    }

    private void resetToDefaultActionPerformed(ActionEvent e) {
        String message = "This will reset all heatmap preferences to their default values.\n" +
                "Are you sure you want to continue?";

        boolean status = UIUtilities.showConfirmationDialog(this, message);
        if (status) {

            ((LegendPanel) copyNoCanvas).resetPreferencesToDefault();
            ((LegendPanel) expressionCanvas).resetPreferencesToDefault();
            ((LegendPanel) rnaiPanel).resetPreferencesToDefault();
            ((LegendPanel) lohCanvas).resetPreferencesToDefault();
            ((LegendPanel) methylationCanvas).resetPreferencesToDefault();
            ((LegendPanel) mutationCanvas).resetPreferencesToDefault();
        }

    }

    private void copyNumberButtonActionPerformed(ActionEvent e) {
        ((LegendPanel)copyNoCanvas).edit();
    }

    private void expressionButtonActionPerformed(ActionEvent e) {
        ((LegendPanel) expressionCanvas).edit();
    }

    private void rnaiButtonActionPerformed(ActionEvent e) {
        ((LegendPanel) rnaiPanel).edit();
    }

    private void methylationButtonActionPerformed(ActionEvent e) {
        ((LegendPanel) methylationCanvas).edit();
    }

    private void lohButtonActionPerformed(ActionEvent e) {
        ((LegendPanel) lohCanvas).edit();
    }

    private void mutationButtonActionPerformed(ActionEvent e) {
        ((LegendPanel) mutationCanvas).edit();
    }


    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        jLabel1 = new JLabel();
        panel1 = new JPanel();
        cnLabel = new JLabel();
        copyNoCanvas = new org.broad.igv.ui.legend.HeatmapLegendPanel(TrackType.COPY_NUMBER);
        legendLabel1 = new JLabel();
        methylationCanvas = new HeatmapLegendPanel(TrackType.DNA_METHYLATION);
        lohCanvas = new LohLegendPanel();
        legendLabel2 = new JLabel();
        jLabel2 = new JLabel();
        expressionCanvas = new HeatmapLegendPanel(TrackType.GENE_EXPRESSION);
        mutationCanvas = new MutationLegendPanel();
        legendLabel3 = new JLabel();
        legendLabel4 = new JLabel();
        rnaiPanel = new HeatmapLegendPanel(TrackType.RNAI);
        copyNumberButton = new JButton();
        expressionButton = new JButton();
        rnaiButton = new JButton();
        methylationButton = new JButton();
        lohButton = new JButton();
        mutationButton = new JButton();
        resetToDefault = new JButton();
        okButton = new JButton();
        jLabel3 = new JLabel();

        //======== this ========
        setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
        Container contentPane = getContentPane();
        contentPane.setLayout(null);

        //---- jLabel1 ----
        jLabel1.setFont(new Font("Lucida Sans", Font.BOLD, 18));
        jLabel1.setText("Legends");
        contentPane.add(jLabel1);
        jLabel1.setBounds(new Rectangle(new Point(325, 0), jLabel1.getPreferredSize()));

        //======== panel1 ========
        {
            panel1.setLayout(null);

            //---- cnLabel ----
            cnLabel.setText("Copy Number");
            panel1.add(cnLabel);
            cnLabel.setBounds(new Rectangle(new Point(0, 7), cnLabel.getPreferredSize()));

            //======== copyNoCanvas ========
            {
                copyNoCanvas.setBorder(null);

                GroupLayout copyNoCanvasLayout = new GroupLayout(copyNoCanvas);
                copyNoCanvas.setLayout(copyNoCanvasLayout);
                copyNoCanvasLayout.setHorizontalGroup(
                    copyNoCanvasLayout.createParallelGroup()
                        .add(0, 483, Short.MAX_VALUE)
                );
                copyNoCanvasLayout.setVerticalGroup(
                    copyNoCanvasLayout.createParallelGroup()
                        .add(0, 60, Short.MAX_VALUE)
                );
            }
            panel1.add(copyNoCanvas);
            copyNoCanvas.setBounds(110, 7, 485, 62);

            //---- legendLabel1 ----
            legendLabel1.setText("Expression");
            panel1.add(legendLabel1);
            legendLabel1.setBounds(new Rectangle(new Point(0, 87), legendLabel1.getPreferredSize()));

            //======== methylationCanvas ========
            {
                methylationCanvas.setMinimumSize(new Dimension(100, 100));
                methylationCanvas.setBorder(null);

                GroupLayout methylationCanvasLayout = new GroupLayout(methylationCanvas);
                methylationCanvas.setLayout(methylationCanvasLayout);
                methylationCanvasLayout.setHorizontalGroup(
                    methylationCanvasLayout.createParallelGroup()
                        .add(0, 483, Short.MAX_VALUE)
                );
                methylationCanvasLayout.setVerticalGroup(
                    methylationCanvasLayout.createParallelGroup()
                        .add(0, 60, Short.MAX_VALUE)
                );
            }
            panel1.add(methylationCanvas);
            methylationCanvas.setBounds(110, 247, 485, 62);

            //======== lohCanvas ========
            {
                lohCanvas.setBorder(null);

                GroupLayout lohCanvasLayout = new GroupLayout(lohCanvas);
                lohCanvas.setLayout(lohCanvasLayout);
                lohCanvasLayout.setHorizontalGroup(
                    lohCanvasLayout.createParallelGroup()
                        .add(0, 483, Short.MAX_VALUE)
                );
                lohCanvasLayout.setVerticalGroup(
                    lohCanvasLayout.createParallelGroup()
                        .add(0, 60, Short.MAX_VALUE)
                );
            }
            panel1.add(lohCanvas);
            lohCanvas.setBounds(110, 327, 485, 62);

            //---- legendLabel2 ----
            legendLabel2.setText("LOH");
            panel1.add(legendLabel2);
            legendLabel2.setBounds(0, 327, 45, legendLabel2.getPreferredSize().height);

            //---- jLabel2 ----
            jLabel2.setText("Methylation");
            panel1.add(jLabel2);
            jLabel2.setBounds(new Rectangle(new Point(0, 247), jLabel2.getPreferredSize()));

            //======== expressionCanvas ========
            {
                expressionCanvas.setBorder(null);

                GroupLayout expressionCanvasLayout = new GroupLayout(expressionCanvas);
                expressionCanvas.setLayout(expressionCanvasLayout);
                expressionCanvasLayout.setHorizontalGroup(
                    expressionCanvasLayout.createParallelGroup()
                        .add(0, 483, Short.MAX_VALUE)
                );
                expressionCanvasLayout.setVerticalGroup(
                    expressionCanvasLayout.createParallelGroup()
                        .add(0, 60, Short.MAX_VALUE)
                );
            }
            panel1.add(expressionCanvas);
            expressionCanvas.setBounds(110, 87, 485, 62);

            //======== mutationCanvas ========
            {
                mutationCanvas.setBorder(null);

                GroupLayout mutationCanvasLayout = new GroupLayout(mutationCanvas);
                mutationCanvas.setLayout(mutationCanvasLayout);
                mutationCanvasLayout.setHorizontalGroup(
                    mutationCanvasLayout.createParallelGroup()
                        .add(0, 485, Short.MAX_VALUE)
                );
                mutationCanvasLayout.setVerticalGroup(
                    mutationCanvasLayout.createParallelGroup()
                        .add(0, 62, Short.MAX_VALUE)
                );
            }
            panel1.add(mutationCanvas);
            mutationCanvas.setBounds(110, 407, 485, 62);

            //---- legendLabel3 ----
            legendLabel3.setText("Mutation");
            panel1.add(legendLabel3);
            legendLabel3.setBounds(new Rectangle(new Point(0, 430), legendLabel3.getPreferredSize()));

            //---- legendLabel4 ----
            legendLabel4.setText("RNAi");
            panel1.add(legendLabel4);
            legendLabel4.setBounds(new Rectangle(new Point(0, 167), legendLabel4.getPreferredSize()));

            //======== rnaiPanel ========
            {
                rnaiPanel.setBorder(null);

                GroupLayout rnaiPanelLayout = new GroupLayout(rnaiPanel);
                rnaiPanel.setLayout(rnaiPanelLayout);
                rnaiPanelLayout.setHorizontalGroup(
                    rnaiPanelLayout.createParallelGroup()
                        .add(0, 483, Short.MAX_VALUE)
                );
                rnaiPanelLayout.setVerticalGroup(
                    rnaiPanelLayout.createParallelGroup()
                        .add(0, 60, Short.MAX_VALUE)
                );
            }
            panel1.add(rnaiPanel);
            rnaiPanel.setBounds(110, 167, 485, 62);

            //---- copyNumberButton ----
            copyNumberButton.setText("Edit");
            copyNumberButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    copyNumberButtonActionPerformed(e);
                }
            });
            panel1.add(copyNumberButton);
            copyNumberButton.setBounds(new Rectangle(new Point(620, 7), copyNumberButton.getPreferredSize()));

            //---- expressionButton ----
            expressionButton.setText("Edit");
            expressionButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    expressionButtonActionPerformed(e);
                }
            });
            panel1.add(expressionButton);
            expressionButton.setBounds(new Rectangle(new Point(620, 87), expressionButton.getPreferredSize()));

            //---- rnaiButton ----
            rnaiButton.setText("Edit");
            rnaiButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    rnaiButtonActionPerformed(e);
                }
            });
            panel1.add(rnaiButton);
            rnaiButton.setBounds(new Rectangle(new Point(620, 167), rnaiButton.getPreferredSize()));

            //---- methylationButton ----
            methylationButton.setText("Edit");
            methylationButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    methylationButtonActionPerformed(e);
                }
            });
            panel1.add(methylationButton);
            methylationButton.setBounds(new Rectangle(new Point(620, 247), methylationButton.getPreferredSize()));

            //---- lohButton ----
            lohButton.setText("Edit");
            lohButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    lohButtonActionPerformed(e);
                }
            });
            panel1.add(lohButton);
            lohButton.setBounds(new Rectangle(new Point(620, 327), lohButton.getPreferredSize()));

            //---- mutationButton ----
            mutationButton.setText("Edit");
            mutationButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    mutationButtonActionPerformed(e);
                }
            });
            panel1.add(mutationButton);
            mutationButton.setBounds(new Rectangle(new Point(620, 424), mutationButton.getPreferredSize()));

            { // compute preferred size
                Dimension preferredSize = new Dimension();
                for(int i = 0; i < panel1.getComponentCount(); i++) {
                    Rectangle bounds = panel1.getComponent(i).getBounds();
                    preferredSize.width = Math.max(bounds.x + bounds.width, preferredSize.width);
                    preferredSize.height = Math.max(bounds.y + bounds.height, preferredSize.height);
                }
                Insets insets = panel1.getInsets();
                preferredSize.width += insets.right;
                preferredSize.height += insets.bottom;
                panel1.setMinimumSize(preferredSize);
                panel1.setPreferredSize(preferredSize);
            }
        }
        contentPane.add(panel1);
        panel1.setBounds(10, 30, 720, 495);

        //---- resetToDefault ----
        resetToDefault.setText("Reset to default");
        resetToDefault.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                resetToDefaultActionPerformed(e);
            }
        });
        contentPane.add(resetToDefault);
        resetToDefault.setBounds(new Rectangle(new Point(475, 575), resetToDefault.getPreferredSize()));

        //---- okButton ----
        okButton.setText("OK");
        okButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                okButtonActionPerformed(e);
            }
        });
        contentPane.add(okButton);
        okButton.setBounds(new Rectangle(new Point(635, 575), okButton.getPreferredSize()));

        //---- jLabel3 ----
        jLabel3.setText("<html>* click any item to bring up its editor");
        contentPane.add(jLabel3);
        jLabel3.setBounds(15, 550, 396, jLabel3.getPreferredSize().height);

        { // compute preferred size
            Dimension preferredSize = new Dimension();
            for(int i = 0; i < contentPane.getComponentCount(); i++) {
                Rectangle bounds = contentPane.getComponent(i).getBounds();
                preferredSize.width = Math.max(bounds.x + bounds.width, preferredSize.width);
                preferredSize.height = Math.max(bounds.y + bounds.height, preferredSize.height);
            }
            Insets insets = contentPane.getInsets();
            preferredSize.width += insets.right;
            preferredSize.height += insets.bottom;
            contentPane.setMinimumSize(preferredSize);
            contentPane.setPreferredSize(preferredSize);
        }
        setSize(755, 665);
        setLocationRelativeTo(getOwner());
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JLabel jLabel1;
    private JPanel panel1;
    private JLabel cnLabel;
    private JPanel copyNoCanvas;
    private JLabel legendLabel1;
    private JPanel methylationCanvas;
    private JPanel lohCanvas;
    private JLabel legendLabel2;
    private JLabel jLabel2;
    private JPanel expressionCanvas;
    private JPanel mutationCanvas;
    private JLabel legendLabel3;
    private JLabel legendLabel4;
    private JPanel rnaiPanel;
    private JButton copyNumberButton;
    private JButton expressionButton;
    private JButton rnaiButton;
    private JButton methylationButton;
    private JButton lohButton;
    private JButton mutationButton;
    private JButton resetToDefault;
    private JButton okButton;
    private JLabel jLabel3;
    // JFormDesigner - End of variables declaration  //GEN-END:variables

    public static void main(String args[]) {
        java.awt.EventQueue.invokeLater(new Runnable() {

            public void run() {
                LegendDialog  dialog = new LegendDialog(new javax.swing.JFrame());
                dialog.addWindowListener(new java.awt.event.WindowAdapter() {

                    public void windowClosing(java.awt.event.WindowEvent e) {
                        System.exit(0);
                    }
                });
                dialog.setVisible(true);
            }
        });
    }

}
