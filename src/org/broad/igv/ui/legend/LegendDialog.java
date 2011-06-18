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


    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        jLabel1 = new JLabel();
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
        jLabel1.setBounds(new Rectangle(new Point(392, 0), jLabel1.getPreferredSize()));

        //---- cnLabel ----
        cnLabel.setText("Copy Number");
        contentPane.add(cnLabel);
        cnLabel.setBounds(new Rectangle(new Point(10, 40), cnLabel.getPreferredSize()));

        //======== copyNoCanvas ========
        {
            copyNoCanvas.setBorder(null);

            GroupLayout copyNoCanvasLayout = new GroupLayout(copyNoCanvas);
            copyNoCanvas.setLayout(copyNoCanvasLayout);
            copyNoCanvasLayout.setHorizontalGroup(
                copyNoCanvasLayout.createParallelGroup()
                    .add(0, 588, Short.MAX_VALUE)
            );
            copyNoCanvasLayout.setVerticalGroup(
                copyNoCanvasLayout.createParallelGroup()
                    .add(0, 60, Short.MAX_VALUE)
            );
        }
        contentPane.add(copyNoCanvas);
        copyNoCanvas.setBounds(120, 38, 590, 62);

        //---- legendLabel1 ----
        legendLabel1.setText("Expression");
        contentPane.add(legendLabel1);
        legendLabel1.setBounds(new Rectangle(new Point(10, 120), legendLabel1.getPreferredSize()));

        //======== methylationCanvas ========
        {
            methylationCanvas.setMinimumSize(new Dimension(100, 100));
            methylationCanvas.setBorder(null);

            GroupLayout methylationCanvasLayout = new GroupLayout(methylationCanvas);
            methylationCanvas.setLayout(methylationCanvasLayout);
            methylationCanvasLayout.setHorizontalGroup(
                methylationCanvasLayout.createParallelGroup()
                    .add(0, 588, Short.MAX_VALUE)
            );
            methylationCanvasLayout.setVerticalGroup(
                methylationCanvasLayout.createParallelGroup()
                    .add(0, 60, Short.MAX_VALUE)
            );
        }
        contentPane.add(methylationCanvas);
        methylationCanvas.setBounds(120, 272, 590, 62);

        //======== lohCanvas ========
        {
            lohCanvas.setBorder(null);

            GroupLayout lohCanvasLayout = new GroupLayout(lohCanvas);
            lohCanvas.setLayout(lohCanvasLayout);
            lohCanvasLayout.setHorizontalGroup(
                lohCanvasLayout.createParallelGroup()
                    .add(0, 588, Short.MAX_VALUE)
            );
            lohCanvasLayout.setVerticalGroup(
                lohCanvasLayout.createParallelGroup()
                    .add(0, 60, Short.MAX_VALUE)
            );
        }
        contentPane.add(lohCanvas);
        lohCanvas.setBounds(120, 350, 590, 62);

        //---- legendLabel2 ----
        legendLabel2.setText("LOH");
        contentPane.add(legendLabel2);
        legendLabel2.setBounds(10, 375, 45, legendLabel2.getPreferredSize().height);

        //---- jLabel2 ----
        jLabel2.setText("Methylation");
        contentPane.add(jLabel2);
        jLabel2.setBounds(new Rectangle(new Point(10, 285), jLabel2.getPreferredSize()));

        //======== expressionCanvas ========
        {
            expressionCanvas.setBorder(null);

            GroupLayout expressionCanvasLayout = new GroupLayout(expressionCanvas);
            expressionCanvas.setLayout(expressionCanvasLayout);
            expressionCanvasLayout.setHorizontalGroup(
                expressionCanvasLayout.createParallelGroup()
                    .add(0, 588, Short.MAX_VALUE)
            );
            expressionCanvasLayout.setVerticalGroup(
                expressionCanvasLayout.createParallelGroup()
                    .add(0, 60, Short.MAX_VALUE)
            );
        }
        contentPane.add(expressionCanvas);
        expressionCanvas.setBounds(120, 116, 590, 62);

        //======== mutationCanvas ========
        {
            mutationCanvas.setBorder(null);

            GroupLayout mutationCanvasLayout = new GroupLayout(mutationCanvas);
            mutationCanvas.setLayout(mutationCanvasLayout);
            mutationCanvasLayout.setHorizontalGroup(
                mutationCanvasLayout.createParallelGroup()
                    .add(0, 588, Short.MAX_VALUE)
            );
            mutationCanvasLayout.setVerticalGroup(
                mutationCanvasLayout.createParallelGroup()
                    .add(0, 60, Short.MAX_VALUE)
            );
        }
        contentPane.add(mutationCanvas);
        mutationCanvas.setBounds(120, 428, 590, 62);

        //---- legendLabel3 ----
        legendLabel3.setText("Mutation");
        contentPane.add(legendLabel3);
        legendLabel3.setBounds(new Rectangle(new Point(10, 455), legendLabel3.getPreferredSize()));

        //---- legendLabel4 ----
        legendLabel4.setText("RNAi");
        contentPane.add(legendLabel4);
        legendLabel4.setBounds(new Rectangle(new Point(10, 200), legendLabel4.getPreferredSize()));

        //======== rnaiPanel ========
        {
            rnaiPanel.setBorder(null);

            GroupLayout rnaiPanelLayout = new GroupLayout(rnaiPanel);
            rnaiPanel.setLayout(rnaiPanelLayout);
            rnaiPanelLayout.setHorizontalGroup(
                rnaiPanelLayout.createParallelGroup()
                    .add(0, 588, Short.MAX_VALUE)
            );
            rnaiPanelLayout.setVerticalGroup(
                rnaiPanelLayout.createParallelGroup()
                    .add(0, 60, Short.MAX_VALUE)
            );
        }
        contentPane.add(rnaiPanel);
        rnaiPanel.setBounds(120, 194, 590, 62);

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
