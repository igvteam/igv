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
 * Created by JFormDesigner on Mon Nov 26 16:20:15 EST 2012
 */

package org.broad.igv.plugin.mongovariant;

import com.mongodb.WriteResult;
import org.broad.igv.track.Track;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.variant.VariantTrack;
import htsjdk.tribble.util.ParsingUtils;
import org.broadinstitute.gatk.tools.walkers.na12878kb.core.MongoVariantContext;
import org.broadinstitute.gatk.tools.walkers.na12878kb.core.NA12878DBArgumentCollection;
import org.broadinstitute.gatk.tools.walkers.na12878kb.core.NA12878KnowledgeBase;
import org.broadinstitute.gatk.tools.walkers.na12878kb.core.TruthStatus;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.io.IOException;

/**
 * Dialog for reviewing a variant from a VCF file.
 * The users opinion is submitted to a MongoDB, specified by
 * system property {@link VariantReviewPlugin#DB_PATH_KEY}. Path can be relative to either the local machine,
 * or NA12878kb package or the GenomeAnalysisTK jar
 */
public class VariantReviewDialog extends JDialog {

    private VariantContext variantContext;
    private String userName;

    //Path to spec file describing database
    private String dbSpecPath;

    public VariantReviewDialog(Frame owner, VariantContext vc, String dbSpecPath) {
        super(owner);
        initComponents();

        truthField.setModel(new DefaultComboBoxModel(TruthStatus.values()));
        genotypeTypeField.setModel(new DefaultComboBoxModel(GenotypeType.values()));
        this.variantContext = vc;
        this.userName = System.getProperty("user.name", "unknown");
        this.dbSpecPath = dbSpecPath;

        initComponentData(vc);
    }

    private void initComponentData(VariantContext variant) {

        callsetField.setText(userName);

        chrField.setText(variant.getChr());
        startField.setText("" + variant.getStart());
        stopField.setText("" + variant.getEnd());

        truthField.setSelectedItem(TruthStatus.UNKNOWN);

        initGenotypeTypeField(variant);

        validate();
    }

    private void cancelButtonActionPerformed(ActionEvent e) {
        setVisible(false);
    }

    /**
     * Save information to Mongo
     *
     * @param e
     */
    private void okButtonActionPerformed(ActionEvent e) {

        String callsetName = callsetField.getText();
        TruthStatus truthStatus = (TruthStatus) truthField.getSelectedItem();
        boolean isComplexEvent = isComplexEventCB.isSelected();

        int allele0 = -1;
        int allele1 = -1;

        GenotypeType gtt = (GenotypeType) genotypeTypeField.getSelectedItem();
        switch (gtt) {
            case HOM_REF:
                allele0 = allele1 = 0;
                break;
            case HOM_VAR:
                allele0 = allele1 = 1;
                break;
            case HET:
                allele0 = 0;
                allele1 = 1;
                break;
        }

        MongoVariantContext mvc = VariantReviewSource.createMVC(allele0, allele1, callsetName, variantContext, truthStatus, isComplexEvent);
        String errorMessage = addCall(this.dbSpecPath, mvc);

        if (errorMessage != null) {
            MessageUtils.showErrorMessage(errorMessage, new IOException(errorMessage));
        } else {
            setVisible(false);
            //Find the track showing results, clear it to force a refresh
            for (Track t : IGV.getInstance().getFeatureTracks()) {
                if (t instanceof VariantTrack) {
                    ((VariantTrack) t).clearPackedFeatures();
                }
            }
            IGV.getInstance().repaintDataPanels();
        }
    }

    /**
     * Add the specified mvc to the specified database
     *
     * @param dbSpecPath
     * @param mvc
     * @return
     */
    static String addCall(String dbSpecPath, MongoVariantContext mvc) {

        NA12878DBArgumentCollection args = new NA12878DBArgumentCollection(dbSpecPath);

        String errorMessage = null;
        NA12878KnowledgeBase kb = null;
        try {
            kb = new NA12878KnowledgeBase(null, args);
            WriteResult wr = kb.addCall(mvc);
            errorMessage = wr.getError();
        } catch (Exception ex) {
            errorMessage = ex.getMessage();
            if (errorMessage == null) errorMessage = "" + ex;
        } finally {
            if (kb != null) kb.close();
        }

        return errorMessage;
    }

    private void truthFieldItemStateChanged(ItemEvent e) {
        if (truthField.getSelectedItem() == TruthStatus.FALSE_POSITIVE) {
            genotypeTypeField.setSelectedItem(GenotypeType.NO_CALL);
            genotypeTypeField.setEnabled(false);
        } else {
            genotypeTypeField.setEnabled(true);
            initGenotypeTypeField(this.variantContext);
        }

    }

    private void initGenotypeTypeField(VariantContext variant) {
        String mutationString = null;
        GenotypeType gtt = GenotypeType.NO_CALL;
        String prefSampleName = VariantReviewPlugin.getPreferentialSampleName();

        //If there is only one sample, or we find the preferential sample,
        //use that data.
        for (String sampleName : variant.getSampleNamesOrderedByName()) {
            boolean isPref = sampleName.equalsIgnoreCase(prefSampleName);
            if (isPref || mutationString == null) {
                mutationString = ParsingUtils.join("/", ParsingUtils.sortList(variant.getAlleles()));
                Genotype genotype = variant.getGenotype(sampleName);
                gtt = genotype.getType();
                if (isPref) break;
            } else {
                //If we have several samples with different mutations, don't know which
                //to pick. Make that obvious to the user
                if (gtt != variant.getGenotype(sampleName).getType()) {
                    mutationString = "./.";
                    gtt = GenotypeType.UNAVAILABLE;
                }
            }
        }

        genotypeTypeField.setSelectedItem(gtt);
        mutField.setText(mutationString);
        mutField.setToolTipText(mutationString);
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        dialogPane = new JPanel();
        contentPanel = new JPanel();
        panel1 = new JPanel();
        label4 = new JLabel();
        callsetField = new JTextField();
        hSpacer1 = new JPanel(null);
        panel2 = new JPanel();
        label5 = new JLabel();
        truthField = new JComboBox();
        panel3 = new JPanel();
        label6 = new JLabel();
        genotypeTypeField = new JComboBox();
        hSpacer2 = new JPanel(null);
        panel7 = new JPanel();
        label10 = new JLabel();
        mutField = new JLabel();
        hSpacer5 = new JPanel(null);
        panel4 = new JPanel();
        label7 = new JLabel();
        chrField = new JLabel();
        hSpacer3 = new JPanel(null);
        panel5 = new JPanel();
        label8 = new JLabel();
        startField = new JLabel();
        hSpacer4 = new JPanel(null);
        panel6 = new JPanel();
        label9 = new JLabel();
        stopField = new JLabel();
        panel8 = new JPanel();
        hSpacer6 = new JPanel(null);
        isComplexEventCB = new JCheckBox();
        hSpacer7 = new JPanel(null);
        buttonBar = new JPanel();
        okButton = new JButton();
        cancelButton = new JButton();

        //======== this ========
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== dialogPane ========
        {
            dialogPane.setBorder(new EmptyBorder(12, 12, 12, 12));
            dialogPane.setPreferredSize(new Dimension(600, 150));
            dialogPane.setLayout(new BorderLayout());

            //======== contentPanel ========
            {
                contentPanel.setLayout(new BoxLayout(contentPanel, BoxLayout.X_AXIS));

                //======== panel1 ========
                {
                    panel1.setLayout(new BoxLayout(panel1, BoxLayout.PAGE_AXIS));

                    //---- label4 ----
                    label4.setText("Callset");
                    label4.setHorizontalAlignment(SwingConstants.CENTER);
                    label4.setMaximumSize(new Dimension(80, 16));
                    label4.setAlignmentX(0.5F);
                    panel1.add(label4);

                    //---- callsetField ----
                    callsetField.setMaximumSize(new Dimension(200, 1000));
                    callsetField.setMinimumSize(new Dimension(100, 28));
                    callsetField.setPreferredSize(new Dimension(100, 28));
                    panel1.add(callsetField);
                }
                contentPanel.add(panel1);
                contentPanel.add(hSpacer1);

                //======== panel2 ========
                {
                    panel2.setLayout(new BoxLayout(panel2, BoxLayout.Y_AXIS));

                    //---- label5 ----
                    label5.setText("Truth");
                    label5.setHorizontalAlignment(SwingConstants.CENTER);
                    label5.setLabelFor(truthField);
                    label5.setHorizontalTextPosition(SwingConstants.CENTER);
                    label5.setAlignmentX(0.5F);
                    panel2.add(label5);

                    //---- truthField ----
                    truthField.addItemListener(new ItemListener() {
                        @Override
                        public void itemStateChanged(ItemEvent e) {
                            truthFieldItemStateChanged(e);
                        }
                    });
                    panel2.add(truthField);
                }
                contentPanel.add(panel2);

                //======== panel3 ========
                {
                    panel3.setLayout(new BoxLayout(panel3, BoxLayout.Y_AXIS));

                    //---- label6 ----
                    label6.setText("Genotype");
                    label6.setHorizontalAlignment(SwingConstants.CENTER);
                    label6.setAlignmentX(0.5F);
                    panel3.add(label6);
                    panel3.add(genotypeTypeField);
                }
                contentPanel.add(panel3);

                //---- hSpacer2 ----
                hSpacer2.setMinimumSize(new Dimension(20, 12));
                hSpacer2.setPreferredSize(new Dimension(20, 10));
                contentPanel.add(hSpacer2);

                //======== panel7 ========
                {
                    panel7.setMaximumSize(new Dimension(46, 1000));
                    panel7.setLayout(new BoxLayout(panel7, BoxLayout.Y_AXIS));

                    //---- label10 ----
                    label10.setText("Mut");
                    label10.setHorizontalAlignment(SwingConstants.LEFT);
                    label10.setMaximumSize(new Dimension(80, 16));
                    panel7.add(label10);

                    //---- mutField ----
                    mutField.setText("A/G");
                    mutField.setMaximumSize(new Dimension(46, 1000));
                    panel7.add(mutField);
                }
                contentPanel.add(panel7);

                //---- hSpacer5 ----
                hSpacer5.setMinimumSize(new Dimension(20, 12));
                hSpacer5.setPreferredSize(new Dimension(20, 10));
                contentPanel.add(hSpacer5);

                //======== panel4 ========
                {
                    panel4.setMaximumSize(new Dimension(46, 1000));
                    panel4.setLayout(new BoxLayout(panel4, BoxLayout.Y_AXIS));

                    //---- label7 ----
                    label7.setText("Chr");
                    label7.setHorizontalAlignment(SwingConstants.LEFT);
                    label7.setVerticalAlignment(SwingConstants.TOP);
                    panel4.add(label7);

                    //---- chrField ----
                    chrField.setText("testChr");
                    chrField.setAlignmentY(0.0F);
                    chrField.setMaximumSize(new Dimension(46, 1000));
                    panel4.add(chrField);
                }
                contentPanel.add(panel4);

                //---- hSpacer3 ----
                hSpacer3.setMinimumSize(new Dimension(20, 12));
                hSpacer3.setPreferredSize(new Dimension(20, 10));
                contentPanel.add(hSpacer3);

                //======== panel5 ========
                {
                    panel5.setMaximumSize(new Dimension(500, 1000));
                    panel5.setMinimumSize(new Dimension(80, 32));
                    panel5.setPreferredSize(new Dimension(80, 32));
                    panel5.setLayout(new BoxLayout(panel5, BoxLayout.Y_AXIS));

                    //---- label8 ----
                    label8.setText("Start");
                    label8.setHorizontalAlignment(SwingConstants.LEFT);
                    label8.setMaximumSize(new Dimension(100, 16));
                    panel5.add(label8);

                    //---- startField ----
                    startField.setText("54321");
                    startField.setMaximumSize(new Dimension(500, 1000));
                    startField.setMinimumSize(new Dimension(80, 16));
                    startField.setPreferredSize(new Dimension(80, 16));
                    startField.setHorizontalTextPosition(SwingConstants.LEFT);
                    startField.setHorizontalAlignment(SwingConstants.LEFT);
                    panel5.add(startField);
                }
                contentPanel.add(panel5);

                //---- hSpacer4 ----
                hSpacer4.setMinimumSize(new Dimension(20, 12));
                hSpacer4.setPreferredSize(new Dimension(20, 10));
                contentPanel.add(hSpacer4);

                //======== panel6 ========
                {
                    panel6.setMaximumSize(new Dimension(500, 1000));
                    panel6.setMinimumSize(new Dimension(100, 32));
                    panel6.setPreferredSize(new Dimension(100, 32));
                    panel6.setLayout(new BoxLayout(panel6, BoxLayout.Y_AXIS));

                    //---- label9 ----
                    label9.setText("Stop");
                    label9.setHorizontalAlignment(SwingConstants.LEFT);
                    label9.setMaximumSize(new Dimension(100, 16));
                    panel6.add(label9);

                    //---- stopField ----
                    stopField.setText("12345");
                    stopField.setMaximumSize(new Dimension(500, 1000));
                    stopField.setHorizontalTextPosition(SwingConstants.LEADING);
                    stopField.setHorizontalAlignment(SwingConstants.LEFT);
                    stopField.setMinimumSize(new Dimension(100, 16));
                    stopField.setPreferredSize(new Dimension(100, 16));
                    panel6.add(stopField);
                }
                contentPanel.add(panel6);
            }
            dialogPane.add(contentPanel, BorderLayout.NORTH);

            //======== panel8 ========
            {
                panel8.setLayout(new BoxLayout(panel8, BoxLayout.X_AXIS));

                //---- hSpacer6 ----
                hSpacer6.setMaximumSize(new Dimension(300, 32767));
                hSpacer6.setPreferredSize(new Dimension(115, 10));
                panel8.add(hSpacer6);

                //---- isComplexEventCB ----
                isComplexEventCB.setText("Complex Event");
                panel8.add(isComplexEventCB);
                panel8.add(hSpacer7);
            }
            dialogPane.add(panel8, BorderLayout.CENTER);

            //======== buttonBar ========
            {
                buttonBar.setBorder(new EmptyBorder(12, 0, 0, 0));
                buttonBar.setLayout(new GridBagLayout());
                ((GridBagLayout) buttonBar.getLayout()).columnWidths = new int[]{0, 85, 80};
                ((GridBagLayout) buttonBar.getLayout()).columnWeights = new double[]{1.0, 0.0, 0.0};

                //---- okButton ----
                okButton.setText("Save");
                okButton.addActionListener(new ActionListener() {
                    @Override
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
                    @Override
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
        setSize(700, 160);
        setLocationRelativeTo(getOwner());
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel dialogPane;
    private JPanel contentPanel;
    private JPanel panel1;
    private JLabel label4;
    private JTextField callsetField;
    private JPanel hSpacer1;
    private JPanel panel2;
    private JLabel label5;
    private JComboBox truthField;
    private JPanel panel3;
    private JLabel label6;
    private JComboBox genotypeTypeField;
    private JPanel hSpacer2;
    private JPanel panel7;
    private JLabel label10;
    private JLabel mutField;
    private JPanel hSpacer5;
    private JPanel panel4;
    private JLabel label7;
    private JLabel chrField;
    private JPanel hSpacer3;
    private JPanel panel5;
    private JLabel label8;
    private JLabel startField;
    private JPanel hSpacer4;
    private JPanel panel6;
    private JLabel label9;
    private JLabel stopField;
    private JPanel panel8;
    private JPanel hSpacer6;
    private JCheckBox isComplexEventCB;
    private JPanel hSpacer7;
    private JPanel buttonBar;
    private JButton okButton;
    private JButton cancelButton;
    // JFormDesigner - End of variables declaration  //GEN-END:variables
}
