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
 * Created by JFormDesigner on Tue Jan 29 15:22:45 EST 2013
 */

package org.broad.igv.tools.motiffinder;

import htsjdk.samtools.util.SequenceUtil;
import org.broad.igv.annotations.ForTesting;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.ui.util.UIUtilities;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.StringUtils;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import javax.swing.event.CaretEvent;
import javax.swing.event.CaretListener;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;
import java.util.regex.PatternSyntaxException;

/**
 * Dialog so the user can enter a pattern that can be used
 * to search a nucleotide sequence.
 *
 * @author Jacob Silterra
 */
public class MotifFinderDialog extends JDialog {

    private static Map<String, String> letterToRegex;
    private static Set<String> validIUPACInputStrings;

    static {
        initLetterToRegex();
    }


    private static void initLetterToRegex() {
        letterToRegex = ParsingUtils.loadIUPACMap();
        validIUPACInputStrings = new HashSet<String>(letterToRegex.size());
        for (String key : letterToRegex.keySet()) {
            validIUPACInputStrings.add(key.toUpperCase());
        }
    }

    /**
     * Returns true if character c is a valid IUPAC Ambiguity code
     *
     * @param c Single character
     * @return
     */
    public static boolean isIUPACChar(String c) {
        return validIUPACInputStrings.contains(c);
    }

    private String[] posTrackNames;
    private String[] negTrackNames;

    private String[] inputPatterns;

    public MotifFinderDialog(Frame owner) {
        super(owner);
        initComponents();
    }

    public String[] getInputPattern() {
        return inputPatterns;
    }

    public String[] getPosTrackName() {
        return posTrackNames;
    }

    public String[] getNegTrackName() {
        return negTrackNames;
    }

    /**
     * Replace the ambiguity codes in the motif
     * with regular expression equivalents
     *
     * @param motif
     * @return
     */
    static String convertMotifToRegex(String motif) {
        String output = motif;
        int outloc = 0;
        for (int inloc = 0; inloc < motif.length(); inloc++) {

            String inchar = motif.substring(inloc, inloc + 1);
            String rep = letterToRegex.get(inchar);

            output = output.substring(0, outloc) + rep + motif.substring(inloc + 1);
            outloc += rep.length();
        }
        return output;
    }

    private void cancelButtonActionPerformed(ActionEvent e) {
        this.setVisible(false);
    }

    void okButtonActionPerformed(ActionEvent e) {

        this.inputPatterns = null;

        //Split into lines, each one will be a new pair of tracks
        String[] lines = patternField.getText().split("[\\r\\n]+");
        String[] patterns = new String[lines.length];
        boolean isMultiMatch = lines.length >= 2;

        this.posTrackNames = new String[lines.length];
        this.negTrackNames = new String[lines.length];

        for(int ii=0; ii < lines.length; ii++){
            String strPattern = lines[ii].toUpperCase();

            boolean isIUPAC = checkIUPACPatternValid(strPattern);
            boolean isRegex = checkNucleotideRegex(strPattern);
            boolean patternIsValid = isIUPAC || isRegex;

            if (!patternIsValid) {
                MessageUtils.showMessage("Please enter a valid pattern.\n" +
                        "Patterns using IUPAC ambiguity codes should contain no special characters.\n" +
                        "Regular expressions should contain only 'ACTGN' in addition to special characters.\n" +
                        strPattern + " is invalid");
                return;
            }

            if (isIUPAC) {
                strPattern = convertMotifToRegex(strPattern);
            }

            if(isMultiMatch){
                String posName = getPosNameFromPattern(strPattern);
                this.posTrackNames[ii] = posName;
                this.negTrackNames[ii] = getNegNameFromPositive(posName);
            }else{
                this.posTrackNames[ii] = posNameField.getText();
                this.negTrackNames[ii] = negNameField.getText();
                if (this.posTrackNames[ii].equalsIgnoreCase(negTrackNames[ii])) {
                    MessageUtils.showMessage("Track names must be different");
                    return;
                }
            }
            patterns[ii] = strPattern;

        }
        this.inputPatterns = patterns;

        this.setVisible(false);
    }

    /**
     * Determines whether this string pattern is interpretable as
     * a set of IUPAC nucleotide characters
     *
     * @param strPattern Upper case string pattern
     * @return
     */
    static boolean checkIUPACPatternValid(String strPattern) {
        for (int ii = 0; ii < strPattern.length(); ii++) {
            String c = strPattern.substring(ii, ii + 1);

            if (!isIUPACChar(c)) {
                return false;
            }
        }
        return true;
    }

    /**
     * Determine whether it's a valid regex.
     * Also, any letters should be one of AGCTN (case insensitive)
     *
     * @param strPattern
     * @return
     */
    static boolean checkNucleotideRegex(String strPattern) {
        try {
            //First check if it's valid regex
            Pattern pattern = Pattern.compile(strPattern);
            byte[] bytes = strPattern.getBytes();
            for (byte c : bytes) {
                if (Character.isLetter(c)) {
                    boolean validBase = SequenceUtil.isValidBase(c);
                    validBase |= c == 'N';
                    if (!validBase) return false;
                }
            }
        } catch (PatternSyntaxException e) {
            return false;
        }

        return true;
    }

    /**
     * Whether we are matching multiple patterns (ie whether there are newlines in the text field)
     * @return
     */
    private boolean isMultiMatch(){
        String patternText = MotifFinderDialog.this.patternField.getText();
        return patternText.contains("\n") || patternText.contains("\r");
    }

    private void updateNegNameFieldFromPattern() {
        UIUtilities.invokeOnEventThread(new Runnable() {
            @Override
            public void run() {
                String posText = MotifFinderDialog.this.posNameField.getText();
                MotifFinderDialog.this.negNameField.setText(getNegNameFromPositive(posText));
                MotifFinderDialog.this.negNameField.setEnabled(!isMultiMatch());
            }
        });
    }

    private void updatePosNameFieldFromPattern() {
        UIUtilities.invokeOnEventThread(new Runnable() {
            @Override
            public void run() {
                String posNameText = "Auto";
                boolean hasNewlines = isMultiMatch();
                if(!hasNewlines){
                    String patternText = MotifFinderDialog.this.patternField.getText();
                    posNameText = StringUtils.checkLength(patternText, MaxTrackNameLength);
                }

                MotifFinderDialog.this.posNameField.setEnabled(!isMultiMatch());
                MotifFinderDialog.this.posNameField.setText(posNameText);
            }
        });
    }

    private String getPosNameFromPattern(String patternText){
        return StringUtils.checkLength(patternText, MaxTrackNameLength);
    }

    private String getNegNameFromPositive(String posText){
        return posText + " Negative";
    }

    static final int MaxTrackNameLength = 100;

    private void patternFieldCaretUpdate(CaretEvent e) {
        updatePosNameFieldFromPattern();
        updateNegNameFieldFromPattern();
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        dialogPane = new JPanel();
        contentPanel = new JPanel();
        label2 = new JLabel();
        label4 = new JLabel();
        patternField = new JTextArea();
        textArea1 = new JTextArea();
        vSpacer1 = new JPanel(null);
        panel1 = new JPanel();
        label1 = new JLabel();
        posNameField = new JTextField();
        panel2 = new JPanel();
        label3 = new JLabel();
        negNameField = new JTextField();
        buttonBar = new JPanel();
        okButton = new JButton();
        cancelButton = new JButton();

        //======== this ========
        setModal(true);
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== dialogPane ========
        {
            dialogPane.setBorder(new EmptyBorder(12, 12, 12, 12));
            dialogPane.setLayout(new BorderLayout());

            //======== contentPanel ========
            {
                contentPanel.setAlignmentX(0.0F);
                contentPanel.setLayout(new BoxLayout(contentPanel, BoxLayout.Y_AXIS));

                //---- label2 ----
                label2.setText("Search Pattern:");
                label2.setLabelFor(patternField);
                label2.setHorizontalTextPosition(SwingConstants.LEFT);
                label2.setHorizontalAlignment(SwingConstants.LEFT);
                label2.setAlignmentX(1.0F);
                label2.setMaximumSize(new Dimension(374, 16));
                label2.setPreferredSize(new Dimension(374, 16));
                contentPanel.add(label2);
                contentPanel.add(label4);

                //---- patternField ----
                patternField.setToolTipText("Enter multiple patterns, separated by newlines");
                patternField.setRows(2);
                patternField.addCaretListener(new CaretListener() {
                    @Override
                    public void caretUpdate(CaretEvent e) {
                        patternFieldCaretUpdate(e);
                    }
                });
                contentPanel.add(patternField);

                //---- textArea1 ----
                textArea1.setText("Enter nucleotide sequence (e.g. ACCGCT),  or nucleotide sequence with IUPAC ambiguity codes (e.g. AAARNR),  or regular expression of nucleotides (e.g. TATAAA(A){3,}). ");
                textArea1.setEditable(false);
                textArea1.setBackground(new Color(238, 238, 238));
                textArea1.setLineWrap(true);
                textArea1.setWrapStyleWord(true);
                textArea1.setFont(textArea1.getFont().deriveFont(textArea1.getFont().getStyle() | Font.ITALIC, textArea1.getFont().getSize() - 2f));
                textArea1.setMargin(new Insets(0, 25, 0, 0));
                textArea1.setFocusable(false);
                contentPanel.add(textArea1);

                //---- vSpacer1 ----
                vSpacer1.setMinimumSize(new Dimension(12, 20));
                vSpacer1.setPreferredSize(new Dimension(10, 20));
                contentPanel.add(vSpacer1);

                //======== panel1 ========
                {
                    panel1.setLayout(new BoxLayout(panel1, BoxLayout.X_AXIS));

                    //---- label1 ----
                    label1.setText("Positive Strand Track Name:");
                    label1.setLabelFor(posNameField);
                    label1.setHorizontalTextPosition(SwingConstants.LEFT);
                    label1.setHorizontalAlignment(SwingConstants.LEFT);
                    label1.setMaximumSize(new Dimension(374, 16));
                    label1.setPreferredSize(new Dimension(200, 16));
                    label1.setAlignmentX(1.0F);
                    panel1.add(label1);
                    panel1.add(posNameField);
                }
                contentPanel.add(panel1);

                //======== panel2 ========
                {
                    panel2.setLayout(new BoxLayout(panel2, BoxLayout.X_AXIS));

                    //---- label3 ----
                    label3.setText("Negative Strand Track Name:");
                    label3.setLabelFor(negNameField);
                    label3.setHorizontalTextPosition(SwingConstants.LEFT);
                    label3.setHorizontalAlignment(SwingConstants.LEFT);
                    label3.setMaximumSize(new Dimension(374, 16));
                    label3.setPreferredSize(new Dimension(200, 16));
                    label3.setAlignmentX(1.0F);
                    panel2.add(label3);
                    panel2.add(negNameField);
                }
                contentPanel.add(panel2);
            }
            dialogPane.add(contentPanel, BorderLayout.NORTH);

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
        setSize(450, 300);
        setLocationRelativeTo(getOwner());
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel dialogPane;
    private JPanel contentPanel;
    private JLabel label2;
    private JLabel label4;
    private JTextArea patternField;
    private JTextArea textArea1;
    private JPanel vSpacer1;
    private JPanel panel1;
    private JLabel label1;
    private JTextField posNameField;
    private JPanel panel2;
    private JLabel label3;
    private JTextField negNameField;
    private JPanel buttonBar;
    private JButton okButton;
    private JButton cancelButton;
    // JFormDesigner - End of variables declaration  //GEN-END:variables

    @ForTesting
    void setPatternFieldText(String text) {
        this.patternField.setText(text);
    }
}
