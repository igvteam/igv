/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

/*
 * Created by JFormDesigner on Tue Jan 29 15:22:45 EST 2013
 */

package org.broad.igv.tools.motiffinder;

import net.sf.samtools.util.SequenceUtil;
import org.broad.igv.ui.util.MessageUtils;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.*;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
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

    static{
        initLetterToRegex();
    }

    private static final String codeFilePath = "resources/iupac_regex_table.txt";
    private static void initLetterToRegex() {
        letterToRegex = loadMap(MotifFinderSource.class.getResourceAsStream(codeFilePath));
        validIUPACInputStrings = new HashSet<String>(letterToRegex.size());
        for(String key: letterToRegex.keySet()){
            validIUPACInputStrings.add(key.toUpperCase());
        }
    }

    /**
     * Returns true if character c is a valid IUPAC Ambiguity code
     * @param c Single character
     * @return
     */
    public static boolean isIUPACChar(String c) {
        return validIUPACInputStrings.contains(c);
    }

    /**
     * TODO Move this to someplace more general, use it wherever we store lots of this kind of data
     * @param inputStream
     * @return
     */
    public static Map<String, String> loadMap(InputStream inputStream){
        BufferedReader reader = null;
        Map<String, String> map = new HashMap<String, String>();
        try {
            reader = new BufferedReader(new InputStreamReader(inputStream));
            String nextLine = null;
            while ((nextLine = reader.readLine()) != null) {
                if(nextLine.startsWith("#")) continue;

                String[] tokens = nextLine.split("=");
                if (tokens.length == 2) {
                    map.put(tokens[0], tokens[1]);
                }else{
                    throw new IllegalArgumentException("Incorrect number of tokens at line: " + nextLine);
                }
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        } finally {
            try {
                if (reader != null) {
                    reader.close();
                }
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        }

        return map;
    }

    private String posTrackName;
    private String negTrackName;

    private String inputPattern;

    public MotifFinderDialog(Frame owner) {
        super(owner);
        initComponents();
    }

    public String getInputPattern() {
        return inputPattern;
    }

    public String getPosTrackName() {
        return posTrackName;
    }

    public String getNegTrackName() {
        return negTrackName;
    }

    /**
     * Replace the ambiguity codes in the motif
     * with regular expression equivalents
     * @param motif
     * @return
     */
    static String convertMotifToRegex(String motif){
        String output = motif;
        int outloc = 0;
        for(int inloc=0; inloc < motif.length(); inloc++){

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

    private void okButtonActionPerformed(ActionEvent e) {

        this.posTrackName = null;
        this.inputPattern = null;

        String strPattern = patternField.getText().toUpperCase();

        boolean isIUPAC = checkIUPACPatternValid(strPattern);
        boolean isRegex = checkNucleotideRegex(strPattern);
        boolean patternIsValid = isIUPAC || isRegex;

        if(!patternIsValid){
            MessageUtils.showMessage("Please enter a valid pattern.\n" +
                    "Patterns using IUPAC ambiguity codes should contain no special characters.\n" +
                    "Regular expressions should contain only 'ACTGN' in addition to special characters");
            return;
        }

        if (isIUPAC) {
            strPattern = convertMotifToRegex(strPattern);
        }

        this.posTrackName = posNameField.getText();
        this.negTrackName = negNameField.getText();
        if(this.posTrackName.equalsIgnoreCase(negTrackName)){
            MessageUtils.showMessage("Track names must be different");
            return;
        }
        this.inputPattern = strPattern;
        this.setVisible(false);

    }

    /**
     * Determines whether this string pattern is interpretable as
     * a set of IUPAC nucleotide characters
     * @param strPattern Upper case string pattern
     * @return
     */
    static boolean checkIUPACPatternValid(String strPattern) {
        for(int ii=0; ii < strPattern.length(); ii++){
            String c = strPattern.substring(ii, ii + 1);
            if(!isIUPACChar(c)){
                return false;
            }
        }
        return true;
    }

    /**
     * Determine whether it's a valid regex.
     * Also, any letters should be one of AGCTN (case insensitive)
     * @param strPattern
     * @return
     */
    static boolean checkNucleotideRegex(String strPattern){
        try {
            //First check if it's valid regex
            Pattern pattern = Pattern.compile(strPattern);
            byte[] bytes = strPattern.getBytes();
            for(byte c: bytes){
                if(Character.isLetter(c)){
                    boolean validBase = SequenceUtil.isValidBase(c);
                    validBase |= c == 'N';
                    if(!validBase) return false;
                }
            }
        } catch (PatternSyntaxException e) {
            return false;
        }

        return true;
    }

    private void posNameFieldFocusLost(FocusEvent e) {
        updateNegNameFieldFromPos();
    }

    private void updateNegNameFieldFromPos(){
        String posText = this.posNameField.getText();
        this.negNameField.setText(posText + " Negative");
    }

    private void posNameFieldInputMethodTextChanged(InputMethodEvent e) {
        updateNegNameFieldFromPos();
    }

    private void posNameFieldKeyReleased(KeyEvent e) {
        if(e.getKeyCode() == KeyEvent.VK_ENTER){
            updateNegNameFieldFromPos();
        }
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        dialogPane = new JPanel();
        contentPanel = new JPanel();
        label2 = new JLabel();
        label4 = new JLabel();
        patternField = new JTextField();
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
                contentPanel.add(patternField);

                //---- textArea1 ----
                textArea1.setText("Enter nucleotide sequence (e.g. ACCGCT), or nucleotide sequence with IUPAC ambiguity codes (e.g. AAARNR), or regular expression of nucleotides (e.g. TATAA(A){3,}+) ");
                textArea1.setEditable(false);
                textArea1.setBackground(new Color(238, 238, 238));
                textArea1.setLineWrap(true);
                textArea1.setWrapStyleWord(true);
                textArea1.setFont(textArea1.getFont().deriveFont(textArea1.getFont().getStyle() | Font.ITALIC, textArea1.getFont().getSize() - 2f));
                textArea1.setMargin(new Insets(0, 25, 0, 5));
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
                    label1.setText("Positive Strand Result Track:");
                    label1.setLabelFor(posNameField);
                    label1.setHorizontalTextPosition(SwingConstants.LEFT);
                    label1.setHorizontalAlignment(SwingConstants.LEFT);
                    label1.setMaximumSize(new Dimension(374, 16));
                    label1.setPreferredSize(new Dimension(200, 16));
                    label1.setAlignmentX(1.0F);
                    panel1.add(label1);

                    //---- posNameField ----
                    posNameField.addFocusListener(new FocusAdapter() {
                        @Override
                        public void focusLost(FocusEvent e) {
                            posNameFieldFocusLost(e);
                        }
                    });
                    posNameField.addInputMethodListener(new InputMethodListener() {
                        @Override
                        public void inputMethodTextChanged(InputMethodEvent e) {
                            posNameFieldInputMethodTextChanged(e);
                        }

                        @Override
                        public void caretPositionChanged(InputMethodEvent e) {
                        }
                    });
                    posNameField.addKeyListener(new KeyAdapter() {
                        @Override
                        public void keyReleased(KeyEvent e) {
                            posNameFieldKeyReleased(e);
                        }
                    });
                    panel1.add(posNameField);
                }
                contentPanel.add(panel1);

                //======== panel2 ========
                {
                    panel2.setLayout(new BoxLayout(panel2, BoxLayout.X_AXIS));

                    //---- label3 ----
                    label3.setText("Negative Strand Result Track:");
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
        setSize(400, 300);
        setLocationRelativeTo(getOwner());
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel dialogPane;
    private JPanel contentPanel;
    private JLabel label2;
    private JLabel label4;
    private JTextField patternField;
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
}
