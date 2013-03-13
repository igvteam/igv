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

package org.broad.igv.tools.sequencematch;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * Dialog so the user can enter a pattern that can be used
 * to search a nucleotide sequence.
 *
 * @author Jacob Silterra
 */
public class SequenceMatchDialog extends JDialog {

    private static Map<String, String> letterToRegex;
    private static Set<String> validInputStrings;

    static{
        initLetterToRegex();
    }

    private static final String codeFilePath = "resources/iupac_regex_table.txt";
    private static void initLetterToRegex() {
        letterToRegex = loadMap(SequenceMatchSource.class.getResourceAsStream(codeFilePath));
        validInputStrings = new HashSet<String>(letterToRegex.size());
        for(String key: letterToRegex.keySet()){
            validInputStrings.add(key.toUpperCase());
        }
    }

    public static boolean isValidString(String c) {
        return validInputStrings.contains(c);
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

    private String trackName;
    private String inputPattern;

    public SequenceMatchDialog(Frame owner) {
        super(owner);
        initComponents();
    }

    public String getInputPattern() {
        return inputPattern;
    }

    public String getTrackName() {
        return trackName;
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

        this.trackName = null;
        this.inputPattern = null;

        String strPattern = patternField.getText().toUpperCase();

        //Don't validate regex here, too hard
        boolean patternIsValid = true;
        if(ambiguityCodeButton.isSelected()){
            patternIsValid = checkPatternValid(strPattern);
            if(patternIsValid){
                strPattern = convertMotifToRegex(strPattern);
            }
        }

        if(!patternIsValid){
            // TODO warn user
        }else{
            this.trackName = nameField.getText();
            this.inputPattern = strPattern;
            this.setVisible(false);
        }
    }

    private boolean checkPatternValid(String strPattern) {
        for(int ii=0; ii < strPattern.length(); ii++){
            String c = strPattern.substring(ii, ii + 1);
            if(!isValidString(c)){
                return false;
            }
        }
        return true;
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        dialogPane = new JPanel();
        contentPanel = new JPanel();
        textArea1 = new JTextArea();
        label1 = new JLabel();
        nameField = new JTextField();
        label2 = new JLabel();
        label4 = new JLabel();
        patternField = new JTextField();
        radioButtonPanel = new JPanel();
        ambiguityCodeButton = new JRadioButton();
        regexButton = new JRadioButton();
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

                //---- textArea1 ----
                textArea1.setText("Enter the pattern for which to search. Matches in the current genome sequence will be displayed in a new track. ");
                textArea1.setEditable(false);
                textArea1.setBackground(new Color(238, 238, 238));
                textArea1.setLineWrap(true);
                contentPanel.add(textArea1);

                //---- label1 ----
                label1.setText("Track Name:");
                label1.setLabelFor(nameField);
                label1.setHorizontalTextPosition(SwingConstants.LEFT);
                label1.setHorizontalAlignment(SwingConstants.LEFT);
                label1.setMaximumSize(new Dimension(374, 16));
                label1.setPreferredSize(new Dimension(374, 16));
                label1.setAlignmentX(1.0F);
                contentPanel.add(label1);
                contentPanel.add(nameField);

                //---- label2 ----
                label2.setText("Pattern:");
                label2.setLabelFor(patternField);
                label2.setHorizontalTextPosition(SwingConstants.LEFT);
                label2.setHorizontalAlignment(SwingConstants.LEFT);
                label2.setAlignmentX(1.0F);
                label2.setMaximumSize(new Dimension(374, 16));
                label2.setPreferredSize(new Dimension(374, 16));
                contentPanel.add(label2);
                contentPanel.add(label4);
                contentPanel.add(patternField);

                //======== radioButtonPanel ========
                {
                    radioButtonPanel.setAlignmentX(1.0F);
                    radioButtonPanel.setMaximumSize(new Dimension(200, 46));
                    radioButtonPanel.setPreferredSize(new Dimension(200, 46));
                    radioButtonPanel.setLayout(new BoxLayout(radioButtonPanel, BoxLayout.Y_AXIS));

                    //---- ambiguityCodeButton ----
                    ambiguityCodeButton.setText("IUPAC Ambiguity Codes");
                    ambiguityCodeButton.setSelected(true);
                    ambiguityCodeButton.setToolTipText("Use IUPAC ambiguity codes");
                    radioButtonPanel.add(ambiguityCodeButton);

                    //---- regexButton ----
                    regexButton.setText("Regular Expressions");
                    regexButton.setToolTipText("Use java regular expression syntax");
                    radioButtonPanel.add(regexButton);
                }
                contentPanel.add(radioButtonPanel);
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

        //---- buttonGroup1 ----
        ButtonGroup buttonGroup1 = new ButtonGroup();
        buttonGroup1.add(ambiguityCodeButton);
        buttonGroup1.add(regexButton);
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel dialogPane;
    private JPanel contentPanel;
    private JTextArea textArea1;
    private JLabel label1;
    private JTextField nameField;
    private JLabel label2;
    private JLabel label4;
    private JTextField patternField;
    private JPanel radioButtonPanel;
    private JRadioButton ambiguityCodeButton;
    private JRadioButton regexButton;
    private JPanel buttonBar;
    private JButton okButton;
    private JButton cancelButton;
    // JFormDesigner - End of variables declaration  //GEN-END:variables
}
