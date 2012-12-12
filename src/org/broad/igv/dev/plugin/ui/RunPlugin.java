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
 * Created by JFormDesigner on Fri Aug 03 11:28:02 EDT 2012
 */

package org.broad.igv.dev.plugin.ui;

import org.broad.igv.dev.plugin.Argument;
import org.broad.igv.dev.plugin.PluginFeatureSource;
import org.broad.igv.dev.plugin.PluginSpecReader;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.Track;
import org.broad.igv.ui.IGV;
import org.w3c.dom.Element;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * @author User #2
 */
public class RunPlugin extends JDialog {

    private List<Argument> argumentList;
    private List<String> cmd = new ArrayList<String>();
    private Map<Argument, ArgumentPanel> argumentComponents;
    private Map<String, String> parsingAttrs;

    private String specPath;


    public RunPlugin(Frame owner, PluginSpecReader pluginSpecReader, Element tool, Element command) {
        super(owner);
        initComponents();

        final String toolPath = pluginSpecReader.getToolPath(tool);
        final String toolName = tool.getAttribute(PluginSpecReader.TOOL_NAME_KEY);
        final String cmdName = command.getAttribute("name");
        final String cmdVal = command.getAttribute("cmd");

        specPath = pluginSpecReader.getSpecPath();
        argumentList = pluginSpecReader.getArguments(tool, command);
        parsingAttrs = pluginSpecReader.getParsingAttributes(tool, command);
        initArgumentComponents(toolName, toolPath, cmdName, cmdVal);
    }

    private void initArgumentComponents(String toolName, String toolPath, String cmdName, String cmdVal) {


        if (toolPath.length() > 0) {
            this.cmd.add(toolPath);
        }
        if (cmdVal != null && cmdVal.trim().length() > 0) {
            this.cmd.add(cmdVal);
        }
        argumentComponents = new LinkedHashMap<Argument, ArgumentPanel>(this.argumentList.size());

        String titleText = toolName;
        if (cmdName.length() > 0) {
            titleText += ": " + cmdName;
        }
        setTitle(titleText);

        Dimension minSize = getMinimumSize();
        for (Argument argument : argumentList) {
            ArgumentPanel panel = ArgumentPanel.create(argument);
            if (panel != null) {
                this.contentPanel.add(panel);
                argumentComponents.put(argument, panel);
            }
        }

        this.validate();

        double minWidth = minSize.getWidth();
        for (ArgumentPanel panel : argumentComponents.values()) {

            //If track names are long, feature track combo box can get big
            //TODO Multi-intersect box is in a scrollpanel so it doesn't enlarge quite the same. Maybe we want that, maybe not
            minWidth = Math.max(minWidth, panel.getMinimumSize().getWidth());

            minSize.setSize(minWidth, minSize.getHeight() + panel.getHeight());
            this.setMinimumSize(minSize);
        }
        this.validate();

        outputName.setText(toolName + " " + cmdName);
    }

    private Track genNewTrack() {
        //Retrieve the actual argument values
        LinkedHashMap<Argument, Object> argumentValues = new LinkedHashMap<Argument, Object>(argumentComponents.size());
        for (Map.Entry<Argument, ArgumentPanel> argComp : argumentComponents.entrySet()) {
            Object value = argComp.getValue().getValue();
            argumentValues.put(argComp.getKey(), value);
        }

        String name = outputName.getText();

        //TODO PluginDataSource is already written, just need to know when to use it
        PluginFeatureSource source = new PluginFeatureSource(cmd, argumentValues, parsingAttrs, specPath);
        FeatureTrack newTrack = new FeatureTrack(name, name, source);
        return newTrack;
    }

    private void cancelButtonActionPerformed(ActionEvent e) {
        setVisible(false);
    }

    void okButtonActionPerformed(ActionEvent e) {
        Track newTrack = genNewTrack();
        IGV.getInstance().getTrackPanel(IGV.FEATURE_PANEL_NAME).addTrack(newTrack);

        this.setVisible(false);

        IGV.getInstance().repaint();
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        dialogPane = new JPanel();
        contentPanel = new JPanel();
        panel1 = new JPanel();
        label1 = new JLabel();
        outputName = new JTextField();
        buttonBar = new JPanel();
        okButton = new JButton();
        cancelButton = new JButton();

        //======== this ========
        setModal(true);
        setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== dialogPane ========
        {
            dialogPane.setBorder(new EmptyBorder(12, 12, 12, 12));
            dialogPane.setPreferredSize(new Dimension(300, 100));
            dialogPane.setLayout(new BorderLayout());

            //======== contentPanel ========
            {
                contentPanel.setMaximumSize(new Dimension(2000000, 16));
                contentPanel.setLayout(new BoxLayout(contentPanel, BoxLayout.Y_AXIS));
            }
            dialogPane.add(contentPanel, BorderLayout.NORTH);

            //======== panel1 ========
            {
                panel1.setMaximumSize(new Dimension(2000, 28));
                panel1.setMinimumSize(new Dimension(200, 28));
                panel1.setPreferredSize(new Dimension(200, 28));
                panel1.setLayout(new BoxLayout(panel1, BoxLayout.X_AXIS));

                //---- label1 ----
                label1.setText("Output Track Name:");
                panel1.add(label1);

                //---- outputName ----
                outputName.setText("result");
                outputName.setMaximumSize(new Dimension(500, 28));
                panel1.add(outputName);
            }
            dialogPane.add(panel1, BorderLayout.CENTER);

            //======== buttonBar ========
            {
                buttonBar.setBorder(new EmptyBorder(12, 0, 0, 0));
                buttonBar.setLayout(new GridBagLayout());
                ((GridBagLayout) buttonBar.getLayout()).columnWidths = new int[]{0, 85, 80};
                ((GridBagLayout) buttonBar.getLayout()).columnWeights = new double[]{1.0, 0.0, 0.0};

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
        pack();
        setLocationRelativeTo(getOwner());
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel dialogPane;
    private JPanel contentPanel;
    private JPanel panel1;
    private JLabel label1;
    private JTextField outputName;
    private JPanel buttonBar;
    private JButton okButton;
    private JButton cancelButton;
    // JFormDesigner - End of variables declaration  //GEN-END:variables
}
