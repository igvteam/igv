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

package org.broad.igv.cli_plugin.ui;

import org.broad.igv.cli_plugin.Argument;
import org.broad.igv.cli_plugin.PluginDataSource;
import org.broad.igv.cli_plugin.PluginFeatureSource;
import org.broad.igv.cli_plugin.PluginSpecReader;
import org.broad.igv.feature.CachingFeatureSource;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.track.DataSourceTrack;
import org.broad.igv.track.FeatureSource;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.Track;
import org.broad.igv.ui.IGV;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.*;
import java.util.List;

/**
 * @author jacob
 */
public class RunPlugin extends JDialog {

    private List<Argument> argumentList;
    private List<String> cmd = new ArrayList<String>();
    private Map<Argument, ArgumentPanel> argumentComponents;

    private List<PluginSpecReader.Output> outputAttrs;
    private Map<PluginSpecReader.Output, ArgumentPanel> outputComponents;

    private String specPath;


    public RunPlugin(Frame owner, PluginSpecReader pluginSpecReader, PluginSpecReader.Tool tool, PluginSpecReader.Command command) {
        super(owner);
        initComponents();

        final String toolPath = pluginSpecReader.getToolPath(tool);
        final String toolName = tool.name;
        final String cmdName = command.name;
        final String cmdVal = command.cmd;

        specPath = pluginSpecReader.getSpecPath();
        argumentList = command.argumentList;
        outputAttrs = command.outputList;
        initArgumentComponents(toolName, toolPath, cmdName, cmdVal);
    }

    /**
     * Initialize the components based on the input arguments. Text inputs become text boxes,
     * track inputs become dropdown boxes, etc.
     * @param toolName
     * @param toolPath
     * @param cmdName
     * @param cmdVal
     */
    private void initArgumentComponents(String toolName, String toolPath, String cmdName, String cmdVal) {

        if (toolPath != null && toolPath.length() > 0) {
            this.cmd.add(toolPath);
        }
        if (cmdVal != null && cmdVal.trim().length() > 0) {
            this.cmd.add(cmdVal);
        }
        argumentComponents = new LinkedHashMap<Argument, ArgumentPanel>(this.argumentList.size());
        outputComponents = new LinkedHashMap<PluginSpecReader.Output, ArgumentPanel>(this.outputAttrs.size());

        String titleText = toolName;
        if (cmdName.length() > 0) {
            titleText += ": " + cmdName;
        }
        setTitle(titleText);

        //Inputs
        Dimension minSize = getMinimumSize();
        for (Argument argument : argumentList) {
            ArgumentPanel panel = ArgumentPanel.create(argument);
            if (panel != null) {
                argumentComponents.put(argument, panel);
                if(argument.isVisible()){
                    this.contentPanel.add(panel);
                }
            }
        }

        //Outputs
        //This is somewhat hacky now, because we expect only one type of output
        //That being a track, where the user just needs to give a name
        for(PluginSpecReader.Output output: outputAttrs){
            TextArgument panel = new TextArgument();
            panel.setArgName(output.name);
            String defValue = output.defaultValue != null ? output.defaultValue : toolName + " " + cmdName;
            panel.setValue(defValue);
            this.contentPanel.add(panel);
            outputComponents.put(output, panel);
        }

        this.validate();

        List<ArgumentPanel> components = new ArrayList<ArgumentPanel>(argumentComponents.values());
        components.addAll(outputComponents.values());
        double minWidth = minSize.getWidth();
        for (ArgumentPanel panel : components) {

            //If track names are long, feature track combo box can get big
            //TODO Multi-intersect box is in a scrollpanel so it doesn't enlarge quite the same. Maybe we want that, maybe not
            minWidth = Math.max(minWidth, panel.getMinimumSize().getWidth());

            minSize.setSize(minWidth, minSize.getHeight() + panel.getHeight());
            this.setMinimumSize(minSize);
        }
        this.validate();
    }

    private List<Track> genNewTrack() {
        //Retrieve the actual argument values
        LinkedHashMap<Argument, Object> argumentValues = new LinkedHashMap<Argument, Object>(argumentComponents.size());
        for (Map.Entry<Argument, ArgumentPanel> argComp : argumentComponents.entrySet()) {
            Object value = argComp.getValue().getValue();
            argumentValues.put(argComp.getKey(), value);
        }

        List<Track> newTracks = new ArrayList<Track>(outputAttrs.size());

        //QueryTracker queryTracker = QueryTracker.get();

        for(PluginSpecReader.Output outputAttr: outputAttrs){
            //TODO Hacky, only works for output components being TextArgument, which they are as of this comment writing
            String name = (String) outputComponents.get(outputAttr).getValue();

            Track newTrack = null;
            switch (outputAttr.type){
                case FEATURE_TRACK:
                    PluginFeatureSource featSource1 = new PluginFeatureSource(cmd, argumentValues, outputAttr, specPath);
                    //featSource1.setQueryTracker(queryTracker);
                    FeatureSource featSource = new CachingFeatureSource(featSource1);
                    newTrack = new FeatureTrack(UUID.randomUUID().toString(), name, featSource);
                    break;
                case DATA_SOURCE_TRACK:
                    PluginDataSource dataSource = new PluginDataSource(GenomeManager.getInstance().getCurrentGenome(), cmd, argumentValues, outputAttr, specPath);
                    //dataSource.setQueryTracker(queryTracker);
                    newTrack = new DataSourceTrack(null, UUID.randomUUID().toString(), name, dataSource);
                    break;
            }
            newTracks.add(newTrack);
        }
        return newTracks;
    }

    private void cancelButtonActionPerformed(ActionEvent e) {
        setVisible(false);
    }

    void okButtonActionPerformed(ActionEvent e) {
        List<Track> newTrack = genNewTrack();
        IGV.getInstance().getTrackPanel(IGV.FEATURE_PANEL_NAME).addTracks(newTrack);

        this.setVisible(false);

        IGV.getInstance().repaint();
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        dialogPane = new JPanel();
        contentPanel = new JPanel();
        vSpacer1 = new JPanel(null);
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
            dialogPane.add(vSpacer1, BorderLayout.CENTER);

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
        pack();
        setLocationRelativeTo(getOwner());
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel dialogPane;
    private JPanel contentPanel;
    private JPanel vSpacer1;
    private JPanel buttonBar;
    private JButton okButton;
    private JButton cancelButton;
    // JFormDesigner - End of variables declaration  //GEN-END:variables
}
