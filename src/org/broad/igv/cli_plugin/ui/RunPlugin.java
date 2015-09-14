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
 * Created by JFormDesigner on Fri Aug 03 11:28:02 EDT 2012
 */

package org.broad.igv.cli_plugin.ui;

import org.broad.igv.PreferenceManager;
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
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.BrowserLauncher;
import org.broad.igv.util.FileUtils;
import org.broad.igv.variant.VariantTrack;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;
import java.util.*;
import java.util.List;

/**
 * @author jacob
 */
public class RunPlugin extends JDialog {

    private List<Argument> argumentList;
    private List<String> cmdList = new ArrayList<String>();
    private Map<Argument, ArgumentPanel> argumentComponents;

    private List<PluginSpecReader.Output> outputAttrs;
    private Map<PluginSpecReader.Output, ArgumentPanel> outputComponents;

    private String specPath;

    private String pluginId;
    private String toolName;
    private String commandName;
    private boolean forbidEmptyOutput;


    public RunPlugin(Frame owner, PluginSpecReader pluginSpecReader, final PluginSpecReader.Tool tool, PluginSpecReader.Command command) {
        super(owner);
        initComponents();

        specPath = pluginSpecReader.getSpecPath();
        argumentList = command.argumentList;
        outputAttrs = command.outputList;
        initArgumentComponents(pluginSpecReader, tool, command);

        if (tool.helpUrl != null) {
            helpButton.setEnabled(true);
            helpButton.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    try {
                        BrowserLauncher.openURL(tool.helpUrl);
                    } catch (IOException e1) {
                        MessageUtils.showErrorMessage(e1.getMessage(), e1);
                    }
                }
            });
        }
    }

    /**
     * Initialize the components based on the input arguments. Text inputs become text boxes,
     * track inputs become dropdown boxes, etc.
     *
     * @param pluginSpecReader
     * @param tool
     * @param command
     */
    private void initArgumentComponents(PluginSpecReader pluginSpecReader, PluginSpecReader.Tool tool, PluginSpecReader.Command command) {

        final String toolPath = pluginSpecReader.getToolPath(tool);

        this.pluginId = pluginSpecReader.getId();
        this.toolName = tool.name;
        this.forbidEmptyOutput = tool.forbidEmptyOutput;
        this.commandName = command.name;

        String[] cmdEls = new String[]{toolPath, command.cmd};
        if(tool.msgList != null && tool.msgList.size() > 0){
            this.cmdList.addAll(tool.msgList);
        }
        for (String cmdEl : cmdEls) {
            if (cmdEl != null && cmdEl.length() > 0) {
                this.cmdList.add(cmdEl);
            }
        }

        argumentComponents = new LinkedHashMap<Argument, ArgumentPanel>(this.argumentList.size());
        outputComponents = new LinkedHashMap<PluginSpecReader.Output, ArgumentPanel>(this.outputAttrs.size());

        String titleText = tool.name;
        if (commandName.length() > 0) {
            titleText += ": " + commandName;
        }
        setTitle(titleText);

        //Inputs
        Dimension minSize = getMinimumSize();
        for (Argument argument : argumentList) {

            if(argument.getType() == Argument.InputType.TEXT || argument.getType() == Argument.InputType.LONGTEXT){

                String defValue = argument.getDefaultValue();
                if(argument.isRemembered()){
                    String lastValue = PreferenceManager.getInstance().getArgumentValue(pluginId, toolName, commandName, argument.getId());
                    defValue = lastValue != null ? lastValue : defValue;
                }

                if(defValue != null && defValue.contains(Argument.TOOL_DIR_KEY)){
                    String toolDir = FileUtils.getParent(toolPath);
                    defValue = defValue.replace(Argument.TOOL_DIR_KEY, toolDir);
                }

                argument.setDefaultValue(defValue);
            }

            ArgumentPanel panel = ArgumentPanel.create(argument);
            if (panel != null) {
                argumentComponents.put(argument, panel);
                if (argument.isVisible()) {
                    this.contentPanel.add(panel);
                }
            }
        }

        //Outputs
        //This is somewhat hacky now, because we expect only one type of output
        //That being a track, where the user just needs to give a name
        for (PluginSpecReader.Output output : outputAttrs) {
            TextArgument panel = new TextArgument();
            panel.setArgName(output.name);
            String defValue = output.defaultValue != null ? output.defaultValue : tool.name + " " + commandName;
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

    /**
     * Retrieve the argument values from the UI components, and save whichever
     * are appropriate to preferences (for remembering)
     * @return
     */
    LinkedHashMap<Argument, Object> getArgumentValues(){
        LinkedHashMap<Argument, Object> argumentValues = new LinkedHashMap<Argument, Object>(argumentComponents.size());
        for (Map.Entry<Argument, ArgumentPanel> argComp : argumentComponents.entrySet()) {
            Object value = argComp.getValue().getValue();
            argumentValues.put(argComp.getKey(), value);

            //Save to preferences
            if(value instanceof String && argComp.getKey().isRemembered()){
                PreferenceManager.getInstance().putArgumentValue(pluginId, toolName, commandName, argComp.getKey().getId(), (String) value);
            }
        }
        return argumentValues;
    }

    private List<Track> genNewTracks() {
        //Retrieve the actual argument values
        LinkedHashMap<Argument, Object> argumentValues = getArgumentValues();

        List<Track> newTracks = new ArrayList<Track>(outputAttrs.size());

        for (PluginSpecReader.Output outputAttr : outputAttrs) {
            //TODO Hacky, only works for output components being TextArgument, which they are as of this comment writing
            String name = (String) outputComponents.get(outputAttr).getValue();

            Track newTrack = null;
            switch (outputAttr.type) {
                case FEATURE_TRACK:
                    PluginFeatureSource featSource1 = new PluginFeatureSource(cmdList, argumentValues, outputAttr, specPath, this.forbidEmptyOutput);
                    //featSource1.setQueryTracker(queryTracker);
                    FeatureSource featSource = new CachingFeatureSource(featSource1);
                    newTrack = new FeatureTrack(UUID.randomUUID().toString(), name, featSource);
                    break;
                case DATA_SOURCE_TRACK:
                    PluginDataSource dataSource = new PluginDataSource(GenomeManager.getInstance().getCurrentGenome(), cmdList, argumentValues, outputAttr, specPath);
                    //dataSource.setQueryTracker(queryTracker);
                    newTrack = new DataSourceTrack(null, UUID.randomUUID().toString(), name, dataSource);
                    break;
                case VARIANT_TRACK:
                    PluginFeatureSource VfeatSource1 = new PluginFeatureSource(cmdList, argumentValues, outputAttr, specPath, this.forbidEmptyOutput);
                    FeatureSource VfeatSource = new CachingFeatureSource(VfeatSource1);
                    newTrack = new VariantTrack(name, VfeatSource);
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
        List<Track> newTrack = genNewTracks();
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
        helpButton = new JButton();
        hSpacer1 = new JPanel(null);
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
                ((GridBagLayout) buttonBar.getLayout()).columnWidths = new int[]{0, 0, 85, 80};
                ((GridBagLayout) buttonBar.getLayout()).columnWeights = new double[]{0.0, 1.0, 0.0, 0.0};

                //---- helpButton ----
                helpButton.setText("Help");
                helpButton.setEnabled(false);
                buttonBar.add(helpButton, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0,
                        GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                        new Insets(0, 0, 0, 5), 0, 0));
                buttonBar.add(hSpacer1, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0,
                        GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                        new Insets(0, 0, 0, 5), 0, 0));

                //---- okButton ----
                okButton.setText("OK");
                okButton.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        okButtonActionPerformed(e);
                    }
                });
                buttonBar.add(okButton, new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0,
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
                buttonBar.add(cancelButton, new GridBagConstraints(3, 0, 1, 1, 0.0, 0.0,
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
    private JButton helpButton;
    private JPanel hSpacer1;
    private JButton okButton;
    private JButton cancelButton;
    // JFormDesigner - End of variables declaration  //GEN-END:variables
}
