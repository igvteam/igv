/*
 * Created by JFormDesigner on Wed Apr 11 16:44:25 EDT 2012
 */

package org.broad.igv.feature.genome;

import java.awt.event.*;

import org.broad.igv.PreferenceManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.FileDialogUtils;
import org.broad.igv.util.FileUtils;

import java.awt.*;
import java.io.File;
import javax.swing.*;
import javax.swing.border.*;

/**
 * @author Stan Diamond
 */
public class GenomeBuilderDialog extends JDialog {

    boolean canceled;
    File archiveFile;

    public GenomeBuilderDialog(Frame owner, IGV igv) {
        super(owner);
        setModal(true);
        initComponents();
        this.panel1.setIgv(igv);
    }

    public boolean isCanceled() {
        return canceled;
    }

    public File getArchiveFile() {
        return archiveFile;
    }

    private void okButtonActionPerformed(ActionEvent e) {

        archiveFile = chooseArchiveFile();
        canceled = (archiveFile == null);
        setVisible(false);
        dispose();
    }

    private void cancelButtonActionPerformed(ActionEvent e) {
        canceled = true;
        setVisible(false);
        dispose();
    }

    private File chooseArchiveFile() {
        File dir = PreferenceManager.getInstance().getLastGenomeImportDirectory();
        File initFile = new File(getGenomeId() + ".genome");
        return FileDialogUtils.chooseFile("Save .genome file", dir, initFile, FileDialog.SAVE);
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        dialogPane = new JPanel();
        contentPanel = new JPanel();
        panel1 = new GenomeBuilderPane();
        buttonBar = new JPanel();
        okButton = new JButton();
        cancelButton = new JButton();

        //======== this ========
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== dialogPane ========
        {
            dialogPane.setBorder(new EmptyBorder(12, 12, 12, 12));
            dialogPane.setLayout(new BorderLayout());

            //======== contentPanel ========
            {
                contentPanel.setLayout(new BorderLayout());
                contentPanel.add(panel1, BorderLayout.CENTER);
            }
            dialogPane.add(contentPanel, BorderLayout.CENTER);

            //======== buttonBar ========
            {
                buttonBar.setBorder(new EmptyBorder(12, 0, 0, 0));
                buttonBar.setLayout(new GridBagLayout());
                ((GridBagLayout) buttonBar.getLayout()).columnWidths = new int[]{0, 85, 80};
                ((GridBagLayout) buttonBar.getLayout()).columnWeights = new double[]{1.0, 0.0, 0.0};

                //---- okButton ----
                okButton.setText("OK");
                okButton.addActionListener(new ActionListener() {
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
        setSize(870, 430);
        setLocationRelativeTo(getOwner());
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel dialogPane;
    private JPanel contentPanel;
    private GenomeBuilderPane panel1;
    private JPanel buttonBar;
    private JButton okButton;
    private JButton cancelButton;
    // JFormDesigner - End of variables declaration  //GEN-END:variables

    public String getCytobandFileName() {
        return panel1.getCytobandFileName();
    }

    public String getRefFlatFileName() {
        return panel1.getRefFlatFileName();
    }

    public String getFastaFileName() {
        return panel1.getFastaFileName();
    }

    public String getChrAliasFileName() {
        return panel1.getChrAliasFileName();
    }


    public String getGenomeDisplayName() {
        return panel1.getGenomeDisplayName();
    }

    public String getGenomeId() {
        return panel1.getGenomeId();
    }

    public String getArchiveFileName() {
        return panel1.getArchiveFileName();
    }
}
