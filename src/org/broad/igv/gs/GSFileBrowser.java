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
 * Created by JFormDesigner on Sun Jun 05 19:25:45 EDT 2011
 */

package org.broad.igv.gs;

import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.gs.dm.DMUtils;
import org.broad.igv.gs.dm.GSDirectoryListing;
import org.broad.igv.gs.dm.GSFileMetadata;
import org.broad.igv.ui.util.MessageUtils;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.*;
import java.io.IOException;
import java.net.URL;
import java.util.List;

/**
 * @author Jim Robinson
 */
public class GSFileBrowser extends JDialog {

    private GSDirectoryListing dirListing;

    public enum Mode {OPEN, SAVE}

    private static Logger log = Logger.getLogger(GSFileBrowser.class);

    static ImageIcon folderIcon;
    static ImageIcon fileIcon;
    static GSFileMetadata selectedFile;

    Mode mode = Mode.OPEN;

    public static void main(String[] args) throws IOException {
        GSFileBrowser fb = (new GSFileBrowser(null, GSFileBrowser.Mode.SAVE));
        fb.setVisible(true);
        System.out.println(fb.getPath());
    }

    String userRootUrl = null;

    public GSFileBrowser(Frame owner) throws IOException {
        this(owner, Mode.OPEN);
    }

    public GSFileBrowser(Frame owner, Mode mode) throws IOException {
        super(owner);
        setModal(true);
        initComponents();
        init(mode);
        // For now single selection
        this.fileList.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        this.mode = mode;
        savePanel.setVisible(mode == Mode.SAVE);
        newFolderButton.setVisible(mode == Mode.SAVE);
        openButton.setText(mode == Mode.OPEN ? "Open" : "Save");
        getRootPane().setDefaultButton(openButton);


    }

    void init(Mode mode) throws IOException {

        if (folderIcon == null) {
            folderIcon = new ImageIcon(getClass().getResource("/images/Folder-icon.png"));
            fileIcon = new ImageIcon(getClass().getResource("/images/file-document.png"));
        }
        fileList.setCellRenderer(new CellRenderer());

        MouseListener mouseListener = new MouseAdapter() {
            public void mouseClicked(MouseEvent e) {
                int index = fileList.locationToIndex(e.getPoint());
                GSFileMetadata md = (GSFileMetadata) fileList.getModel().getElementAt(index);
                setSelectedFile(md);
                if (e.getClickCount() == 2) {
                    if (md.isDirectory()) {
                        try {
                            fetchContents(new URL(md.getUrl()));
                        } catch (IOException e1) {
                            e1.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                        }
                    } else {
                        setVisible(false);
                    }
                }
            }
        };
        fileList.addMouseListener(mouseListener);

        String rootdirectory = mode == Mode.OPEN ? DMUtils.DEFAULT_DIRECTORY : DMUtils.PERSONAL_DIRECTORY;
        URL defaultURL = new URL(PreferenceManager.getInstance().get(PreferenceManager.GENOME_SPACE_DM_SERVER) +
                rootdirectory);
        fetchContents(defaultURL);
    }

    private void setSelectedFile(GSFileMetadata md) {
        selectedFile = md;
        if (md.isDirectory()) {
            selectedFileTextField.setText(null);
        } else {
            selectedFileTextField.setText(md.getName());
        }
    }

    public String getFileURL() {
        return selectedFile == null ? null : selectedFile.getUrl();
    }

    public String getPath() {
        if (selectedFile == null) {
            return null;
        }
        if (selectedFile.isDirectory()) {
            if (mode == Mode.OPEN) {
                return null;
            } else {
                return selectedFile.getPath() + "/" + selectedFileTextField.getText();
            }
        } else {
            return selectedFile.getPath();
        }
    }


    private void fetchContents(URL url) throws IOException {

        dirListing = DMUtils.getDirectoryListing(url);
        String dirUrlString = dirListing.getDirectory().getUrl();

        setTitle(dirUrlString);
        if (userRootUrl == null) {
            userRootUrl = dirUrlString;
            // A little trick to get the root meta data
            int idx = userRootUrl.lastIndexOf("/");
            String user = userRootUrl.substring(idx).replace("/", "");
            GSFileMetadata rootMD = new GSFileMetadata(".", "/users/" + user, userRootUrl, "", "", true);
            setSelectedFile(rootMD);
        }

        List<GSFileMetadata> elements = dirListing.getContents();
        //Unless this is the root directory create a "up-one-level" entry
        if (!dirUrlString.equals(userRootUrl)) {
            int lastSlashIdx = dirUrlString.lastIndexOf("/");
            String parentURL = dirUrlString.substring(0, lastSlashIdx);
            elements.add(0, new GSFileMetadata("Parent Directory", "", parentURL, "", "", true));
        }

        ListModel model = new ListModel(dirListing.getContents());
        fileList.setModel(model);

    }

    private void cancelButtonActionPerformed(ActionEvent e) {
        selectedFile = null;
        setVisible(false);
        dispose();
    }

    private void loadButtonActionPerformed(ActionEvent e) {

        try {
            Object[] selections = fileList.getSelectedValues();

            // If there is a single selection, and it is a directory,  load that directory
            if (selections.length == 1) {
                GSFileMetadata md = (GSFileMetadata) selections[0];
                if (md.isDirectory()) {
                    fetchContents(new URL(md.getUrl()));
                    return;
                }
            }

            for (Object obj : selections) {
                if (obj instanceof GSFileMetadata) {
                    GSFileMetadata md = (GSFileMetadata) obj;
                    if (mode == Mode.OPEN) {
                        if (!md.isDirectory()) {
                            selectedFile = md;
                        }
                    }
                }
            }

            setVisible(false);
            dispose();

        } catch (Exception e1) {
            log.error("Error loading GS files", e1);
            MessageUtils.showMessage("Error: " + e1.toString());
        }
    }

    private void newFolderButtonActionPerformed(ActionEvent e) {
        String folderName = MessageUtils.showInputDialog("Name of new folder:");
        if (folderName != null && folderName.length() > 0) {
            String dirurl = selectedFile.getUrl();
            if (!selectedFile.isDirectory()) {
                // Strip off file part
                int idx = dirurl.lastIndexOf("/");
                dirurl = dirurl.substring(0, idx);
            }
            String putURL = dirurl + "/" + folderName;
            try {
                GSFileMetadata metaData = DMUtils.createDirectory(putURL);
                if(metaData != null) {
                    setSelectedFile(metaData);
                }
                // Refresh
                fetchContents(new URL(putURL));
            } catch (IOException e1) {
                log.error("Error creating directory: " + putURL, e1);
                MessageUtils.showMessage("<html>Error creating directory: " + e1 + "<br>" + e1.getMessage());
            }
        }

    }


    static class ListModel extends AbstractListModel {

        List<GSFileMetadata> elements;

        ListModel(List<GSFileMetadata> elements) {
            this.elements = elements;
        }

        public int getSize() {
            return elements.size();
        }

        public Object getElementAt(int i) {
            return elements.get(i);
        }

    }

    static class CellRenderer extends JLabel implements ListCellRenderer {
        // This is the only method defined by ListCellRenderer.
        // We just reconfigure the JLabel each time we're called.

        public Component getListCellRendererComponent(
                JList list,
                Object value,            // value to display
                int index,               // cell index
                boolean isSelected,      // is the cell selected
                boolean cellHasFocus)    // the list and the cell have the focus
        {
            GSFileMetadata fileElement = (GSFileMetadata) value;

            String s = value.toString();
            setText(s);
            setIcon(fileElement.isDirectory() ? folderIcon : fileIcon);
            if (isSelected) {
                setBackground(list.getSelectionBackground());
                setForeground(list.getSelectionForeground());
            } else {
                setBackground(list.getBackground());
                setForeground(list.getForeground());
            }
            setEnabled(list.isEnabled());
            setFont(list.getFont());
            setOpaque(true);
            return this;
        }
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        dialogPane = new JPanel();
        buttonBar = new JPanel();
        hSpacer2 = new JPanel(null);
        newFolderButton = new JButton();
        hSpacer1 = new JPanel(null);
        cancelButton = new JButton();
        openButton = new JButton();
        savePanel = new JPanel();
        label2 = new JLabel();
        selectedFileTextField = new JTextField();
        splitPane1 = new JPanel();
        scrollPane1 = new JScrollPane();
        fileList = new JList();
        label1 = new JLabel();

        //======== this ========
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== dialogPane ========
        {
            dialogPane.setBorder(new EmptyBorder(12, 12, 12, 12));
            dialogPane.setLayout(new BorderLayout());

            //======== buttonBar ========
            {
                buttonBar.setBorder(new EmptyBorder(12, 0, 0, 0));
                buttonBar.setLayout(new BoxLayout(buttonBar, BoxLayout.X_AXIS));
                buttonBar.add(hSpacer2);

                //---- newFolderButton ----
                newFolderButton.setText("New Folder");
                newFolderButton.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        newFolderButtonActionPerformed(e);
                    }
                });
                buttonBar.add(newFolderButton);
                buttonBar.add(hSpacer1);

                //---- cancelButton ----
                cancelButton.setText("Cancel");
                cancelButton.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        cancelButtonActionPerformed(e);
                    }
                });
                buttonBar.add(cancelButton);

                //---- openButton ----
                openButton.setText("Open");
                openButton.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        loadButtonActionPerformed(e);
                    }
                });
                buttonBar.add(openButton);
            }
            dialogPane.add(buttonBar, BorderLayout.SOUTH);

            //======== savePanel ========
            {
                savePanel.setVisible(false);
                savePanel.setLayout(new BoxLayout(savePanel, BoxLayout.X_AXIS));

                //---- label2 ----
                label2.setText("Save As: ");
                savePanel.add(label2);
                savePanel.add(selectedFileTextField);
            }
            dialogPane.add(savePanel, BorderLayout.NORTH);

            //======== splitPane1 ========
            {
                splitPane1.setLayout(new BoxLayout(splitPane1, BoxLayout.Y_AXIS));

                //======== scrollPane1 ========
                {
                    scrollPane1.setViewportView(fileList);
                }
                splitPane1.add(scrollPane1);

                //---- label1 ----
                label1.setHorizontalAlignment(SwingConstants.CENTER);
                label1.setIcon(new ImageIcon(getClass().getResource("/images/genomespacelogo.png")));
                splitPane1.add(label1);
            }
            dialogPane.add(splitPane1, BorderLayout.CENTER);
        }
        contentPane.add(dialogPane, BorderLayout.CENTER);
        setSize(490, 530);
        setLocationRelativeTo(getOwner());
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel dialogPane;
    private JPanel buttonBar;
    private JPanel hSpacer2;
    private JButton newFolderButton;
    private JPanel hSpacer1;
    private JButton cancelButton;
    private JButton openButton;
    private JPanel savePanel;
    private JLabel label2;
    private JTextField selectedFileTextField;
    private JPanel splitPane1;
    private JScrollPane scrollPane1;
    private JList fileList;
    private JLabel label1;
    // JFormDesigner - End of variables declaration  //GEN-END:variables
}
