/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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
 * Created by JFormDesigner on Sun Jun 05 19:25:45 EDT 2011
 */

package org.broad.igv.gs;

import java.awt.event.*;

import org.apache.log4j.Logger;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.IGVHttpUtils;
import org.broad.igv.util.LongRunningTask;
import org.broad.igv.util.NamedRunnable;
import org.broad.igv.util.ResourceLocator;
import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;
import org.json.JSONTokener;

import java.awt.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.BufferedInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.*;
import java.util.List;
import javax.swing.*;
import javax.swing.border.*;

/**
 * @author Stan Diamond
 */
public class GSFileBrowser extends JDialog {

    private static Logger log = Logger.getLogger(GSFileBrowser.class);

    static ImageIcon folderIcon;
    static ImageIcon fileIcon;

    public static void main(String[] args) throws IOException, JSONException {
        (new GSFileBrowser(null)).setVisible(true);
    }

    URL baseUrl;
    String userRootUrl = null;

    public GSFileBrowser(Frame owner) throws IOException, JSONException {
        super(owner);
        initComponents();
        baseUrl = new URL("https://dmtest.genomespace.org:8444/datamanager/defaultdirectory");
        init(baseUrl);
    }


    void init(URL url) throws JSONException, IOException {

        if (folderIcon == null) {
            folderIcon = new ImageIcon(getClass().getResource("/images/Folder-icon.png"));
            fileIcon = new ImageIcon(getClass().getResource("/images/file-document.png"));
        }
        fileList.setCellRenderer(new CellRenderer());

        MouseListener mouseListener = new MouseAdapter() {
            public void mouseClicked(MouseEvent e) {
                if (e.getClickCount() == 2) {
                    int index = fileList.locationToIndex(e.getPoint());
                    GSMetaData md = (GSMetaData) fileList.getModel().getElementAt(index);
                    if (md.isDirectory) {
                        try {
                            fetchContents(new URL(md.url));
                        } catch (IOException e1) {
                            e1.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                        } catch (JSONException e1) {
                            e1.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                        }
                    } else {
                        load(Arrays.asList(new ResourceLocator(md.url)));
                    }
                }
            }
        };
        fileList.addMouseListener(mouseListener);


        fetchContents(url);
    }

    private void fetchContents(URL url) throws IOException, JSONException {
        StringBuffer buf = new StringBuffer();
        InputStream is = null;
        try {
            is = IGVHttpUtils.openConnectionStream(url);
            BufferedInputStream bis = new BufferedInputStream(is);
            int b;
            while ((b = bis.read()) >= 0) {
                buf.append((char) b);
            }

            // System.out.println(buf.toString());

            JSONTokener tk = new JSONTokener(buf.toString());
            JSONObject obj = new JSONObject(tk);

            Object c = obj.get("contents");
            List<JSONObject> contents = new ArrayList();
            if (c instanceof JSONObject) {
                contents.add((JSONObject) c);
            } else {
                JSONArray tmp = (JSONArray) c;
                int l = tmp.length();
                for (int i = 0; i < l; i++) {
                    contents.add((JSONObject) tmp.get(i));
                }
            }

            JSONObject directory = (JSONObject) obj.get("directory");

            String dirUrlString = directory.get("url").toString();
            setTitle(dirUrlString);
            if (userRootUrl == null) {
                userRootUrl = dirUrlString;
            }


            ArrayList<GSMetaData> dirElements = new ArrayList();
            ArrayList<GSMetaData> fileElements = new ArrayList();

            //Unless this is the root directory create a "up-one-level" entry
            if (!dirUrlString.equals(userRootUrl)) {
                int lastSlashIdx = dirUrlString.lastIndexOf("/");
                String parentURL = dirUrlString.substring(0, lastSlashIdx);
                dirElements.add(new GSMetaData("Parent Directory", parentURL, "", "", true));
            }

            int contentsLength = contents.size();
            for (int i = 0; i < contentsLength; i++) {
                JSONObject o = contents.get(i);
                String name = (String) o.get("name");
                String objurl = (String) o.get("url");
                if (o.get("directory").equals("true")) {
                    dirElements.add(new GSMetaData(name, objurl, "", "", true));
                } else {
                    JSONObject dataFormat = o.has("dataFormat") ? (JSONObject) o.get("dataFormat") : null;
                    String format = dataFormat == null ? "" : dataFormat.getString("name");
                    String size = (String) o.get("size");
                    fileElements.add(new GSMetaData(name, objurl, format, size, false));
                }
            }
            ArrayList elements = new ArrayList(dirElements.size() + fileElements.size());
            elements.addAll(dirElements);
            elements.addAll(fileElements);

            ListModel model = new ListModel(elements);
            fileList.setModel(model);

        } finally {
            is.close();
        }
    }

    private void cancelButtonActionPerformed(ActionEvent e) {
        setVisible(false);
        dispose();
    }

    private void loadButtonActionPerformed(ActionEvent e) {

        try {
            Object[] selections = fileList.getSelectedValues();

            // If there is a single selection, and it is a directory,  load that directory
            if (selections.length == 1) {
                GSMetaData md = (GSMetaData) selections[0];
                if (md.isDirectory) {
                    fetchContents(new URL(md.url));
                    return;
                }
            }

            ArrayList<ResourceLocator> toLoad = new ArrayList(selections.length);
            for (Object obj : selections) {
                if (obj instanceof GSMetaData) {
                    GSMetaData md = (GSMetaData) obj;
                    if (!md.isDirectory) {
                        toLoad.add(new ResourceLocator(md.url));
                    }
                }
            }
            if (toLoad.size() > 0) {
                load(toLoad);
            }
        } catch (Exception e1) {
            log.error("Error loading GS files", e1);
            MessageUtils.showMessage("Error: " + e1.toString());
        }
    }

    private void load(final java.util.List<ResourceLocator> toLoad) {
        setVisible(false);
        dispose();
        NamedRunnable runnable = new NamedRunnable() {
            public void run() {
                IGV.getInstance().loadTracks(toLoad);
            }

            public String getName() {
                return "GS Load";
            }
        };

        LongRunningTask.submit(runnable);
    }


    static class GSMetaData {
        boolean isDirectory;
        String name;
        String url;
        String format;
        String size;

        GSMetaData(String name, String url, String format, String size, boolean isDirectory) {
            this.isDirectory = isDirectory;
            this.name = name;
            this.url = url;
            this.format = format;
            this.size = size;
        }

        public String toString() {
            return name;
        }
    }


    static class ListModel extends AbstractListModel {

        ArrayList elements;

        ListModel(ArrayList elements) {
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
            GSMetaData metaData = (GSMetaData) value;

            String s = value.toString();
            setText(s);
            setIcon(metaData.isDirectory ? folderIcon : fileIcon);
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
        cancelButton = new JButton();
        openButton = new JButton();
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
                buttonBar.setLayout(new GridBagLayout());
                ((GridBagLayout) buttonBar.getLayout()).columnWidths = new int[]{0, 85, 85, 0};
                ((GridBagLayout) buttonBar.getLayout()).columnWeights = new double[]{1.0, 0.0, 0.0, 0.0};

                //---- cancelButton ----
                cancelButton.setText("Cancel");
                cancelButton.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        cancelButtonActionPerformed(e);
                    }
                });
                buttonBar.add(cancelButton, new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0,
                        GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                        new Insets(0, 0, 0, 5), 0, 0));

                //---- openButton ----
                openButton.setText("Open");
                openButton.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        loadButtonActionPerformed(e);
                    }
                });
                buttonBar.add(openButton, new GridBagConstraints(3, 0, 1, 1, 0.0, 0.0,
                        GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                        new Insets(0, 0, 0, 0), 0, 0));
            }
            dialogPane.add(buttonBar, BorderLayout.SOUTH);

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
    private JButton cancelButton;
    private JButton openButton;
    private JPanel splitPane1;
    private JScrollPane scrollPane1;
    private JList fileList;
    private JLabel label1;
    // JFormDesigner - End of variables declaration  //GEN-END:variables
}
