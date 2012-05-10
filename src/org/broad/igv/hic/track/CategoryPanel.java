/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not
 * responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL), Version 2.1 which is
 * available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.hic.track;

import com.jidesoft.swing.ButtonStyle;
import com.jidesoft.swing.JideButton;
import org.broad.igv.util.ResourceLocator;

import javax.swing.*;
import javax.swing.border.Border;
import javax.swing.border.LineBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.geom.Path2D;
import java.awt.geom.Rectangle2D;
import java.awt.geom.RoundRectangle2D;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Set;

/**
 * @author Jim Robinson
 * @date 5/8/12
 */
public class CategoryPanel extends JPanel {

    boolean expanded;
    JideButton toggleButton;
    int nColumns = 5;
    private JPanel listPanel;
    private JPanel labelBar;

    public CategoryPanel(String name, List<ResourceLocator> locatorList, Set<String> loadedTrackNames) {

        expanded = true;

        setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
        setAlignmentX(LEFT_ALIGNMENT);
        //setLayout(null);

        labelBar = new JPanel();
        //labelBar.setBackground(Color.blue);
        labelBar.setLayout(new BoxLayout(labelBar, BoxLayout.X_AXIS));
        labelBar.setBorder(BorderFactory.createRaisedBevelBorder()); //  new LabelBorder(Color.black));
        labelBar.setAlignmentX(LEFT_ALIGNMENT);
        toggleButton = new JideButton(expanded ? " - " : " + ");
        toggleButton.setButtonStyle(ButtonStyle.HYPERLINK_STYLE);
        labelBar.add(toggleButton);

        labelBar.add(new JLabel(name));
        this.add(labelBar);


        listPanel = new JPanel();
        listPanel.setLayout(new GridLayout(0, 4));
        for (ResourceLocator loc : locatorList) {
            JCheckBox cb = new JCheckBox(loc.getTrackName());
            listPanel.add(cb);
        }
        this.add(listPanel);

        toggleButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent actionEvent) {
                expanded = !expanded;
                listPanel.setVisible(expanded);
            }
        });
        labelBar.addMouseListener(new MouseAdapter() {
            @Override
            public void mouseClicked(MouseEvent mouseEvent) {
                expanded = !expanded;
                listPanel.setVisible(expanded);
            }
        });

    }

    Collection<String> getSelectedTracks() {
        List<String> selectedTracks = new ArrayList<String>();
        for (Component c : listPanel.getComponents()) {
            if (c instanceof JCheckBox && ((JCheckBox) c).isSelected()) {
                selectedTracks.add(((JCheckBox) c).getText());
            }
        }
        return selectedTracks;
    }


    /**
     * If the <code>preferredSize</code> has been set to a
     * non-<code>null</code> value just returns it.
     * If the UI delegate's <code>getPreferredSize</code>
     * method returns a non <code>null</code> value then return that;
     * otherwise defer to the component's layout manager.
     *
     * @return the value of the <code>preferredSize</code> property
     * @see #setPreferredSize
     * @see javax.swing.plaf.ComponentUI
     */
    @Override
    public Dimension getPreferredSize() {
        if (listPanel == null) {
            return super.getPreferredSize();
        } else {

            Dimension d = listPanel.getPreferredSize();
            Component p = getRootPane();
            int h = listPanel.isVisible() ? d.height : 0;
            int w = p == null ? d.width : getParent().getWidth();
            return new Dimension(w, 3 + 3 + 30 + h);

        }
    }

    /**
     * If the minimum size has been set to a non-<code>null</code> value
     * just returns it.  If the UI delegate's <code>getMinimumSize</code>
     * method returns a non-<code>null</code> value then return that; otherwise
     * defer to the component's layout manager.
     *
     * @return the value of the <code>minimumSize</code> property
     * @see #setMinimumSize
     * @see javax.swing.plaf.ComponentUI
     */
    @Override
    public Dimension getMinimumSize() {
        return getPreferredSize();
    }

    @Override
    public void doLayout() {
        if (labelBar == null || listPanel == null) return;

        Dimension d = listPanel.getPreferredSize();
        Component p = getParent();
        int w = p == null ? d.width : getParent().getWidth();


        int y = 0;
        labelBar.setBounds(0, y, w, 30);
        y += 30;

        listPanel.setBounds(y, 33, w, d.height);

    }


}
