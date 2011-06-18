/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
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

package org.broad.igv.ui.util;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

/**
 * Created by IntelliJ IDEA.
 * User: nazaire
 * Date: Apr 21, 2011
 * Time: 11:06:07 AM
 * To change this template use File | Settings | File Templates.
 */
public class MagetabSignalDialog  extends JDialog
{
    private String quantitationColumn = null;
    private boolean canceled = true;

    public MagetabSignalDialog(Frame parent, String[] elements)
    {
        super(parent, true);


        JPanel panel = new JPanel(new GridLayout(0, 1));
        JLabel messageLabel = new JLabel();
        messageLabel.setText("Please select a column to obtain signal values from");
        panel.add(messageLabel);
        
        JButton btnCancel = new JButton("Cancel");
        btnCancel.addActionListener(new ActionListener()
        {
            public void actionPerformed(ActionEvent e)
            {
                canceled = true;
                quantitationColumn = null;
                dispose();
            }
        });

        ButtonGroup group = new ButtonGroup();

        JButton btnOk = new JButton("OK");
        btnOk.addActionListener(new ActionListener()
        {
            public void actionPerformed(ActionEvent e)
            {
                canceled = false;
                dispose();
            }
        });

        for(int i =0; i < elements.length; i++)
        {
            JRadioButton b = new JRadioButton(elements[i]);
            b.addActionListener(new ActionListener()
            {
                public void actionPerformed(ActionEvent event)
                {
                    JRadioButton btn = (JRadioButton)event.getSource();
                    quantitationColumn = btn.getText();
                }
            });
            group.add(b);
            panel.add(b);
        }
        panel.setBorder(BorderFactory.createEmptyBorder(3,8,4,8));
        JPanel bp = new JPanel();
        bp.add(btnOk);
        bp.add(btnCancel);

        JScrollPane scrollPane = new JScrollPane(panel);
        getContentPane().add(scrollPane, BorderLayout.CENTER);
        getContentPane().add(bp, BorderLayout.PAGE_END);

        pack();
        setResizable(false);
        setLocationRelativeTo(parent);
    }

    public String getQuantitationColumn()
    {
        return quantitationColumn;
    }

    public boolean isCanceled()
    {
        return canceled;
    }
}

