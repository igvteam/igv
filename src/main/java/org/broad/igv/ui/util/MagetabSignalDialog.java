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

