package org.igv.ui.util;

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
public class MagetabSignalDialog  extends org.igv.ui.IGVDialog
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

