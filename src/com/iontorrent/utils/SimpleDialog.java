/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package com.iontorrent.utils;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.Toolkit;
import javax.swing.JDialog;
import javax.swing.JPanel;
import org.apache.log4j.Logger;

/**
 *
 * @author Chantal Roth
 */
public class SimpleDialog extends JDialog {
  
    private static Logger log = Logger.getLogger(SimpleDialog.class);
     
    
    
    public SimpleDialog(String title, JPanel mainpanel, int width, int height) {
        setLocationRelativeTo(null);
        this.setUndecorated(false);
        
        JPanel main = new JPanel(new BorderLayout());
        super.setTitle(title);
        this.add(main);
        main.add(mainpanel, BorderLayout.CENTER);
        Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
        
        int x = (int) Math.max(100, screen.getWidth() / 2 - 400);
        int y = (int) Math.max(100, screen.getHeight() / 2 - 200);
        this.setLocation(x, y);
        this.setVisible(true);
        this.setSize(width, height);
        this.toFront();       
    }
     private void p(String msg) {
        log.info(msg);
    }
}
