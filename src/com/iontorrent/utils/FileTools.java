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
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package com.iontorrent.utils;

import java.io.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JOptionPane;

/**
 *
 * @author Chantal Roth
 */
public class FileTools {

    public static String getFile(String title, String ext, String val) {
        return getFile(title, ext, val, false);
    }

    public static String getFile(String title, String ext, String val, boolean toSave) {
        JFileChooser cc = new JFileChooser();
        cc.setDialogTitle(title);


        if (val != null) {
            File f = new File(val);
            // if (!dir.isDirectory()) dir 
            cc.setSelectedFile(f);
            if (f.isDirectory()) {
                cc.setCurrentDirectory(f);
            } else if (f.getParentFile() != null) {
                cc.setCurrentDirectory(f.getParentFile());
            }
        }
        cc.setVisible(true);
        String[] Ext = new String[]{ext};
        if (ext.indexOf(",") > 0) {
            Ext = StringTools.parseList(ext, ", ").toArray(Ext);

        }
        ExtensionFileFilter filter1 = new ExtensionFileFilter(Ext[0] + " files", Ext);
        cc.setFileFilter(filter1);

        String res = val;
        int ans = 0;
        if (!toSave) {
            ans = cc.showOpenDialog(null);
        } else {
            ans = cc.showSaveDialog(null);
        }
        if (ans == JOptionPane.OK_OPTION) {
            File f = cc.getSelectedFile();
            if (toSave && f != null && ext != null && ext.length() > 0 && !ext.endsWith("*")) {
                String file = f.getAbsolutePath();
                if (ext.startsWith("*")) {
                    ext = ext.substring(1);
                }
                if (!file.endsWith(ext)) {
                    if (!ext.startsWith(".")) {
                        ext += ".";
                    }
                    p("Attaching " + ext + " to " + file);
                    file = file + ext;
                    f = new File(file);
                }
            }

            if (!f.exists()) {
                if (!toSave) {
                    JOptionPane.showMessageDialog(new JFrame(), f + " does not exist - please select an existing file");
                }
            } else {
                if (toSave) {

                    int ok = JOptionPane.showConfirmDialog(new JFrame(), "Would you want to overwrite file " + f + "?", "File exists", JOptionPane.YES_NO_OPTION);
                    if (ok != JOptionPane.YES_OPTION) {
                        return null;
                    }
                }
            }
            if (f.isDirectory()) {
                JOptionPane.showMessageDialog(null, f + " is a directory - please select a file");
                //f.get
            } else {
                res = f.getAbsolutePath();

            }

        } else {
            res = null;
        }
        return res;
    }
    public static boolean writeStringToFile(File f, String content, boolean append) {
        PrintWriter fout = null;
        try {
            fout = new PrintWriter(new BufferedWriter(new FileWriter(f, append)));
            fout.print(content);
            fout.flush();
            fout.close();
            return true;
        } catch (FileNotFoundException e) {
            err("File " + f + " not found");
        } catch (IOException e) {
            err("IO Exception");
        } finally {
            if (fout != null) {
                fout.flush();
                fout.close();
            }
        }
        return false;
    }

    private static void err(String msg, Exception ex) {
        Logger.getLogger(FileTools.class.getName()).log(Level.SEVERE, msg, ex);
    }

    private static void err(String msg) {
        Logger.getLogger(FileTools.class.getName()).log(Level.SEVERE, msg);
    }
    
    private static void p(String msg) {
        Logger.getLogger(FileTools.class.getName()).log(Level.INFO, msg);
    }
}
