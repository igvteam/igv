/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.util;

import java.io.*;
import java.lang.instrument.Instrumentation;

/**
 * @author jrobinso
 */
public class RuntimeUtils {

    private static Instrumentation instrumentation;

    public static long getAvailableMemory() {
        Runtime runtime = Runtime.getRuntime();
        long maxMemory = runtime.maxMemory();
        long allocatedMemory = runtime.totalMemory();
        long freeMemory = runtime.freeMemory();
        return freeMemory + (maxMemory - allocatedMemory);
    }

    public static double getAvailableMemoryFraction() {
        Runtime runtime = Runtime.getRuntime();
        long maxMemory = runtime.maxMemory();
        long allocatedMemory = runtime.totalMemory();
        long freeMemory = runtime.freeMemory();
        return (double) ((freeMemory + (maxMemory - allocatedMemory))) / maxMemory;

    }

    public static void premain(String args, Instrumentation inst) {
        instrumentation = inst;
    }

    /**
     * Must invoke with instrumentation on command line
     * Note: Not recursive, so of somewhat limited utility
     *
     * @param o
     * @return
     */
    public static long getObjectSize(Object o) {
        if (instrumentation == null) {
            throw new IllegalStateException("No instrumentation available. Need to launch width -javaagent:path/to/RuntimeUtils.jar");
        }
        return instrumentation.getObjectSize(o);
    }


    public static String executeShellCommand(String cmd, String[] envp, File dir) throws IOException {

        Runtime run = Runtime.getRuntime();
        Process pr = null;
        try {
            pr = run.exec(cmd, envp, dir);
        } catch (IOException e) {
            e.printStackTrace();
        }
        try {
            pr.waitFor();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        InputStream inputStream = null;
        String line = "";
        try {
            inputStream = pr.getInputStream();
            BufferedReader buf = new BufferedReader(new InputStreamReader(inputStream));
            StringWriter writer = new StringWriter();
            PrintWriter pw = new PrintWriter(writer);
            while ((line = buf.readLine()) != null) {
                pw.println(line);
            }
            pw.close();
            return writer.toString();
        } finally {
            if (inputStream != null) {
                inputStream.close();

            }
        }
    }


}
