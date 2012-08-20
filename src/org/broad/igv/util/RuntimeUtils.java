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

import org.apache.log4j.Logger;

import java.io.*;
import java.lang.instrument.Instrumentation;

/**
 * @author jrobinso
 */
public class RuntimeUtils {

    private static Logger log = Logger.getLogger(RuntimeUtils.class);

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


    /**
     * Start an external process with the provided message.
     * Also starts a separate thread to read back error stream
     * <p/>
     * See {@link Runtime#exec(String, String[], File)} for explanation of arguments
     *
     * @return
     */
    public static Process startExternalProcess(String[] msg, String[] envp, File dir) throws IOException {
        Process pr = Runtime.getRuntime().exec(msg, envp, dir);
        startErrorReadingThread(pr);
        return pr;
    }

    private static Process startErrorReadingThread(Process pr) {
        final BufferedReader err = new BufferedReader(new InputStreamReader(pr.getErrorStream()));

        //Supposed to read error stream on separate thread to prevent blocking
        Thread runnable = new Thread() {
            @Override
            public void run() {
                String line;
                try {
                    while ((line = err.readLine()) != null) {
                        log.error(line);
                    }
                    err.close();
                } catch (IOException e) {
                    log.error(e);
                    e.printStackTrace();
                    throw new RuntimeException(e);
                }
            }
        };
        runnable.start();
        return pr;
    }

    /**
     * @param cmd
     * @param envp
     * @param dir
     * @return
     * @throws IOException
     * @deprecated Use {@link #executeShellCommand(String[], String[], File)}
     */
    @Deprecated
    public static String executeShellCommand(String cmd, String[] envp, File dir) throws IOException {
        return executeShellCommand(new String[]{cmd}, envp, dir);
    }


    public static String executeShellCommand(String cmd[], String[] envp, File dir) throws IOException {
        Process pr = startExternalProcess(cmd, envp, dir);

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
