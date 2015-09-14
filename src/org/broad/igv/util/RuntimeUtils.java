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
package org.broad.igv.util;

import org.apache.log4j.Logger;
import org.broad.igv.DirectoryManager;
import org.broad.igv.ui.util.MessageUtils;

import java.io.*;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLClassLoader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author jrobinso
 */
public class RuntimeUtils {

    private static Logger log = Logger.getLogger(RuntimeUtils.class);

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

    /**
     * Start an external process with the provided message.
     * Also starts a separate thread to read back error stream
     * <p/>
     * See {@link Runtime#exec(String, String[], java.io.File)} for explanation of arguments
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

            private boolean messageDisplayed = false;

            @Override
            public void run() {
                String line;
                try {
                    while ((line = err.readLine()) != null) {
                        log.error(line);
                        if (!messageDisplayed && line.toLowerCase().contains("error")) {
                            MessageUtils.showMessage(line + "<br>See igv.log for more details");
                            messageDisplayed = true;
                        }
                    }
                    err.close();
                } catch (IOException e) {
                    log.error(e.getMessage(), e);
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
     * @throws java.io.IOException
     * @deprecated Use {@link #executeShellCommand(String[], String[], java.io.File)}
     */
    @Deprecated
    public static String executeShellCommand(String cmd, String[] envp, File dir) throws IOException {
        return executeShellCommand(new String[]{cmd}, envp, dir);
    }


    public static String executeShellCommand(String cmd[], String[] envp, File dir) throws IOException {
        return executeShellCommand(cmd, envp, dir, true);
    }

    public static String executeShellCommand(String cmd[], String[] envp, File dir, boolean waitFor) throws IOException {
        Process pr = startExternalProcess(cmd, envp, dir);

        if(waitFor){
            try {
                pr.waitFor();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
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
            OutputStream os = pr.getOutputStream();
            if(os != null){
                os.close();
            }
        }
    }

    /**
     * Converts input files to file URLs
     * @param files
     * @return
     */
    private static URL[] filesToURLs(File[] files) {
        URL[] urls = new URL[files.length];
        for (int pp = 0; pp < files.length; pp++) {
            try {
                urls[pp] = new URL("file://" + files[0].getAbsolutePath());
            } catch (MalformedURLException e) {
                log.error(e);
            }
        }
        return urls;
    }

    /**
     * Returns an array of builtin library URL locations, as well as those passed in with {@code libURLs}
     * @param libURLs
     * @return
     */
    private static URL[] addBuiltinURLs(URL[] libURLs) {
        File[] builtin_paths = getPluginJars();
        URL[] builtin_urls = filesToURLs(builtin_paths);

        List<URL> outURLs = new ArrayList<URL>(Arrays.asList(libURLs != null ? libURLs : new URL[0]));
        outURLs.addAll(Arrays.asList(builtin_urls));

        return outURLs.toArray(new URL[outURLs.size()]);
    }

    private static URL[] getClassURLs(URL[] libURLs){
        return addBuiltinURLs(libURLs);
    }

    /**
     *  Create {@link java.lang.Class} object with the desired name,
     *  looking in {@code libURLs} as well as built-in locations
     * @param className
     * @param libURLs
     * @return
     * @throws IllegalAccessException
     * @throws InstantiationException
     * @throws ClassNotFoundException
     */
    public static Class loadClassForName(String className, URL[] libURLs) throws IllegalAccessException, InstantiationException, ClassNotFoundException {

        Class clazz = null;
        //Easy way
        try {
            clazz = Class.forName(className);
        } catch (ClassNotFoundException e) {
            //Try with custom loader below
        }
        if (clazz != null) return clazz;

        //If not found, check other locations
        URL[] allURLs = getClassURLs(libURLs);
        ClassLoader loader = URLClassLoader.newInstance(
                allURLs, ClassLoader.getSystemClassLoader()
        );
        clazz = loader.loadClass(className);
        return clazz;
    }

    /**
     * Create an instance of the specified class. Must have no-arg constructor
     * @param className
     * @param libURLs
     * @return
     * @throws IllegalAccessException
     * @throws InstantiationException
     * @throws ClassNotFoundException
     */
    public static Object loadInstanceForName(String className, URL[] libURLs) throws IllegalAccessException, InstantiationException, ClassNotFoundException, NoSuchMethodException, InvocationTargetException {

        Class clazz = loadClassForName(className, libURLs);
        Constructor constructor = clazz.getConstructor();
        return constructor.newInstance();
    }

    /**
     * Return an array of jar files in the plugin directory.
     * Never returns null, only empty array
     * @return
     */
    public static File[] getPluginJars(){
        File builtinDir = new File(DirectoryManager.getIgvDirectory(), "plugins/");
        File[] pluginJars =  builtinDir.listFiles(new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                return name.toLowerCase().endsWith(".jar");
            }
        });

        if(pluginJars == null){
            pluginJars = new File[0];
        }
        return pluginJars;
    }

    /**
     * Add jars in plugin directory to system classloader
     * This may not work on all systems, for security reasons
     */
    public static void loadPluginJars(){
        File[] pluginJars = getPluginJars();
        if(pluginJars != null){
            for(File jar: pluginJars){
                loadLibrary(jar);
            }
        }

    }

    //Not sure about the security implications of this method
    private static boolean loadLibrary(java.io.File jar){
        try {
            //We are using reflection here to circumvent encapsulation; addURL is not public
            java.net.URLClassLoader loader = (java.net.URLClassLoader)ClassLoader.getSystemClassLoader();
            java.net.URL url = jar.toURI().toURL();

            //Skip if already loaded
            for (java.net.URL it : java.util.Arrays.asList(loader.getURLs())){
                if (it.equals(url)){
                    return false;
                }
            }

            java.lang.reflect.Method method = java.net.URLClassLoader.class.getDeclaredMethod(
                    "addURL",
                    new Class[]{java.net.URL.class}
            );
            method.setAccessible(true);
            method.invoke(loader, new Object[]{url});
            return true;
        } catch (Exception e){
            log.warn("Error loading jar: " + e.getMessage());
            return false;
        }
    }


}
