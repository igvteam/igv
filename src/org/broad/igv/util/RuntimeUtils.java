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

import com.google.common.primitives.Primitives;
import org.apache.log4j.Logger;
import org.broad.igv.DirectoryManager;

import java.io.*;
import java.lang.instrument.Instrumentation;
import java.lang.reflect.Array;
import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLClassLoader;
import java.util.*;

/**
 * @author jrobinso
 */
public class RuntimeUtils {

    private static Logger log = Logger.getLogger(RuntimeUtils.class);

    private static Instrumentation instrumentation;

    /**
     * Size, in BITS, of each primitive
     */
    private static Map<Class, Long> primitiveMemMap = new HashMap<Class, Long>(9);
    static{
        primitiveMemMap.put(Byte.TYPE, 8l);
        primitiveMemMap.put(Short.TYPE, 16l);
        primitiveMemMap.put(Integer.TYPE, 32l);
        primitiveMemMap.put(Long.TYPE, 64l);
        primitiveMemMap.put(Float.TYPE, 32l);
        primitiveMemMap.put(Double.TYPE, 64l);
        //May be an overestimate
        primitiveMemMap.put(Boolean.TYPE, 1l);
        primitiveMemMap.put(Character.TYPE, 16l);
        primitiveMemMap.put(Void.TYPE, 0l);
        assert primitiveMemMap.size() == Primitives.allPrimitiveTypes().size();
    }


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
     * Note: Not recursive
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

    static boolean isPrimitiveOrWrapper(Object o){
        return o.getClass().isPrimitive() || Primitives.isWrapperType(o.getClass());
    }

    public static long getObjectSizeRecursive(Object o, Set<Object> completedObjs){

        long fullSize = getObjectSize(o);
        completedObjs.add(o);

        //May be off a bit for wrapper classes,
        //but we let that go for now
        if(isPrimitiveOrWrapper(o)){
            return fullSize;
        }

        if(o.getClass().isArray()){
            for(int ii=0; ii < Array.getLength(o); ii++){
                Object el = Array.get(o, ii);
                fullSize += getObjectSizeRecursive(el, completedObjs);
            }
        }

        Set<Field> fields = new HashSet<Field>();
        getAllFields(o.getClass(), fields);

        for(Field field: fields){
            if(Modifier.isStatic(field.getModifiers())){
                continue;
            }

            field.setAccessible(true);

            //Field.get boxes primitives as wrapper classes,
            //which would artificially inflate the total
            Class fieldType = field.getType();
            if(fieldType.isPrimitive()){
                fullSize += primitiveMemMap.get(fieldType);
                continue;
            }

            try {
                Object fieldValue = field.get(o);
                fullSize += getObjectSizeRecursive(fieldValue, completedObjs);
            } catch (IllegalAccessException e) {
                throw new RuntimeException(e);
            }
        }
        return fullSize;
    }

    /**
     * Get all fields an object contains, including public, private, protected,
     * and inherited
     * @param clazz
     * @param fields
     */
    static void getAllFields(Class clazz, Set<Field> fields){
        fields.addAll(Arrays.asList(clazz.getDeclaredFields()));
        Class supClass = clazz.getSuperclass();
        if(supClass != null){
            getAllFields(supClass, fields);
        }
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
            OutputStream os = pr.getOutputStream();
            if(os != null){
                os.close();
            }
        }
    }

    private static URL[] getClassURLs() {
        String[] paths = new String[2];
        paths[0] = (new File(DirectoryManager.getIgvDirectory(), "plugins/")).getAbsolutePath();

        URL[] urls = new URL[paths.length];
        for (int pp = 0; pp < paths.length; pp++) {
            try {
                urls[pp] = new URL("file://" + paths[pp]);
            } catch (MalformedURLException e) {
                log.error(e);
            }
        }
        return urls;
    }

    public static Object loadClassForName(String className) throws IllegalAccessException, InstantiationException, ClassNotFoundException {

        Object object = null;
        //Easy way
        try {
            object = Class.forName(className).newInstance();
        } catch (ClassNotFoundException e) {
            //Try with custom loader below
        }
        if (object != null) return object;

        //If not found, check other locations
        ClassLoader loader = URLClassLoader.newInstance(
                getClassURLs(),
                ClassLoader.getSystemClassLoader()
        );
        return loader.loadClass(className);
    }


}
