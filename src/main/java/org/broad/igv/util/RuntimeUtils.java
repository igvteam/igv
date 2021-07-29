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
import org.broad.igv.ui.util.MessageUtils;

import java.io.*;
import java.util.List;
import java.util.concurrent.Executors;
import java.util.function.Consumer;

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

    public static String exec(String command) throws IOException, InterruptedException {

        List<String> commandArgs = StringUtils.breakQuotedString(command, ' ');

        Process process;
        if(commandArgs.size() == 1) {
             process = Runtime.getRuntime().exec(command);
        } else {
            String [] args = commandArgs.toArray(new String []{});
            process = Runtime.getRuntime().exec(args);
        }

        StringBuilder results = new StringBuilder();
        BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
        String line;
        while ((line = reader.readLine()) != null) {
            results.append(line);
            results.append('\n');
        }
        reader.close();

        StringBuilder errors = new StringBuilder();
        BufferedReader errorReader = new BufferedReader(
                new InputStreamReader(process.getErrorStream()));
        while ((line = errorReader.readLine()) != null) {
            errors.append(line);
            errors.append('\n');
        }
        errorReader.close();

        int exitValue = process.waitFor();
        if (exitValue != 0) {
            throw new RuntimeException(errorReader.toString());
        }
        process.destroy();

        return results.toString();
    }

    private static class StreamGobbler implements Runnable {
        private InputStream inputStream;
        private Consumer<String> consumer;

        public StreamGobbler(InputStream inputStream, Consumer<String> consumer) {
            this.inputStream = inputStream;
            this.consumer = consumer;
        }

        @Override
        public void run() {
            new BufferedReader(new InputStreamReader(inputStream)).lines()
                    .forEach(consumer);
        }
    }
}
