package org.igv.util;

import org.igv.logging.*;
import org.igv.ui.util.MessageUtils;

import java.io.*;
import java.util.List;
import java.util.concurrent.Executors;
import java.util.function.Consumer;

/**
 * @author jrobinso
 */
public class RuntimeUtils {

    private static Logger log = LogManager.getLogger(RuntimeUtils.class);

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
