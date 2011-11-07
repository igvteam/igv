package org.broad.igv;

import java.io.*;
import java.net.Socket;
import java.net.UnknownHostException;
import java.util.Arrays;
import java.util.Comparator;

/**
 * @author Jim Robinson
 * @date 11/3/11
 */
public class Launcher {

    public static void main(String[] args) throws IOException {

        String file = null;
        String locus = null;
        String genome = null;
        String user = null;
        String memory = null;
        String index = null;
        boolean newSession = false;

        for (String a : args) {
            String tokens[] = a.split("=");
            if (tokens.length == 2) {
                String key = tokens[0].toLowerCase();
                if (key.equals("file")) {
                    file = tokens[1];
                } else if (key.equals("locus")) {
                    locus = tokens[1];
                } else if (key.equals("genome")) {
                    genome = tokens[1];
                } else if (key.equals("user")) {
                    user = tokens[1];
                } else if (key.equals("maxheapsize")) {
                    memory = tokens[1];
                } else if (key.equals("index")) {
                    index = tokens[1];
                } else if (key.equals("newsession")) {
                    newSession = tokens[1].equals("true");
                }
            }
        }


        // TODO -- read port from igv preferences
        int port = 60151;
        boolean igvIsRunning = loadDirectly(port, file, locus, genome, newSession);

        if (!igvIsRunning) {

            // Launch a new IGV

            StringBuffer buf = new StringBuffer("http://www.broadinstitute.org/igv/projects/dev/igv.php");

            boolean firstArg = true;
            if (file != null) {
                String delim = firstArg ? "?" : "&";
                buf.append(delim);
                buf.append("sessionURL=");
                buf.append(file);
                firstArg = false;
            }
            if (locus != null) {
                String delim = firstArg ? "?" : "&";
                buf.append(delim);
                buf.append("locus=");
                buf.append(locus);
                firstArg = false;

            }
            if (genome != null) {
                String delim = firstArg ? "?" : "&";
                buf.append(delim);
                buf.append("genome=");
                buf.append(genome);
                firstArg = false;

            }
            if (user != null) {

                String delim = firstArg ? "?" : "&";
                buf.append(delim);
                buf.append("user=");
                buf.append(user);
                firstArg = false;
            }

            File jnlpFile = createJNLP(file, locus, genome, memory, index);
            System.out.println(jnlpFile.getAbsolutePath());

            ProcessBuilder pb = new ProcessBuilder("javaws", jnlpFile.getAbsolutePath());
            Process p = pb.start();
            // TODO -- read from stderr and report any errors to user


        }

        // Cleanup jnlp files
        File desktop = new File(System.getProperty("user.home") + "/Desktop");
        if (desktop.exists() && desktop.isDirectory()) {
            cleanupJnlpFiles(desktop);
        }
        File downloads = new File(System.getProperty("user.home") + "/Downloads");
        if (downloads.exists() && downloads.isDirectory()) {
            cleanupJnlpFiles(downloads);
        }


        System.exit(1);
    }


    private static boolean loadDirectly(int port, String file, String locus, String genome, boolean newSession) throws IOException {
        boolean success;
        Socket socket = null;
        PrintWriter out = null;
        BufferedReader in = null;
        try {
            socket = new Socket("127.0.0.1", port);
            out = new PrintWriter(socket.getOutputStream(), true);
            in = new BufferedReader(new InputStreamReader(socket.getInputStream()));

            if (genome != null) {
                out.println("genome " + genome);
                String response = in.readLine();
            }
            if (file != null) {
                String cmd = "load " + file;
                if(newSession) cmd += " clearSession";
                out.println(cmd);
                String response = in.readLine();
            }
            if (locus != null) {
                out.println("goto " + locus);
                String response = in.readLine();
            }
            success = true;

        } catch (UnknownHostException e) {
            success = false;
        } catch (IOException e) {
            success = false;
        } finally {
            if (in != null) in.close();
            if (out != null) out.close();
            if (socket != null) socket.close();
        }
        return success;
    }

    private static File createJNLP(String file, String locus, String genome, String memory, String index) throws IOException {

        String tmp = System.getProperty("java.io.tmpdir");
        if (tmp == null) tmp = ".";
        File tmpDir = new File(tmp);

        File f = new File(tmpDir, "igv" + System.currentTimeMillis() + ".jnlp");

        InputStream is = null;
        PrintWriter pw = null;

        String maxMem = memory == null ? "1050m" : memory;
        try {
            is = Launcher.class.getResourceAsStream("jnlpTemplate.txt");
            BufferedReader br = new BufferedReader(new InputStreamReader(is));
            pw = new PrintWriter(new BufferedWriter(new FileWriter(f)));

            String nextLine;
            while ((nextLine = br.readLine()) != null) {
                if (nextLine.contains("$ARGUMENTS")) {
                    if (file != null) pw.println("<argument>" + file + "</argument>");
                    if (locus != null) pw.println("<argument>" + locus + "</argument>");
                    if (genome != null) pw.println("<argument>-g</argument><argument>" + genome + "</argument>");
                    if (index != null) pw.println("<argument>-i</argument><argument>" + index + "</argument>");
                } else {
                    pw.println(nextLine.replace("$MAXHEAP", maxMem));
                }
            }
        } finally {
            if (is != null) is.close();
            if (pw != null) pw.close();
        }

        return f;
    }


    /**
     * Cleanup extra jnlp files.  This method is written specifcally for Mac OS.
     */
    public static void cleanupJnlpFiles(File dir) {

        if (dir.exists() && dir.isDirectory()) {
            File[] jnlpFiles = dir.listFiles(new FileFilter() {

                public boolean accept(File arg0) {
                    return arg0.getName().startsWith("igvLaunch") && arg0.getName().endsWith(".jnlp");
                }
            });

            // Delete all but the highest version (newest) jnlp file
            for (File f : jnlpFiles) {
                System.out.println("Delete " + f.getAbsolutePath());
                f.delete();
            }


        }
    }
}
