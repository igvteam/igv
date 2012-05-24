package org.broad.igv.hic;

import org.broad.igv.Globals;
import org.broad.igv.lists.GeneList;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.WaitCursorManager;
import org.broad.igv.util.LongRunningTask;
import org.broad.igv.util.NamedRunnable;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.net.Socket;
import java.net.UnknownHostException;
import java.util.Arrays;
import java.util.concurrent.*;

/**
 * @author Jim Robinson
 * @date 1/5/12
 */
public class IGVUtils {

    private static ExecutorService threadExecutor = Executors.newFixedThreadPool(1);

    static SocketHelper helper = null;

    static class SocketHelper {
        Socket socket = null;
        PrintWriter out = null;
        BufferedReader in = null;

        SocketHelper(BufferedReader in, PrintWriter out, Socket socket) {
            this.in = in;
            this.out = out;
            this.socket = socket;
        }
    }

    private static void createSocketHelper() {
        if (helper == null) {
            Socket socket = null;
            PrintWriter out = null;
            BufferedReader in = null;
            try {
                socket = new Socket("127.0.0.1", 60151);
                out = new PrintWriter(socket.getOutputStream(), true);
                in = new BufferedReader(new InputStreamReader(socket.getInputStream()));
                helper = new SocketHelper(in, out, socket);
            } catch (IOException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                helper = null;
            }

        }
    }


    /**
     * Send instructions to IGV to open or adjust views on the 2 loci.
     *
     * @param locus1
     * @param locus2
     */
    public static void sendToIGV(final String locus1, final String locus2) {

        Runnable runnable = new Runnable() {
            public void run() {

                if (locus1.startsWith(Globals.CHR_ALL) || locus1.startsWith("chrAll")) {
                    return;
                }

                // Same JVM?
                if (IGV.hasInstance()) {
                    IGV.getInstance().goToLociList(Arrays.asList(locus1, locus2));

                } else {

                    if (helper == null) createSocketHelper();

                    if (helper != null) {
                        String cmd = "gotoimmediate " + locus1 + " " + locus2;
                        helper.out.println(cmd);
                        String response = null;
                        try {
                            response = helper.in.readLine();
                        } catch (IOException e) {
                            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                        }
                        System.out.println(cmd + " " + response);
                    }
                }
            }
        };

        threadExecutor.submit(runnable);

    }

}
