package org.broad.igv.batch;

import java.io.*;

/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: Apr 29, 2010
 * Time: 1:29:50 PM
 * To change this template use File | Settings | File Templates.
 */
public class Exec {

    public static void main(String args[]) throws IOException {

        Process p = Runtime.getRuntime().exec("cat");

        OutputStream os = p.getOutputStream();

        PrintWriter pw = new PrintWriter(new OutputStreamWriter(os));
        pw.println("abcd") ;
        pw.flush();
        os.close();

        InputStream is = p.getInputStream();
        
        int b;
        while (((b = is.read())) > 0) {
            System.out.print((char) b);
        }
    }
}
