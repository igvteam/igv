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
