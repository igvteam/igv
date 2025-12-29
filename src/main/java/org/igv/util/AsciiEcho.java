/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.igv.util;

import java.io.IOException;

/**
 * @author jrobinso
 */
public class AsciiEcho {

    public static void main(String[] args) throws IOException {

        while (true) {
            System.out.println((int) System.in.read());
        }

    }
}
