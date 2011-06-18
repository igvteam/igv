/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.feature;

import org.broad.igv.util.IGVHttpUtils;

import java.io.DataInputStream;
import java.io.IOException;
import java.net.URL;
import java.net.URLConnection;

/**
 * @author jrobinso
 */
public class GetSequenceTest {

    static String sequenceURL = "http://www.broadinstitute.org/igv/sequence/hg18";


    public static void main(String[] args) {
        byte[] sequence = getSequence("hg18", "chr1", 1000000, 1000020);

        // Convert bytes to a string for convenience
        String sequenceString = new String(sequence);

        System.out.println(sequenceString);

        // Returned value should be  ACGTGGCTGCTCTCACACAT
    }


    public static byte[] getSequence(String genome, String chr, int start, int end) {


        try {

            URL url = new URL(sequenceURL + "?chr=" + chr + "&start=" + start + "&end=" + end);

            // Get the url content as an array of bytes
            byte[] sequence = new byte[end - start];
            DataInputStream dis = new DataInputStream(IGVHttpUtils.openConnectionStream(url));
            dis.readFully(sequence);
            dis.close();

            return sequence;

        } catch (IOException ex) {


            return null;
        }
    }

}
