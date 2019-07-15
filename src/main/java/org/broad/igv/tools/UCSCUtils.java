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
package org.broad.igv.tools;

import org.apache.log4j.Logger;

import java.io.*;

/**
 * @author jrobinso
 */
public class UCSCUtils {

    static private Logger log = Logger.getLogger(UCSCUtils.class);

    /**
     * 0 bin smallint(5) unsigned not null,	# bin scheme for indexing
     * 1  chrom varchar(255) not null,        # Human chromosome or FPC contig
     * 2  chromStart int unsigned not null,   # Start position in chromosome
     * 3  chromEnd int unsigned not null,     # End position in chromosome
     * 4  name varchar(255) not null, 	# Name of item
     * 5  span int unsigned not null, 	# each value spans this many bases
     * 6  count int unsigned not null,        # number of values in this block
     * 7  offset int unsigned not null,       # offset in File to fetch data
     * 8  file varchar(255) not null, 	# path name to .wib data file
     * 9  lowerLimit double not null, 	# lowest data value in this block
     * 10 dataRange double not null,  	# lowerLimit + dataRange = upperLimit
     * 11 validCount int unsigned not null,   # number of valid data values in this block
     * 12 sumData double not null,    	# sum of the data points, for average and stddev calculation
     * 13sumSquares double not null, 	# sum of data points squared, used for stddev calculation
     * 14 INDEX(chrom(8),bin)
     * );
     * <p/>
     * <p/>
     * 0       1       2       3       4       5       6       7        8                              9       10      11
     * 607     chr1    3000000 3005120 chr1.0  5       1024    0       /gbdb/mm8/wib/gc5Base.wib       0       100     1024    37240   1866400
     * <p/>
     * variableStep chrom=chr1 span=5
     * 3000001 19.685
     *
     * @param txtFile
     * @param wibFile
     * @param wigFile
     */
    public static void convertWIBFile(File txtFile, File wibFile, File wigFile, String trackLine) {

        BufferedReader txtReader = null;
        FileInputStream dataInputStream = null;
        PrintWriter wigWriter = null;

        String lastChr = "";
        int lastSpan = -1;
        try {
            txtReader = new BufferedReader(new FileReader(txtFile));
            dataInputStream = new FileInputStream(wibFile);
            wigWriter = new PrintWriter(new BufferedWriter(new FileWriter(wigFile)));

            if (trackLine != null) {
                wigWriter.println(trackLine);
            }

            String nextLine;
            while ((nextLine = txtReader.readLine()) != null) {
                String[] tokens = nextLine.split("\t");
                String chr = tokens[1];
                // .wig files are 1 based, database is 0 based
                int start = Integer.parseInt(tokens[2]) + 1;
                int span = Integer.parseInt(tokens[5]);

                if (!chr.equals(lastChr) || (span != lastSpan)) {
                    log.info("variableStep chrom=" + chr + " span=" + span);
                    wigWriter.println("variableStep chrom=" + chr + " span=" + span);
                    lastChr = chr;
                    lastSpan = span;
                }

                int count = Integer.parseInt(tokens[6]);
                byte[] data = new byte[count];

//seek to position "offset"
//read a block of "count" bytes into array: unsigned char data[]
//The first data byte corresponds to genome position "chromStart"
//The next data byte genome position is: "chromStart" + "span"
//data byte value of 128 indicates "NO DATA" at this position
//data byte values 0 to 127 are expanded to double types via:
                Integer offset = Integer.parseInt(tokens[7]);
                dataInputStream.getChannel().position(offset);
                DataInputStream dis = new DataInputStream(dataInputStream);
                dis.readFully(data);

                double lowerLimit = Double.parseDouble(tokens[9]);
                double dataRange = Double.parseDouble(tokens[10]);

                int chromPosition = start;
                for (int i = 0; i < count; ++i, chromPosition += span) {
                    if (data[i] < 128) {
                        double value = lowerLimit + (dataRange * ((double) data[i] / 127.0));
                        wigWriter.println(chromPosition + "\t" + value);
                    }
                }


            }

        } catch (IOException ex) {
            ex.printStackTrace();
        } finally {
            try {
                txtReader.close();
                dataInputStream.close();
                wigWriter.close();
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        }
    }

    public static void main(String[] args) {

        if (args.length < 3) {
            System.out.println("Usage convertWIGVile <txtFile> <wibFile> <wigFile> [trackLine]");
            System.exit(-1);
        }
        File txtFile = new File(args[0]);
        File wibFile = new File(args[1]);
        File wigFile = new File(args[2]);
        String trackLine = args.length > 3 ? args[3] : null;
        convertWIBFile(txtFile, wibFile, wigFile, trackLine);
    }
}
