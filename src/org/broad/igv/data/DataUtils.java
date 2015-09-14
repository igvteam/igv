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
package org.broad.igv.data;

//~--- JDK imports ------------------------------------------------------------

import java.io.*;

/**
 * @author jrobinso
 */
public class DataUtils {

    public static int getIndexBefore(int[] values, int x) {
        return getIndexBefore(values, x, 0, values.length);
    }

    public static int getIndexBefore(int[] values, int x, int leftBound, int rightBound) {
        int idx = (leftBound + rightBound) / 2;

        if ((idx == 0) || (idx == values.length - 1)) {
            return idx;
        }
        if (values[idx] == x) {
            return idx;
        }

        if (values[idx] < x) {
            if (values[idx + 1] >= x) {
                return idx;
            } else {
                leftBound = idx;
                return getIndexBefore(values, x, leftBound, rightBound);
            }
        } else {    // values[idx] > x
            if (values[idx - 1] <= x) {
                return idx - 1;
            } else {
                rightBound = idx;
                return getIndexBefore(values, x, leftBound, rightBound);
            }
        }
    }


    /**
     * Estimate the number of rows in an ascii data file.  Estimate is based on
     * the first 100 lines, and assumes the line length is approximately
     * constant.
     *
     * @param textFile
     * @return
     */
    public static AsciiFileMetrics estimateFileMetrics(String textFile) {
        int estRowCount = 0;
        try {

            BufferedReader reader = null;

            File file = new File(textFile);

            reader = new BufferedReader(new FileReader(file));
            String nextLine = reader.readLine();

            double lineCount = 0;
            double nChars = 0;
            while ((nextLine = reader.readLine()) != null && (lineCount < 100)) {
                nChars += nextLine.length();
                lineCount++;
            }

            int columnCount = nextLine.split("\t").length;

            double charsPerLine = ((lineCount > 0) ? nChars / lineCount : 0);
            estRowCount = (int) (file.length() / charsPerLine);
            return new AsciiFileMetrics(estRowCount, columnCount, charsPerLine);

        }
        catch (FileNotFoundException ex) {
            ex.printStackTrace();

        }
        catch (IOException ex) {
            ex.printStackTrace();

        }
        return null;
    }

    /**
     * Method description
     *
     * @param fileMetrics
     * @return
     */
    public static int estimatePreprocessingTime(AsciiFileMetrics fileMetrics) {

        return 8 + (int) ((0.0036 * fileMetrics.getEstRowCount() * fileMetrics.getColumnCount()) / 100);
    }

    /**
     * This class has some useful metrics for optimizing reading of large ascii files
     *
     * @author jrobinso
     */
    public static class AsciiFileMetrics {

        private int estRowCount;

        private int columnCount;

        private double estBytesPerLine;

        public AsciiFileMetrics(int estRowCount, int columnCount, double estBytesPerLine) {
            this.estRowCount = estRowCount;
            this.columnCount = columnCount;
            this.estBytesPerLine = estBytesPerLine;
        }


        public double getEstBytesPerLine() {
            return estBytesPerLine;
        }


        public int getEstRowCount() {
            return estRowCount;
        }


        public int getColumnCount() {
            return columnCount;
        }
    }
}
