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

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.ui.WaitCursorManager;
import org.broad.igv.util.NamedRunnable;
import org.broad.igv.util.ParsingUtils;

import java.io.BufferedReader;
import java.io.IOException;

public class BatchRunner implements NamedRunnable {
    private static Logger log = Logger.getLogger(BatchRunner.class);

    String inputFile;

    public BatchRunner(final String inputFile) {

        this.inputFile = inputFile;

    }

    public String getName() {
        return "batchExecution";  //To change body of implemented methods use File | Settings | File Templates.
    }

    public static void setIsBatchMode(boolean isBatchMode){
        Globals.setSuppressMessages(isBatchMode);
        Globals.setBatch(isBatchMode);
    }

    public void run() {
        String inLine;
        setIsBatchMode(true);

        CommandExecutor cmdExe = new CommandExecutor();

        WaitCursorManager.CursorToken cursorToken = null;
        BufferedReader reader = null;
        try {
            cursorToken = WaitCursorManager.showWaitCursor();
            reader = ParsingUtils.openBufferedReader(inputFile);

            while ((inLine = reader.readLine()) != null) {
                if (!(inLine.startsWith("#") || inLine.startsWith("//"))) {
                    log.info("Executing Command: " + inLine);
                    cmdExe.execute(inLine);
                }
            }


        } catch (IOException ioe) {
            throw new DataLoadException(ioe.getMessage(), inputFile);
        } finally {
            setIsBatchMode(false);
            if (cursorToken != null) WaitCursorManager.removeWaitCursor(cursorToken);
            if (reader != null) {
                try {
                    reader.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }

    /**
     * Returns true if this is "provably" a batch file.  Proof in this instance means the first line of the file
     * is #batch
     *
     * @param resource path to a file or URL
     */
    public static boolean isBatchFile(String resource) {

        BufferedReader reader = null;

        try {
            reader = ParsingUtils.openBufferedReader(resource);
            String firstLine = reader.readLine();
            return firstLine.startsWith("#batch");
        }
        catch (IOException e) {
            return false;
        }
        finally {
            if (reader != null) {
                try {
                    reader.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }

    }

}
