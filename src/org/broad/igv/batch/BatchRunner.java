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

    public void run() {
        String inLine;
        CommandExecutor cmdExe = new CommandExecutor();
        Globals.setSuppressMessages(true);
        Globals.setBatch(true);

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
            Globals.setSuppressMessages(false);
            Globals.setBatch(false);
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
