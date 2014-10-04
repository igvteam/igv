/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
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
