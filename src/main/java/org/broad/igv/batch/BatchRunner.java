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

import org.broad.igv.logging.*;
import org.broad.igv.Globals;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.WaitCursorManager;
import org.broad.igv.ui.util.SnapshotUtilities;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.NamedRunnable;
import org.broad.igv.util.ParsingUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;

public class BatchRunner implements NamedRunnable {

    private static Logger log = LogManager.getLogger(BatchRunner.class);

    private String inputFile;
    private String rootPath;
    private IGV igv;

    public BatchRunner(final String inputFile, IGV igv) {
        this.inputFile = inputFile;
        this.igv = igv;

        // Find a root path for relative file paths
        if(FileUtils.isRemote(inputFile)) {
            rootPath = inputFile;
        } else {
            File f = new File(inputFile);
            if(f.exists()) {
                rootPath = f.getParent();
            }
        }
    }

    public String getName() {
        return "batchExecution";
    }

    public static void setIsBatchMode(boolean isBatchMode) {
        Globals.setSuppressMessages(isBatchMode);
        Globals.setBatch(isBatchMode);
    }

    public void run() {
        runWithDefaultGenome(null);
    }

    public void runWithDefaultGenome(String genomeId) {

        String inLine;
        setIsBatchMode(true);

        CommandExecutor cmdExe = new CommandExecutor(igv, rootPath);
        WaitCursorManager.CursorToken cursorToken = null;
        BufferedReader reader = null;
        try {
            cursorToken = WaitCursorManager.showWaitCursor();
            reader = ParsingUtils.openBufferedReader(inputFile);

            boolean firstCommand = true;
            while ((inLine = reader.readLine()) != null) {
                if (!(inLine.startsWith("#") || inLine.startsWith("//"))) {

                    if (firstCommand && genomeId != null && !inLine.toLowerCase().startsWith("genome")) {
                        log.debug("Loading genome " + genomeId);
                        GenomeManager.getInstance().loadGenomeById(genomeId);
                    }

                    log.debug("Executing Command: " + inLine);
                    String result = cmdExe.execute(inLine);
                    if(!result.equalsIgnoreCase("ok")) {
                        log.warn(inLine + " => " + result);
                    }
                    firstCommand = false;
                }
            }


        } catch (IOException ioe) {
            throw new DataLoadException(ioe.getMessage(), inputFile);
        } finally {
            setIsBatchMode(false);
            SnapshotUtilities.resetMaxPanelHeight();
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
}
