package org.igv.batch;

import org.igv.logging.*;
import org.igv.Globals;
import org.igv.exceptions.DataLoadException;
import org.igv.feature.genome.GenomeManager;
import org.igv.ui.IGV;
import org.igv.ui.WaitCursorManager;
import org.igv.ui.util.SnapshotUtilities;
import org.igv.util.FileUtils;
import org.igv.util.NamedRunnable;
import org.igv.util.ParsingUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;

public class BatchRunner  {

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

    public void run() throws IOException {
            runWithDefaultGenome(null);
    }

    public void runWithDefaultGenome(String genomeId) throws IOException {

        log.info("Executing batch script: " + inputFile);
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
        }  finally {
            setIsBatchMode(false);
            SnapshotUtilities.resetMaxPanelHeight();
            if (cursorToken != null) WaitCursorManager.removeWaitCursor(cursorToken);
            if (reader != null) {
                try {
                    reader.close();
                } catch (IOException e) {
                    log.error("Error closing reader", e);
                }
            }
        }
    }
}
