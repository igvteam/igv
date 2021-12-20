package org.broad.igv.logging;

import org.broad.igv.DirectoryManager;
import java.io.File;
import java.io.IOException;
import java.util.logging.FileHandler;
import java.util.logging.Handler;
import java.util.logging.LogRecord;

/**
 * The purpose of this class it to enable moving IGV directories during a session and reinitializing logging
 * file handlers accordingly.  If the log file was in a stable location we could just use java.util.logging.FileHandler
 * directly.
 */

public class LogFileHandler extends Handler {

    static LogFileHandler instance = new LogFileHandler();;
    private FileHandler fileHandler;

    public static synchronized LogFileHandler getInstance() {
        return instance;
    }

    private LogFileHandler() {

    }

    /**
     * Update the fileHandler with the (possibly new) IGV directory.
     */
    public void updateHandler() {
        if (this.fileHandler != null) {
            this.fileHandler.flush();
            this.fileHandler.close();
        }

        String logFilename = new File(DirectoryManager.getIgvDirectory(), "igv%g.log").getAbsolutePath();
        try {
            this.fileHandler = new FileHandler(logFilename, 100000, 2, true);
            this.fileHandler.setFormatter(new LogFormatter());
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void closeHandler() {
        fileHandler.flush();
        fileHandler.close();
        fileHandler = null;
    }

    @Override
    public void publish(LogRecord record) {
        if (fileHandler != null) {
            fileHandler.publish(record);
        }
    }

    @Override
    public void flush() {
        if (fileHandler != null) {
            fileHandler.flush();
        }
    }

    @Override
    public void close() throws SecurityException {
        if (fileHandler != null) {
            fileHandler.close();
        }
    }
}
