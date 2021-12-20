package org.broad.igv.logging;

import org.broad.igv.DirectoryManager;

import java.io.File;
import java.io.IOException;
import java.util.logging.FileHandler;
import java.util.logging.Handler;
import java.util.logging.LogRecord;

public class LogFileHandler extends Handler {

    static LogFileHandler instance = new LogFileHandler();;
    private FileHandler currentHandler;

    private LogFileHandler() {

    }

    public static synchronized LogFileHandler getInstance() {
        return instance;
    }

    public void updateHandler() {
        if (this.currentHandler != null) {
            this.currentHandler.flush();
            this.currentHandler.close();
        }

        String logFilename = new File(DirectoryManager.getIgvDirectory(), "igv%g.log").getAbsolutePath();
        try {
            this.currentHandler = new FileHandler(logFilename, 100000, 2, true);
            this.currentHandler.setFormatter(new LogFormatter());

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void closeHandler() {
        currentHandler.flush();
        currentHandler.close();
        currentHandler = null;
    }

    @Override
    public void publish(LogRecord record) {
        if (currentHandler != null) {
            currentHandler.publish(record);
        }
    }

    @Override
    public void flush() {
        if (currentHandler != null) {
            currentHandler.flush();
        }
    }

    @Override
    public void close() throws SecurityException {
        if (currentHandler != null) {
            currentHandler.close();
        }
    }
}
