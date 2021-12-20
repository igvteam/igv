package org.broad.igv.logging;

public class LogManager {

    public static Logger getLogger(Class classz) {
        return new Logger(classz.getSimpleName());
    }

    public static void shutdown() {
        LogFileHandler.getInstance().closeHandler();
    }
}
