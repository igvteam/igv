package org.broad.igv.logging;

import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.List;
import java.util.function.Supplier;
import java.util.logging.ConsoleHandler;

/**
 * Emulates some aspects of the log4j Logger API (the aspects IGV actually uses).
 */
public class Logger {


    java.util.logging.Logger wrappedLogger;

    /**
     * Instantiate a logger and add console and file handlers
     *
     * @param name
     */
    public Logger(String name) {

        wrappedLogger = java.util.logging.Logger.getLogger(name);
        wrappedLogger.setUseParentHandlers(false);

        ConsoleHandler handler = new ConsoleHandler();
        handler.setFormatter(new LogFormatter());
        wrappedLogger.addHandler(handler);
        wrappedLogger.addHandler(LogFileHandler.getInstance());

    }

    public void info(String message) {
        wrappedLogger.info(message);
    }

    public void error(String message) {
        wrappedLogger.severe(message);
    }

    public void error(String message, Throwable e) {
        wrappedLogger.severe(message);
        logThrowable(e);
    }

    public void error(Throwable e) {
        logThrowable(e);
    }

    public void log(Level level, String message) {

        switch (level) {
            case ERROR:
                error(message);
                break;
            case WARN:
                warn(message);
                break;
            default:
                info(message);
        }
    }

    public void warn(String message) {
        wrappedLogger.warning(message);
    }

    public void warn(String message, Throwable e1) {
        wrappedLogger.warning(message);
        wrappedLogger.warning(e1.toString());
    }


    public void debug(String message) {
        if (isDebugEnabled()) {
            wrappedLogger.fine(message);
        }
    }

    public void debug(List<String> message) {
        if (isDebugEnabled()) {
            StringBuffer buf = new StringBuffer();
            for (String m : message) {
                if (buf.length() > 0) buf.append(", ");
                buf.append(m);
            }
            wrappedLogger.fine(buf.toString());
        }
    }

    public void trace(Object message) {
        if(isTraceEnabled()) {
            wrappedLogger.finest(message.toString());
        }
    }

    // TODO -- implement if needed
    public boolean isDebugEnabled() {
        return false;
    }

    // TODO -- implement if needed
    public boolean isTraceEnabled() {
        return false;
    }

    private void logThrowable(Throwable e) {
        Supplier<String> s = () -> {
            StringWriter writer = new StringWriter();
            PrintWriter pw = new PrintWriter(writer);
            e.printStackTrace(pw);
            pw.flush();
            return writer.toString();
        };
        wrappedLogger.severe(s);
    }
}
