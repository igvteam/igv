package org.igv.session;

import java.io.IOException;
import java.io.InputStream;

/**
 * @author Jim Robinson
 * @date 1/12/12
 */
public interface SessionReader {

    static boolean isSessionFile(String f) {
        return f.endsWith(".xml") || f.endsWith(".php") || f.endsWith(".php3") || f.endsWith(".session");
    }

    void loadSession(InputStream inputStream, Session session, String sessionPath) throws IOException;
}
