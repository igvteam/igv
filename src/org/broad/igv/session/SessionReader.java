package org.broad.igv.session;

import java.io.IOException;
import java.io.InputStream;

/**
 * @author Jim Robinson
 * @date 1/12/12
 */
public interface SessionReader {
    void loadSession(InputStream inputStream, Session session, String sessionName) throws IOException;
}
