package org.igv.ui;

import java.util.ArrayList;
import java.util.List;

/**
 * Class to hold a collection of messages for display to the user
 *
 * @author jrobinso
 */
public class MessageCollection {

    String header;

    List<String> messages = new ArrayList();


    public void setHeader(String header) {
        this.header = header;
    }

    public String getHeader() {
        return header;
    }

    public void append(String message) {
        messages.add(message);
    }

    public void prepend(String message) {
        messages.add(0, message);
    }

    public List<String> getMessages() {
        return messages;
    }

    public boolean isEmpty() {
        return messages.isEmpty();
    }

    public String getFormattedMessage() {
        StringBuffer buffer = new StringBuffer(80 * messages.size());

        for (String msg : messages) {
            buffer.append(msg);
            buffer.append("\n");
        }

        return buffer.toString();
    }
}
