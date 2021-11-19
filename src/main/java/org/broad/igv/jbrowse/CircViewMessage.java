package org.broad.igv.jbrowse;

class CircViewMessage {
    String message;
    Object data;
    public CircViewMessage(String message, Object data) {
        this.message = message;
        this.data = data;
    }
}
