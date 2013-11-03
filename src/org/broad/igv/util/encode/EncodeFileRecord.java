package org.broad.igv.util.encode;

import java.util.Collection;
import java.util.Map;

/**
* @author jrobinso
*         Date: 10/31/13
*         Time: 10:11 PM
*/
public class EncodeFileRecord {

    boolean selected = false;
    String path;
    Map<String, String> attributes;

    EncodeFileRecord(String path, Map<String, String> attributes) {
        this.path = path;
        this.attributes = attributes;
    }

    public String getPath() {
        return path;
    }

    public String getFileType() {
        //String off trailing gz, if present
        String filetype = path;
        if (filetype.endsWith(".gz")) {
            filetype = filetype.substring(0, filetype.length() - 3);
        }
        int idx = filetype.lastIndexOf(".");
        return filetype.substring(idx + 1);
    }

    public String getAttributeValue(String name) {
        String value = attributes.get(name);
        if(name.equals("type") && value == null) value = getFileType();
        return value;
    }

    public Collection<String> getAttributeNames() {
        return attributes.keySet();
    }

    public boolean containsText(String filter) {
        for(String value : attributes.values()) {
            if(value.contains(filter)) return true;
        }
        return false;
    }

    boolean isSelected() {
        return selected;
    }

    void setSelected(boolean selected) {
        this.selected = selected;
    }
}
