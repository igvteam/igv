package org.broad.igv.ui.genome;

import java.util.HashMap;
import java.util.Map;

/**
 * A record representing a genome hosted on a remote server, or more specifically a row in the IGV hosted genome table.
 */
public class GenomeListItem {

    String id;
    String path;
    String displayableName;

    Map<String, String> attributes;

    public GenomeListItem(String displayableName, String path, String id, Map<String, String> attributes) {
        this.id = id;
        this.path = path;
        this.displayableName = displayableName;
        this.attributes = attributes;
    }

    /**
     * Alternate constructor for legacy genome list format.
     *
     * @param displayableName The name that can be shown to a user.
     * @param path            The location of the genome archive, can be a file path or URL
     * @param id              The id of the genome.
     */
    public GenomeListItem(String displayableName, String path, String id) {
        this.id = id;
        this.path = path;
        this.displayableName = displayableName;
        attributes = new HashMap<>();
        attributes.put("common name", displayableName);
        attributes.put("url", path);
        attributes.put("assembly", id);
    }

    public String getAttributeValue(String name) {
        return attributes.get(name);
    }

    public String getId() {
        return id;
    }

    public String getPath() {
        return path;
    }

    public String getDisplayableName() {
        return displayableName;
    }

    public String toString() {
        return displayableName;
    }

@Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof GenomeListItem)) return false;
        GenomeListItem that = (GenomeListItem) o;
        return id != null && id.equals(that.id) &&
               path != null && path.equals(that.path);
    }

    @Override
    public int hashCode() {
        int result = id != null ? id.hashCode() : 0;
        result = 31 * result + (path != null ? path.hashCode() : 0);
        return result;
    }
}
