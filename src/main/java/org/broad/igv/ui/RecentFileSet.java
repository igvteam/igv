package org.broad.igv.ui;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;

public class RecentFileSet extends LimitedLinkedSet<String> {

    public RecentFileSet(int maxSize) {
        super(maxSize);
    }

    public RecentFileSet(Collection<String> c, int maxSize) {
        super(c, maxSize);
    }

    public String asString() {
        return String.join(";", this);
    }

    public static RecentFileSet fromString(String string, int maxSize) {
        if(string == null || string.isBlank()){
            return new RecentFileSet(maxSize);
        }
        String[] files = string.split(";");
        List<String> fileList = Arrays.stream(files)
                .filter(s -> !s.isBlank())
                // "null" was previously accounted for in older code so it's handled here
                // it doesn't seem like it should be possible to produce now though
                .filter(s -> !s.equals("null"))
                .map(String::strip)
                .toList();
        return new RecentFileSet(fileList, maxSize);
    }
}
