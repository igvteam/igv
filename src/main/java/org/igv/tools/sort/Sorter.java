package org.igv.tools.sort;

import java.io.File;
import java.io.IOException;
import java.util.Comparator;

/**
 * Created by jrobinson on 5/10/16.
 */
public interface Sorter {

    void run() throws IOException;

    void setComparator(Comparator<SortableRecord> readNameComparator);

    void setTmpDir(File tmpDir);

    void setMaxRecords(int maxRecords);

}
