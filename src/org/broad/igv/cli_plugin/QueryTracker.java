/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.cli_plugin;

import java.util.HashMap;
import java.util.Map;
import java.util.UUID;

/**
 * When there are multiple outputs, we don't want to re-run the command
 * each time. So we store what the last query was, and if it's the same
 * just use the old file.
 *
 * @author jacob
 * @date 2013-May-02
 */
public class QueryTracker {

    static Map<String, QueryTracker> trackerMap = new HashMap();

    String chr;
    int start;
    int end;
    int zoom;

    private QueryTracker(){}

    public void setLastQuery(String chr, int start, int end, int zoom){
        this.chr = chr;
        this.start = start;
        this.end = end;
        this.zoom = zoom;
    }

    public boolean isQuerySame(String chr, int start, int end, int zoom){
        return zoom == this.zoom && end == this.end && start == this.start && chr.equals(this.chr);
    }

    public static QueryTracker get() {
        return get(UUID.randomUUID().toString());
    }

    /**
     * Get a QueryTracker with the provided id.
     * One will be created if it doesn't exist.
     * @param id
     * @return
     */
    public static QueryTracker get(String id){
        if(trackerMap.containsKey(id)){
            return trackerMap.get(id);
        }else{
            QueryTracker qt = new QueryTracker();
            trackerMap.put(id, qt);
            return qt;
        }

    }
}
