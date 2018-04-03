/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
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
