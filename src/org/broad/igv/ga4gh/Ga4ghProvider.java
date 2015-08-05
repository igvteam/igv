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

package org.broad.igv.ga4gh;

import java.util.List;

/**
 * Created by jrobinso on 8/25/14.
 */
public class Ga4ghProvider {

    String name;
    String baseURL;
    String authKey;
    List<Ga4ghDataset> datasets;
    private Object id;

    public Ga4ghProvider(String name, String baseURL, String authKey, List<Ga4ghDataset> datasets) {
        this.name = name;
        this.baseURL = baseURL;
        this.authKey = authKey;
        this.datasets = datasets;
    }

    public String getName() {
        return name;
    }

    public String getBaseURL() {
        return baseURL;
    }

    public String getAuthKey() {
        return authKey;
    }

    public List<Ga4ghDataset> getDatasets() {
        return datasets;
    }

    public String toString() {
        return name;
    }

    public Object getId() {
        return id;
    }
}
