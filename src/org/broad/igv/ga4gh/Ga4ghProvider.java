/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
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
