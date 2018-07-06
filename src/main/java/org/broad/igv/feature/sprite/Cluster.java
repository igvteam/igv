package org.broad.igv.feature.sprite;

import java.util.*;

/**
 * Created by jrobinso on 6/30/18.
 */
public class Cluster {

    String name;
    Map<String, List<Integer>> posMap;

    public Cluster(String name, Map<String, List<Integer>> posMap) {
        this.name = name;
        this.posMap = posMap;

    }


}
