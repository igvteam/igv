package org.broad.igv.lists;

import java.util.LinkedHashMap;
import java.util.List;

/**
 * @author jrobinso
 * @date Dec 15, 2010
 */
public class GeneListGroup {

    private String name;
    private boolean editable = true;
    private LinkedHashMap<String, GeneList> geneLists = new LinkedHashMap();


    public GeneListGroup(String name, List<GeneList> lists) {
        this.name = name;
        for(GeneList gl : lists) {
            geneLists.put(gl.getName(), gl);
        }
    }

    public boolean isEditable() {
        return editable;
    }

    public void setEditable(boolean editable) {
        this.editable = editable;
    }

    public String getName() {
        return name;
    }

    public GeneList getGeneList(String listID) {
        return geneLists.get(listID);
    }

    public LinkedHashMap<String, GeneList> getGeneLists() {
        return geneLists;
    }
}
