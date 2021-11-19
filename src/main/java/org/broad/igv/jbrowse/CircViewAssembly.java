package org.broad.igv.jbrowse;

class CircViewAssembly {
    String id;
    String name;
    CircViewRegion [] chromosomes;

    public CircViewAssembly(String id, String name, CircViewRegion[] regions) {
        this.id = id;
        this.name = name;
        this.chromosomes = regions;
    }
}
