package org.broad.igv.aws;

public class IGVS3Object {

    private String name;
    private String storageClass;

    private boolean isDir;

    public IGVS3Object(String name) {
        this.name = name;
        this.isDir = false;
        this.storageClass = "STANDARD";
    }

    public IGVS3Object(String name, boolean isDir, String storageClass) {
        this.name = name;
        this.isDir = isDir;
        this.storageClass = storageClass;
    }

    public String toString() {
        return name;
    }

    public String getName() {
        return name;
    }

    public String getStorageClass() { return storageClass; }

    public void setName(String name) {
        this.name = name;
    }

    public boolean isDir() {
        return isDir;
    }

    public void setDir(boolean dir) {
        isDir = dir;
    }
}
