package org.broad.igv.aws;

import javax.swing.tree.DefaultMutableTreeNode;
import java.util.Collection;

public class S3TreeNode extends DefaultMutableTreeNode {


    public S3TreeNode(IGVS3Object userObject) {
        super(userObject);
     }

    public S3TreeNode(IGVS3Object userObject, boolean allowsChildren) {
        super(userObject, allowsChildren);
     }

    public void addS3Children(Collection<IGVS3Object> IGVS3Objects) {
        for (IGVS3Object IGVS3Object : IGVS3Objects) {
            this.add(new S3TreeNode(IGVS3Object));
        }
     }

    public boolean isLeaf() {
        return !((IGVS3Object) this.userObject).isDir();
    }

    public String toString() {
        return ((IGVS3Object) this.userObject).getName();
    }

    public IGVS3Object getUserObject() {
        return (IGVS3Object) super.getUserObject();
    }
}
