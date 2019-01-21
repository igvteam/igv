package org.broad.igv.aws;

import javax.swing.tree.DefaultMutableTreeNode;
import java.util.Collection;

public class S3TreeNode extends DefaultMutableTreeNode {


    public S3TreeNode(S3Object userObject) {
        super(userObject);
     }

    public S3TreeNode(S3Object userObject, boolean allowsChildren) {
        super(userObject, allowsChildren);
     }

    public void addS3Children(Collection<S3Object> s3Objects) {
        for (S3Object s3Object : s3Objects) {
            this.add(new S3TreeNode(s3Object));
        }
     }

    public boolean isLeaf() {
        return !((S3Object) this.userObject).isDir();
    }

    public String toString() {
        return ((S3Object) this.userObject).getName();
    }

    public S3Object getUserObject() {
        return (S3Object) super.getUserObject();
    }
}
