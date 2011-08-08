package org.broad.igv.feature.genome;

/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: 8/8/11
 * Time: 8:18 AM
 * To change this template use File | Settings | File Templates.
 */
public interface Sequence {

    byte[] readSequence(String chr, int start, int end);
}
