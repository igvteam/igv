package org.broad.igv.synteny;

/**
 * @author Jim Robinson
 * @date 3/28/12
 */
public interface Mapping {

    double mapPosition(int paramInt);

    boolean containsFromPosition(int fromPosition);

    int getFromStart();

    int getFromEnd();

    String getToChr();
}
