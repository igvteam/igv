package org.broad.igv.feature;

/**
 * This class mainly exists so we can create a Proxy object of Exons
 * <p/>
 * User: jacob
 * Date: 2012/04/27
 */
public interface IExon extends IGVFeature {

    public int getCdStart();

    public int getCdEnd();

    boolean isNonCoding();
}
