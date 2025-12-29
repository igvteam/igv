package org.broad.igv.session;

import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * Interface for a session persistable object.
 *
 * @author User: jrobinso
 * Date: Feb 23, 2010
 */
public interface Persistable {

    /**
     * Marshal object state in XML element
     * @return
     */

    public default void marshalXML(Document document, Element element) {};

    /**
     * Restore object state from an XML element
     */

    public default void unmarshalXML(Element element, Integer version) {};

}
