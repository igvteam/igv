module org.igv {
    exports org.broad.igv;
    exports org.broad.igv.tools;
    exports org.broad.igv.ui;
    
    requires AbsoluteLayout;
    requires ant;
    requires com.google.common;
    requires commons.io;
    requires commons.math;
    requires goby.io.igv;
    requires gson;
    requires htsjdk;
    requires java.datatransfer;
    requires java.desktop;
    requires java.instrument;
    requires java.management;
    requires java.prefs;
    requires java.sql;
    requires java.xml;
    requires jdk.xml.dom;
    requires jfreechart;
    requires log4j;
    requires org.apache.logging.log4j;
    requires org.apache.logging.log4j.core;
    requires swing.layout;
}
