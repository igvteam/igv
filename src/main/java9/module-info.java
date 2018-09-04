module org.igv {

    opens org.broad.igv.cli_plugin to java.xml.bind;
    opens org.broad.igv.data to java.xml.bind;
    opens org.broad.igv.feature.basepair to java.xml.bind;
    opens org.broad.igv.gwas to java.xml.bind;
    opens org.broad.igv.renderer to java.xml.bind;
    opens org.broad.igv.sam to java.xml.bind;
    opens org.broad.igv.session to java.xml.bind;
    opens org.broad.igv.tools.motiffinder to java.xml.bind;
    opens org.broad.igv.track to java.xml.bind;
    opens org.broad.igv.variant to java.xml.bind;
    
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
    requires java.xml.bind;
    requires javafx.base;
    requires javafx.controls;
    requires javafx.fxml;
    requires javafx.graphics;
    requires javafx.swing;
    requires jdk.xml.dom;
    requires jfreechart;
    requires log4j;
    requires org.apache.logging.log4j;
    requires org.apache.logging.log4j.core;
    requires swing.layout;
}
