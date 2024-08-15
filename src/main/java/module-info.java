module org.igv {
    exports org.broad.igv;
    exports org.broad.igv.tools;
    exports org.broad.igv.ui;
    exports org.broad.igv.event;
    exports org.broad.igv.jbrowse;
    exports org.broad.igv.logging;
    exports org.broad.igv.util.liftover;
    exports org.broad.igv.sam.smrt;
    exports org.broad.igv.ui.supdiagram;
    exports org.broad.igv.feature.genome.load to com.google.gson;
    opens org.broad.igv.feature.genome.load to com.google.gson;
    opens org.broad.igv.feature to com.google.gson;
    exports org.broad.igv.ucsc;

    requires com.google.common;
    requires commons.math3;
    requires com.google.gson;
    requires htsjdk;
    requires java.datatransfer;
    requires java.desktop;
    requires java.instrument;
    requires java.management;
    requires java.prefs;
    requires java.sql;
    requires java.xml;
    requires jdk.xml.dom;
    requires swing.layout;
    requires org.apache.commons.io;
    requires org.apache.commons.lang3;
    requires batik.codec;
    requires batik.svggen;
    requires batik.dom;
    requires AbsoluteLayout.RELEASE110;

    // AWS
    requires software.amazon.awssdk.core;
    requires software.amazon.awssdk.regions;
    requires software.amazon.awssdk.services.s3;
    requires software.amazon.awssdk.auth;
    requires software.amazon.awssdk.services.cognitoidentity;
    requires software.amazon.awssdk.services.sts;
    requires software.amazon.awssdk.http;
    requires software.amazon.awssdk.utils;

    requires jide.oss;
}
