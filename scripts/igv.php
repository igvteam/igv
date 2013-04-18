<?php
    $SESSION_URL=htmlspecialchars($_GET['sessionURL']);
    if($SESSION_URL == null || $SESSION_URL == "")
    {
      $SESSION_URL=htmlspecialchars($_GET['file']);
    }
    $LOCUS=$_GET['locus'];
    $MAXHEAP=$_GET['maxHeapSize'];
    $INITHEAP=$_GET['initialHeapSize'];
    $GENOME=$_GET['genome'];
    $INDEX=$_GET['index'];
    $NAME=$_GET['name'];

    header('Content-type: application/x-java-jnlp-file');
    header('Content-Disposition: attachment; filename="igv.jnlp"');
    header("Expires: Mon, 26 Jul 1997 05:00:00 GMT"); // Date in the past
    header("Last-Modified: " . gmdate("D, d M Y H:i:s") . " GMT");
    header("Cache-Control: no-store, no-cache, must-revalidate");
    header("Cache-control: post-check=0, pre-check=0, false");
    header("Pragma: no-cache");
    header("Content-Type: application/x-java-jnlp-file");    

?>
<?xml version="1.0" encoding="utf-8"?>

<jnlp
  spec="6.0+"
  codebase="http://www.broadinstitute.org/igv/projects/current">
  <information>
    <title>IGV 2.2</title>
    <vendor>The Broad Institute</vendor>
    <homepage href="http://www.broadinstitute.org/igv"/>
    <description>IGV Software</description>
    <description kind="short">IGV</description>
    <icon href="IGV_64.gif"/>
    <icon kind="splash" href="IGV_64.gif"/>
    <offline-allowed/>
	<shortcut/>
  </information>
  <security>
      <all-permissions/>
  </security>
  <update check="always" policy="always"/>
  <resources>
<?php 
   if($MAXHEAP == null || $MAXHEAP == "") {
     $MAXHEAP="900m";
   }
   if($INITHEAP == null || $INITHEAP == "") {
     $INITHEAP="256m";
   }   
   print('<java version="1.6+" initial-heap-size="');
   print($INITHEAP); 
   print('" max-heap-size="');
   print($MAXHEAP);
   print('"/>');
?>
    <jar href="igv.jar" download="eager" main="true"/>
    <jar href="batik-codec.jar" download="eager"/>
    <jar href="goby-io-igv.jar" download="lazy"/>   
    <property name="java.net.preferIPv4Stack" value="true"/> 
    <property name="apple.laf.useScreenMenuBar" value="true"/>
    <property name="com.apple.mrj.application.growbox.intrudes" value="false"/>
    <property name="com.apple.mrj.application.live-resize" value="true"/>
    <property name="com.apple.macos.smallTabs" value="true"/>
    <property name="production" value="true"/>
  </resources>
  <application-desc main-class="org.broad.igv.ui.Main">
<?php
    if($SESSION_URL != null && $SESSION_URL != "")
    {
        print("     <argument>$SESSION_URL</argument>\n");
    }
    if($LOCUS != null && $LOCUS != "")
    {
        print("     <argument>$LOCUS</argument>\n");
    }
    if($GENOME != null && $GENOME != "") 
    {
        print("    <argument>-g</argument>\n");
        print("    <argument>$GENOME</argument>\n");
    }
    if($INDEX != null && $INDEX != "") 
    {
        print("    <argument>-i</argument>\n");
        print("    <argument>$INDEX</argument>\n");
    }
    if($NAME != null && $NAME != "") 
    {
        print("    <argument>-n</argument>\n");
        print("    <argument>$NAME</argument>\n");
    }
    
?>
  </application-desc>
</jnlp>
