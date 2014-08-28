<?php

    $SESSION_URL=$_GET['sessionURL'];
    if($SESSION_URL == null || $SESSION_URL == "")
    {
      $SESSION_URL=$_GET['file'];
    }
    $LOCUS=$_GET['locus'];
    $MAXHEAP=$_GET['maxHeapSize'];
    $INITHEAP=$_GET['initialHeapSize'];
    $GENOME=$_GET['genome'];
    $INDEX=$_GET['index'];
    $NAME=$_GET['name'];
    $PREFS=$_GET['prefs'];
    $DATA_FORMAT=$_GET['dataFormat'];
    $USER_AGENT=$_SERVER['HTTP_USER_AGENT'];
    $iPod = stripos($_SERVER['HTTP_USER_AGENT'],"iPod");
    $iPhone = stripos($_SERVER['HTTP_USER_AGENT'],"iPhone");
    $iPad = stripos($_SERVER['HTTP_USER_AGENT'],"iPad");

    if($iPod || $iPhone || $iPad) {
       $IPAD_URL = "igvipadapp://eval";
       if($SESSION_URL != null && $SESSION_URL != "") {
          $IPAD_URL = $IPAD_URL . "?file=". $SESSION_URL;
          if($LOCUS != null && $LOCUS != "") {
            $IPAD_URL = $IPAD_URL . "&locus=". $LOCUS;
          }
          if($DATA_FORMAT != null && $DATA_FORMAT != "") {
            $IPAD_URL = $IPAD_URL . "&dataFormat=" . $DATA_FORMAT;
          }
          if($INDEX != null && $INDEX != "") 
          {
            $IPAD_URL = $IPAD_URL . "&index=" . $INDEX;
          }
        }
        else  if($LOCUS != null && $LOCUS != "") {
          $IPAD_URL = $IPAD_URL . "?locus=". $LOCUS;
         }
 
        header( 'Location: ' .$IPAD_URL);
    }


    $igv_project = htmlspecialchars($_GET['user']);   
    if($igv_project == null) {
      $igv_project = "launcher";
    }
    $igv_version = "2.3";
    $referred_by = $_SERVER['HTTP_USER_AGENT'];
    $user_ip = $_SERVER['HTTP_X_FORWARDED_FOR'];
    if ($user_ip == '') {
        $user_ip = $_SERVER['REMOTE_ADDR'];
    }

 
    
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
  codebase="http://igv.broadinstitute.org/app/current">
  <information>
    <title>IGV 2.3</title>
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
  <update check="background"/>
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
    <jar href="igv.jar" download="eager" main="true" version="2.3.33"/>
    <jar href="batik-codec.jar" download="eager" version="1.7"/>
    <jar href="goby-io-igv.jar" download="lazy" version="1.0"/>  
    <property name="jnlp.versionEnabled" value="true"/>
    <property name="java.net.preferIPv4Stack" value="true"/> 
    <property name="apple.laf.useScreenMenuBar" value="true"/>
    <property name="com.apple.mrj.application.growbox.intrudes" value="false"/>
    <property name="com.apple.mrj.application.live-resize" value="true"/>
    <property name="com.apple.macos.smallTabs" value="true"/>
    <property name="http.agent" value="IGV"/>
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
    
    print("    <argument>--preferences</argument>\n");
    if($PREFS != null && $PREFS != "") 
    {
        print("    <argument>$PREFS</argument>\n");
    }
    else {
        print("   <argument>http://www.broadinstitute.org/igv/projects/current/genomespace.properties</argument>\n");
    }
    
?>
  </application-desc>
</jnlp>
