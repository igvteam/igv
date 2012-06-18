/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package com.iontorrent.views;

import com.iontorrent.data.FlowDistribution;
import com.iontorrent.utils.SimpleDialog;

/**
 *
 * @author Chantal Roth
 */
public class FlowSignalDistributionView {
    
    /** contains the chart and actions for configure and save */
   FlowSignalDistributionPanel distributionPanel;
   
   public static void showView(FlowDistribution distributions[], String info) {
       FlowSignalDistributionPanel distributionPanel = new FlowSignalDistributionPanel(distributions, info);
       SimpleDialog dia = new SimpleDialog("Flow Signal Distribution", distributionPanel);      
   }
  
}
