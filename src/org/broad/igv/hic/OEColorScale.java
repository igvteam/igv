/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package org.broad.igv.hic;

import org.broad.igv.renderer.ContinuousColorScale;

import java.awt.*;

/**
 * @author Neva Cherniavsky
 * @date 3/22/12
 */
public class OEColorScale implements org.broad.igv.renderer.ColorScale {
    private int maxval = 5; // equivalent to maxdiv = 5 in Erez's Python code
    private int minval = 0;
    private int scheme = 1;

    public Color getColor(float score) {
        int red, green, blue;
        if (scheme == 1) {

            // why would score ever be < 0 or Nan??
            if (score < 0 || Float.isNaN(score)) {
                return new Color(0,0,0);
            }

         //   if (score < 0) {
         //       score = score*-1;
         //       red =  (int)(255*Math.min((Math.pow(score,-1))/(maxval),1));
         //       green = 0;
         //       blue = (int)(255*Math.min((score)/(maxval),1));
         //   }
         //   else {

                red = (int)(255*Math.min((score)/(maxval),1));
                green = 0;
                blue = (int)(255*Math.min((Math.pow(score,-1))/(maxval),1));
         //   }
            System.out.println(score + " " + red+" " + green + " " + blue );

            return new Color(red,green,blue);
        }
        else if (scheme == 0) {
            int R = 255;
            int G = (int)(Math.max(0,255-255*(score-minval)/(maxval-minval)));
            int B = (int)(Math.max(0,255-255*(score-minval)/(maxval-minval)));
            System.out.println(score + " " + R+" " + G+ " " + B );
            if (G < 0) G = 0;
            if (B < 0) B = 0;
            if (G > 255) G = 255;
            if (B > 255) B = 255;
            return new Color(R,G,B);
        }
        else return new Color(0,0,0) ;
    }

    public Color getColor(String symbol) {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }


    public Color getNoDataColor() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public String asString() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public boolean isDefault() {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }

}