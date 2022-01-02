package org.broad.igv.ext.render;

/*
 * Copyright 1999-2004 Carnegie Mellon University.
 * Portions Copyright 2002-2004 Sun Microsystems, Inc.
 * Portions Copyright 2002-2004 Mitsubishi Electric Research Laboratories.
 * All Rights Reserved.  Use is subject to license terms.
 *
 * See the file "README" for information on usage and
 * redistribution of this file, and for a DISCLAIMER OF ALL
 * WARRANTIES.
 *
 */

import java.util.Arrays;
import java.awt.Color;

/**
 * Color map representation - each entry in the map has an RGB value
 * associated with it
 *
 * @author Ron Weiss (ronw@ee.columbia.edu)
 *
 */
public class ColorMap
{
    public int size;
    public byte r[];
    public byte g[];
    public byte b[];
    public Color table[];

    /**
     * Create a color map that looks like Matlab's jet color map
     */
    public static ColorMap getJet()
    {
        return getJet(64);
    }

    /**
     * Create a color map witn n entries that looks like Matlab's jet
     * color map
     */
    public static ColorMap getJet(int n)
    {
        byte r[] = new byte[n];
        byte g[] = new byte[n];
        byte b[] = new byte[n];

        int maxval = 255;
        //Arrays.fill(g, 0, 8, (byte)0);
        Arrays.fill(g, 0, n/8, (byte)0);
        //for(int x = 0; x < 16; x++)
        //    g[x+8] = (byte)(maxval*x/16);
        for(int x = 0; x < n/4; x++)
            g[x+n/8] = (byte)(maxval*x*4/n);
        //Arrays.fill(g, 24, 40, (byte)maxval);
        Arrays.fill(g, n*3/8, n*5/8, (byte)maxval);
        //for(int x = 0; x < 16; x++)
        //    g[x+40] = (byte)(maxval-(maxval*x/16));
        for(int x = 0; x < n/4; x++)
            g[x+n*5/8] = (byte)(maxval-(maxval*x*4/n));
        //Arrays.fill(g, 56, 64, (byte)0);
        Arrays.fill(g, n*7/8, n, (byte)0);

        //for(int x = 0; x < g.length; x++)
        //    b[x] = g[(x+16) % g.length];
        for(int x = 0; x < g.length; x++)
            b[x] = g[(x+n/4) % g.length];
        //Arrays.fill(b, 56, 64, (byte)0);
        Arrays.fill(b, n*7/8, n, (byte)0);
        //Arrays.fill(g, 0, 8, (byte)0);
        Arrays.fill(g, 0, n/8, (byte)0);
        //for(int x = 8; x < g.length; x++)
        //    r[x] = g[(x+48) % g.length];
        for(int x = n/8; x < g.length; x++)
            r[x] = g[(x+n*6/8) % g.length];

        ColorMap cm = new ColorMap();
        cm.size = n;
        cm.r = r;
        cm.g = g;
        cm.b = b;
        cm.table = new Color[n];
        for(int x = 0; x < n; x++)
            //cm.table[x] = new Color((int)r[x]+maxval/2+1,
            //(int)g[x]+maxval/2+1, (int)b[x]+maxval/2+1);
            cm.table[x] = new Color(cm.getColor(x));
        return cm;
    }


    /**
     * Get the RGB value associated with an entry in this ColorMap
     */
    public int getColor(int idx)
    {
        int pixel = ((r[idx] << 16) & 0xff0000)
                | ((g[idx] << 8) & 0xff00)
                | (b[idx] & 0xff);

        return pixel;
    }

    public String toString()
    {
        StringBuffer s = new StringBuffer(500);
        for(int x = 0; x < size; x++)
        {
            s.append(x+": {"+r[x]+",\t"+g[x]+",\t"+b[x]+"}\t");
            if(x%3 == 2)
                s.append("\n");
        }

        return s.toString();
    }

    public static void main(String[] args)
    {
        ColorMap jet = getJet();
        ColorMap jet128 = getJet(128);


        System.out.println("Jet:\n"+jet+"\n\nJet128:\n"+jet128);
    }
}