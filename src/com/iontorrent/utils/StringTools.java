/*
*	Copyright (C) 2011 Life Technologies Inc.
*
*   This program is free software: you can redistribute it and/or modify
*   it under the terms of the GNU General Public License as published by
*   the Free Software Foundation, either version 2 of the License, or
*   (at your option) any later version.
*
*   This program is distributed in the hope that it will be useful,
*   but WITHOUT ANY WARRANTY; without even the implied warranty of
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*   GNU General Public License for more details.
*
*   You should have received a copy of the GNU General Public License
*   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package com.iontorrent.utils;

import java.util.ArrayList;
import java.util.StringTokenizer;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Chantal Roth
 */
public class StringTools {

    /** pattern can be 0,1,2 or 0-1,3,4-7 */
    public static ArrayList<Integer> parseInts(String pattern) {
        // check for -
        ArrayList<String> items = StringTools.parseList(pattern, ",");
        ArrayList<Integer> flows = new ArrayList<Integer>();
        for (String it : items) {
            if (it.indexOf("-") > 0) {
                ArrayList<Integer> nrs = StringTools.parseListtoInt(it, "-");
                if (nrs != null && nrs.size() == 2) {
                    int a = nrs.get(0);
                    int b = nrs.get(1);
                    for (int f = a; f <= b; f++) {
                        flows.add(f);
                    }
                }
            } else {
                try {
                    int f = Integer.parseInt(it);
                    flows.add(f);
                } catch (Exception e) {
                }
            }
        }
        return flows;
    }

    public static String addNL(String desc, String nl, int WIDTH) {
        if (desc == null || desc.length() < WIDTH) {
            return desc;
        }
        int oldnewlinepos = 0;
        int len = desc.length();
        StringBuffer newdesc = new StringBuffer(len + len / WIDTH + 1);
        int sp = desc.indexOf(" ");
        int oldsp = -1;
        int pos = 0;
        int newlinepos = 0;
        for (; pos < desc.length();) {
            int linelen = sp - oldnewlinepos;
            if (linelen > WIDTH + 5) {
                // ideally, cut at previous space
                newlinepos = oldsp;
                // if no other space before that, cut at next space
                if (newlinepos < oldnewlinepos) {
                    newlinepos = sp;
                }

                linelen = newlinepos - oldnewlinepos;

                if (linelen > 1.5 * WIDTH) {
                    // just cut in the middle of the word
                    newlinepos = oldnewlinepos + WIDTH;
                }
                newdesc = newdesc.append(desc.substring(oldnewlinepos, newlinepos));
                newdesc = newdesc.append(nl);
                oldnewlinepos = newlinepos;
                pos = newlinepos + 1;
            } else {
                pos = sp + 1;
            }
            oldsp = sp;
            if (sp + 1 > desc.length()) {
                break;
            }
            sp = desc.indexOf(" ", pos);

            if (sp < 0) {
                sp = desc.length() - 1;
            }

        }
        pos = Math.max(oldnewlinepos, newlinepos);
        while (len - pos > WIDTH) {
            newdesc = newdesc.append(desc.substring(pos, pos + WIDTH));
            newdesc = newdesc.append(nl);
            pos += WIDTH;
        }
        newdesc = newdesc.append(desc.substring(pos, len));
        return newdesc.toString();
    }

    public static String replace(String source, String tag, String with) {
        if (source == null || tag == null || tag.length() == 0 || with == null
                || tag.equals(with)) {
            return source;
        }
        if (tag.indexOf(with) >= 0) {
            String s = source;
            int tagpos = -1;
            while ((tagpos = s.indexOf(tag)) >= 0) {
                StringBuffer result = new StringBuffer();
                result.append(s.subSequence(0, tagpos));
                result.append(with);
                result.append(s.subSequence(tagpos + tag.length(), s.length()));
                s = result.toString();
            }
            return s;
        }
        StringBuffer result = new StringBuffer();
        int pos = 0;
        while (pos < source.length()) {
            int tagpos = source.indexOf(tag, pos);
            if (tagpos != -1) {
                if (tagpos > pos) {
                    result.append(source.substring(pos, tagpos));
                }
                result.append(with);
                pos = tagpos + tag.length();
            } else {
                result.append(source.substring(pos));
                break;
            }
        }
        return result.toString();
    }

    public static String getEnumeration(String[] list) {
        String result = "";
        if (list != null) {
            for (int i = 0; i < list.length; i++) {
                String element = list[i];
                result += element;
                if (i + 1 < list.length) {
                    result += ", ";
                }

            }
        }
        return result;
    }

    public static ArrayList<String> parseList(String list, String sep) {
        if (list == null) {
            return null;
        }

        list = list.trim();
        if (list.startsWith("[")) {
            list = list.substring(1);
        }
        if (list.endsWith("]")) {
            list = list.substring(0, list.length() - 1);
        }
        // out("input: "+list+", sep: "+sep);
        ArrayList<String> items = splitString(list, sep);

        return items;
    }

    public static ArrayList<Double> parseListToDouble(String list, String sep) {
        if (list == null) {
            return null;
        }

        list = list.trim();
        if (list.startsWith("[")) {
            list = list.substring(1);
        }
        if (list.endsWith("]")) {
            list = list.substring(0, list.length() - 1);
        }
        // p("input: "+list+", sep: "+sep);
        ArrayList<Double> res = new ArrayList<Double>();
        ArrayList<String> items = splitString(list, sep);
        for (String it : items) {
            double d = 0;        
            try {
                d= Double.parseDouble(it);
            }            
            catch (Exception e) {}
            res.add(d);
        }
        return res;
    }

    public static ArrayList<Long> parseListToLong(String list, String sep) {
        if (list == null) {
            return null;
        }

        list = list.trim();
        if (list.startsWith("[")) {
            list = list.substring(1);
        }
        if (list.endsWith("]")) {
            list = list.substring(0, list.length() - 1);
        }
        // p("input: "+list+", sep: "+sep);
        ArrayList<Long> res = new ArrayList<Long>();

        ArrayList<String> items = splitString(list, sep);
        for (String it : items) {
            long d = 0;        
            try {
                d= Long.parseLong(it);
            }            
            catch (Exception e) {}
        }
        return res;
    }

    public static ArrayList<Integer> parseListtoInt(String list, String sep) {
        return parseListtoInt(list, sep, Integer.MAX_VALUE);
    }

    public static ArrayList<Integer> parseListtoInt(String list, String sep, int max) {
        if (list == null) {
            return null;
        }

        list = list.trim();
        if (list.startsWith("[")) {
            list = list.substring(1);
        }
        if (list.endsWith("]")) {
            list = list.substring(0, list.length() - 1);
        }
        // p("input: "+list+", sep: "+sep);
        ArrayList<Integer> res = new ArrayList<Integer>();

        ArrayList<String> items = splitString(list, sep);
        for (String it : items) {
            int val = 0;        
            try {
                val= Integer.parseInt(it);
            }            
            catch (Exception e) {}
          
            if (val > max) {
                warn("Value too large:" + val + ", should be < " + max);
                val = max;
            }
            res.add(val);
        }
        return res;
    }

    public static double[] parseListTodouble(String list, String sep) {
        ArrayList<Double> res = parseListToDouble(list, sep);
        if (res == null) {
            return null;
        }
        double[] dd = new double[res.size()];
        for (int i = 0; i < res.size(); i++) {
            dd[i] = res.get(i);
        }
        return dd;
    }

    public static ArrayList<String> splitString(String line, String delim) {
        if (line == null) {
            err("splitString: line is null");
            return null;
        }
        ArrayList<String> res = new ArrayList<String>();
        StringTokenizer tokenizer = new StringTokenizer(line, delim);

        while (tokenizer.hasMoreElements()) {
            String next = tokenizer.nextToken().trim();
            if (next.startsWith("\"")) {
                next = next.substring(1, next.length());
            }
            if (next.endsWith("\"")) {
                next = next.substring(0, next.length() - 1);
            }
            res.add(next);
            // p("token "+res.size()+" is:"+next);
        }
        return res;
    }

    public static ArrayList<Long> splitStringToLongs(String line, String delim) {
        if (line == null) {
            err("splitString: line is null");
            return null;
        }
        ArrayList<Long> res = new ArrayList<Long>();
        StringTokenizer tokenizer = new StringTokenizer(line, delim);

        while (tokenizer.hasMoreElements()) {
            String next = tokenizer.nextToken().trim();
            if (next.startsWith("\"")) {
                next = next.substring(1, next.length());
            }
            if (next.endsWith("\"")) {
                next = next.substring(0, next.length() - 1);
            }
            long d = 0;        
            try {
                d= Long.parseLong(next);
            }            
            catch (Exception e) {}
            res.add(d);
            // p("token "+res.size()+" is:"+next);
        }
        return res;
    }

    public static String toStringList(ArrayList v) {
        StringBuffer sb = new StringBuffer();
        for (int i = 0; v != null && i < v.size(); i++) {
            String s = v.get(i).toString();
            sb = sb.append(s);
            if (i + 1 < v.size()) {
                sb = sb.append(";");
            }
        }
        return sb.toString();
    }

    public static String toStringList(String[] v) {
        StringBuffer sb = new StringBuffer();
        for (int i = 0; v != null && i < v.length; i++) {
            String s = v[i];
            sb = sb.append(s);
            if (i + 1 < v.length) {
                sb = sb.append(";");
            }
        }
        return sb.toString();
    }

    public static String[] toString(ArrayList v) {
        if (v == null || v.size() < 1) {
            return null;
        }
        String ss[] = new String[v.size()];
        for (int i = 0; v != null && i < v.size(); i++) {
            ss[i] = v.get(i).toString();
        }
        return ss;
    }

    public static ArrayList<String> toArrayList(String[] ss) {
        if (ss == null || ss.length < 1) {
            return null;
        }
        ArrayList<String> v = new ArrayList<String>(ss.length);
        for (int i = 0; ss != null && i < ss.length; i++) {
            v.add(ss[i]);
        }
        return v;
    }

    public static String getEnglishEnumeration(ArrayList<String> list) {
        String result = "";
        if (list != null) {
            for (int i = 0; i < list.size(); i++) {
                String element = (String) list.get(i);
                result += element;
                if (i + 2 < list.size()) {
                    result += ", ";
                }
                if (i + 2 == list.size()) {
                    result += " and ";
                }
            }
        }
        return result;
    }

    public static ArrayList<String> parseList(String list) {
        list = replace(list, " ", ",");
        return parseList(list, ",");
    }

    private void err(String msg, Exception ex) {
        Logger.getLogger(StringTools.class.getName()).log(Level.SEVERE, msg, ex);
    }
    private static void err(String msg) {
        Logger.getLogger(StringTools.class.getName()).log(Level.SEVERE, msg);
    }
    private static void warn(String msg) {
        Logger.getLogger(StringTools.class.getName()).log(Level.WARNING, msg);
    }
    private void p(String msg) {
       // System.out.println("StringTools: " + msg);
        //Logger.getLogger( StringTools.class.getName()).log(Level.INFO, msg, ex);
    }
}
