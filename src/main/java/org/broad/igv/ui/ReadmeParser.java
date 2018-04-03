/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broad.igv.ui;

import org.broad.igv.util.FileUtils;

import java.io.*;
import java.util.regex.Pattern;

/**
 * Class used for examining our readme file, and presenting documentation to the user
 * either at the command line or elsewhere.
 * User: jacob
 * Date: 2012/02/27
 */
public class ReadmeParser {

    /**
     * When packaged in a jar, this is where the readme goes.
     */
    private static final String DEFAULT_README_PATH = "/docs/igvtools_readme.txt";
    
    public static final String END_STR = "See http://www.broadinstitute.org/software/igv/igvtools_commandline, or igvtools_readme.txt, for more help" + FileUtils.LINE_SEPARATOR;

    
    private String path = DEFAULT_README_PATH;

    public ReadmeParser(){}

    public ReadmeParser(String path){
        this.path = path;
    }

    Pattern getCmdPattern(String command){
        return Pattern.compile("command.{1,4}\"" + command, Pattern.CASE_INSENSITIVE);
    }
    
    Pattern getStartEndPattern(){
        return Pattern.compile("[-,=,_]{8}");
    }

    /**
     * Get reader for the given path. To get compliance 
     * in both source and jar form, first we try 
     * getResourceAsStream, then normal file handling.
     * @return
     */
    BufferedReader getReader() throws FileNotFoundException{
        InputStream is = ReadmeParser.class.getResourceAsStream(path);
        if(is != null){
            return new BufferedReader(new InputStreamReader(is));
        }else{
            return new BufferedReader(new FileReader(this.path));
        }

    }
    
    public String getDocForCommand(String cmd){
        String out = "Command " + cmd + " not found" + FileUtils.LINE_SEPARATOR;
        BufferedReader reader = null;
        try{
            reader = getReader();

            String line = null;
            Pattern cmdPattern = getCmdPattern(cmd);
            Pattern stEndPattern = getStartEndPattern();
            int startEndLinesRem = -1;

            while((line = reader.readLine()) != null){
                //Command found, set counter
                if(startEndLinesRem < 0 && cmdPattern.matcher(line).find()){
                    startEndLinesRem = 2;
                    out = "";
                }else if(startEndLinesRem < 0){
                    //No match for command, not within section
                    continue;
                }
                //Start/end line found, decrement number
                if(stEndPattern.matcher(line).find()){
                    startEndLinesRem--;
                }

                if(startEndLinesRem >= 0){
                    out += line + FileUtils.LINE_SEPARATOR;
                }

                if(startEndLinesRem == 0){
                    break;
                }

            }
        }catch(FileNotFoundException e){
            out = "Could not find readme file " + this.path + FileUtils.LINE_SEPARATOR;
        }catch (IOException e){
            out = "Error reading readme file " + this.path + FileUtils.LINE_SEPARATOR;
        }finally{
            if(reader != null){
                try{
                    reader.close();
                }catch(Exception e){
                    return null;
                }
            }
        }
        out += END_STR;
        return out;
    }
}
