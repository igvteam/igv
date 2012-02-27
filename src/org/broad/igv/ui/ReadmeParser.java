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
    
    private static final String DEFAULT_README_PATH = "igvtools_readme.txt";
    
    public static final String END_STR = "See http://www.broadinstitute.org/software/igv/igvtools_commandline, or igvtools_readme, for more help" + FileUtils.LINE_SEPARATOR;

    
    private String path;

    public ReadmeParser(){
        File userFile = new File(DEFAULT_README_PATH);
        if(userFile.exists()){
            this.path = userFile.getAbsolutePath();
            return;
        }
        
        File srcFile = new File("scripts/", DEFAULT_README_PATH);
        if(srcFile.exists()){
            this.path = srcFile.getAbsolutePath();
        }

        throw new IllegalStateException("Unable to find readme file");
    }

    public ReadmeParser(String path){
        this.path = path;
    }

    Pattern getCmdPattern(String command){
        return Pattern.compile("command.{1,10}\"" + command, Pattern.CASE_INSENSITIVE);
    }
    
    Pattern getStartEndPattern(){
        return Pattern.compile("[-,=,_]{4}");
    }
    
    public String getDocForCommand(String cmd){
        String out = "Command " + cmd + " not found";
        BufferedReader reader = null;
        try{
            reader = new BufferedReader(new FileReader(this.path));
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
