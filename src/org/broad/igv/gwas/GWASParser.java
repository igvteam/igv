/*
 * Copyright (c) 2007-2009 by The Broad Institute, Inc. and the Massachusetts Institute of Technology.
 * All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL), Version 2.1 which
 * is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR WARRANTIES OF
 * ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT
 * OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR
 * RESPECTIVE TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES OF
 * ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES, ECONOMIC
 * DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER THE BROAD OR MIT SHALL
 * BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE
 * FOREGOING.
 */
package org.broad.igv.gwas;

import org.apache.log4j.Logger;
import org.broad.igv.exceptions.ParserException;
//import org.broad.igv.feature.ParsingUtils;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.util.ParsingUtils;
//import org.broad.igv.util.AsciiLineReader;
import org.broad.tribble.readers.AsciiLineReader;

import org.broad.igv.util.ResourceLocator;

import java.io.FileInputStream;
import java.io.IOException;
import static java.lang.Math.log10;

/**
 * Parses GWAS PLINK result files
 *
 * @author paananen
 */
public class GWASParser {

    private static Logger log = Logger.getLogger(GWASParser.class);
    private ResourceLocator locator;
    private final int indexSize = 10000;

    private int locationCol = -1;
    private int chrCol = -1;
    private int pCol = -1;
    private int SNPCol = -1;


    public ResourceLocator getLocator() {
        return locator;
    }

    public void setLocator(ResourceLocator locator) {
        this.locator = locator;
    }


    public GWASParser(ResourceLocator locator) {
        this.locator = locator;

    }


    public String parse(long startingBytes, int row) throws IOException {
        FileInputStream fs = null;

        AsciiLineReader reader = null;
        String nextLine = null;

        String returnString = "";
        Genome genome = GenomeManager.getInstance().getCurrentGenome();


        try {
            fs = new FileInputStream(locator.getPath());
            fs.getChannel().position(startingBytes);
            reader = new AsciiLineReader(fs);

            int rowCounter = 0;

            while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {

                nextLine = nextLine.trim();

                String[] tokens = new String[100];
                int size = ParsingUtils.splitSpaces(nextLine, tokens);


                if (tokens.length > 1) {
                    //log.info(nextLine);

                    if (rowCounter >= indexSize - 1)
                        break;

                    


                    //String chr = ParsingUtils.convertChrString(tokens[chrCol].trim());

                    //if (genome != null) {
                    //   String chr = genome.getChromosomeAlias(tokens[chrCol].trim());
                    // }

                    int start;

                    try {
                        start = Integer.parseInt(tokens[locationCol].trim());
                    } catch (NumberFormatException e) {
                        throw new ParserException("Column " + locationCol + " must be a numeric value.", reader.getCurrentLineNumber(), nextLine);
                    }

                    // TODO: Select transformation (-log10, ln, etc.) from preferences
                    // Check if the p-value is NA
                    if (!tokens[pCol].trim().equals("NA")) {
                        float p;

                        try {
                            p = Float.parseFloat(tokens[pCol].trim());
                            // Transform to -log10
                            p = (float) -log10((double) p);


                        } catch (NumberFormatException e) {
                            throw new ParserException("Column " + pCol + " must be a numeric value.", reader.getCurrentLineNumber(), nextLine);
                        }

                        rowCounter++;

                    }

                }

            }

            return returnString;
        }

        catch (
                ParserException e
                )

        {
            throw e;
        }

        catch (
                Exception e
                )

        {
            if (nextLine != null && reader.getCurrentLineNumber() != 0) {
                throw new ParserException(e.getMessage(), e, reader.getCurrentLineNumber(), nextLine);
            } else {
                throw new RuntimeException(e);
            }
        }

        finally

        {
            reader.close();
            fs.close();
        }

    }



    public GWASData parse() throws IOException {


          FileInputStream fs = null;

          AsciiLineReader reader = null;
          String nextLine = null;
          //int byteCounter = 0;
         Genome genome = GenomeManager.getInstance().getCurrentGenome();


          try {
              //reader = ParsingUtils.openAsciiReader(locator);
              fs = new FileInputStream(locator.getPath());
              fs.getChannel().position(0);
              reader = new AsciiLineReader(fs);

              String headerLine = reader.readLine();
              //byteCounter += headerLine.getBytes().length;
              headerLine = headerLine.trim();

              String[] headers = new String[100];
              int headersSize = ParsingUtils.splitSpaces(headerLine, headers);
              if (headersSize < 4)
                  throw new ParserException("Incorrect amount of column headers.", reader.getCurrentLineNumber(), nextLine);

              int colCounter = 0;
              for (int i = 0; i < headersSize; i++) {
                  String header = headers[i];
                  header = header.toLowerCase();
                  if (header.equals("chr") || header.equals("chromosome"))
                      chrCol = colCounter;
                  if (header.equals("bp") || header.equals("pos") || header.equals("position"))
                      locationCol = colCounter;
                  if (header.equals("p") || header.equals("pval") || header.equals("p-value") || header.equals("pvalue"))
                      pCol = colCounter;

                  if (header.equals("snp") || header.equals("rs") || header.equals("rsid") || header.equals("rsnum") || header.equals("id"))
                      this.SNPCol = colCounter;

                  colCounter++;

              }

              if (this.locationCol < 0 || this.chrCol < 0 || this.pCol < 0 || this.SNPCol < 0)
                  throw new ParserException("Did not find required column headers.", reader.getCurrentLineNumber(), nextLine);


              GWASData gData = new GWASData();

              int rowCounter = 0;
              int indexCounter = 0;
              int addedValuesCounter = 0;

              while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {


                  nextLine = nextLine.trim();
                  rowCounter++;

                  String[] tokens = new String[100];
                  int size = ParsingUtils.splitSpaces(nextLine, tokens);

                  if (tokens.length > 1) {

                      //String chr = ParsingUtils.convertChrString(tokens[chrCol].trim());
                      String chr = genome.getChromosomeAlias(tokens[chrCol].trim());

                      int start;

                      try {
                          start = Integer.parseInt(tokens[locationCol].trim());
                      } catch (NumberFormatException e) {
                          throw new ParserException("Column " + locationCol + " must be a numeric value.", reader.getCurrentLineNumber(), nextLine);
                      }

                      // Check if the p-value is NA
                      if (!tokens[pCol].trim().equals("NA")) {
                          float p;

                          try {
                              p = Float.parseFloat(tokens[pCol].trim());
                              // Transform to -log10
                              p = (float) -log10((double) p);


                          } catch (NumberFormatException e) {
                              throw new ParserException("Column " + pCol + " must be a numeric value.", reader.getCurrentLineNumber(), nextLine);
                          }


                          gData.addLocation(chr, start);
                          gData.addValue(chr, p);
                          indexCounter++;
                          addedValuesCounter++;

                          if (indexCounter == indexSize) {
                              //log.info(nextLine.trim());
                             //log.info("bytes: " + reader.getBytesRead() + " readrows: " + rowCounter + " added locations: " + addedValuesCounter + " tot.locations: " + gData.countTotalLocations());
                              //log.info("bytes: " + reader.getBytesRead() + " readrows: " + rowCounter + " added locations: " + addedValuesCounter);
                              //gData.getFileIndex().add(reader.getBytesRead());
                              gData.getFileIndex().add((int) reader.getPosition());

                              indexCounter = 0;

                          }


                      }

                  }


              }


              return gData;
          }

          catch (
                  ParserException e
                  )

          {
              throw e;
          }

          catch (
                  Exception e
                  )

          {
              if (nextLine != null && reader.getCurrentLineNumber() != 0) {
                  throw new ParserException(e.getMessage(), e, reader.getCurrentLineNumber(), nextLine);
              } else {
                  throw new RuntimeException(e);
              }
          }

          finally

          {
              reader.close();
              fs.close();
          }


      }



}