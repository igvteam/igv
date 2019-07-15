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

package org.broad.igv.feature.tribble;

import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.readers.AsciiLineReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.feature.Mutation;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.collections.MultiMap;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;

/**
 * Codec for .mut and .maf mutation files
 *
 * @author jrobinso
 * Date: 10/24/12
 * Time: 12:05 PM
 */
public class MUTCodec extends AsciiFeatureCodec<Mutation> {

    private static Logger log = Logger.getLogger(MUTCodec.class);

    private String path;  // for error messages
    private boolean isMAF;
    private String[] headers;
    private String[] samples;
    private Genome genome;
    private int chrColumn;
    private int startColumn;
    private int endColumn;
    private int sampleColumn = -1;
    private int typeColumn = -1;
    private int refAlleleColumn = -1;
    private int tumorAllele1Column = -1;
    private int tumorAllele2Column = -1;
    private int errorCount = 0;


    public MUTCodec(String path, Genome genome) {
        super(Mutation.class);
        this.path = path;
        this.genome = genome;
        try {
            LineIterator reader = new LineIteratorImpl(new AsciiLineReader(ParsingUtils.openInputStream(path)));
            readActualHeader(reader);
        } catch (IOException e) {
            log.error(e.getMessage(), e);
        }
    }

    public Void readActualHeader(LineIterator reader) {
        String nextLine = null;

        while (reader.hasNext()) {
            nextLine = reader.peek();
            if (nextLine.startsWith("#")) {
                if (nextLine.startsWith("#samples")) {
                    String[] tokens = Globals.whitespacePattern.split(nextLine, 2);
                    if (tokens.length < 2) {
                        log.error("Error parsing sample header in mutation file: " + path);
                    } else {
                        samples = Globals.commaPattern.split(tokens[1]);
                        for (int i = 0; i < samples.length; i++) {
                            samples[i] = samples[i].trim();
                        }
                    }
                }
                reader.next();
                continue;
            }

            String[] tokens = Globals.tabPattern.split(nextLine);
            if (tokens.length >= 5) {
                reader.next();
                headers = tokens;
                isMAF = headers[0].trim().equalsIgnoreCase("Hugo_Symbol") || isMafLiteFile(headers);
                setColumns(isMAF, tokens);
                return null;
            } else {
                throw new RuntimeException(String.format("Not enough columns in header line found in %s: %s", path, nextLine));
            }


        }
        throw new RuntimeException("Unexpected end-of-file (no header line): " + path);
    }

    public String[] getSamples() {
        return samples;
    }

    public int getChrColumn() {
        return chrColumn;
    }

    public int getStartColumn() {
        return startColumn;
    }

    @Override
    public Mutation decode(String line) {

        if (line.startsWith("#") || line.startsWith("Hugo_Symbol")) {
            return null;
        } else {

            String[] tokens = Globals.tabPattern.split(line);

            String chr = genome == null ? tokens[chrColumn].trim() : genome.getCanonicalChrName(tokens[chrColumn].trim());

            int start;
            try {
                start = Integer.parseInt(tokens[startColumn].trim());
            } catch (NumberFormatException e) {
                errorCount++;
                if (errorCount > 100) {
                    throw new DataLoadException("Column " + (startColumn + 1) + " must be a numeric value.", path);
                } else {
                    log.info("Error parsing line: " + line);
                    return null;
                }
            }

            int end;
            try {
                end = Integer.parseInt(tokens[endColumn].trim());
            } catch (NumberFormatException e) {
                errorCount++;
                if (errorCount > 100) {
                    throw new DataLoadException("Column " + (endColumn + 1) + " must be a numeric value.", path);
                } else {
                    log.info("Error parsing line: " + line);
                    return null;
                }
            }


            // MAF files use the 1-based inclusive convention for coordinates.  The convention is not
            // specified for MUT files, and it appears both conventions have been used.  We can detect
            // the convention used for single base mutations by testing start == end.
            if (isMAF || (start == end)) {
                start--;
            }

            String sampleId = "Unknown";
            if(sampleColumn >= 0) {
                sampleId = tokens[sampleColumn].trim();
            }
            String type = "Unknown";
            if(typeColumn >= 0) {
                type =
                        tokens[typeColumn].trim();
            }

            MultiMap<String, String> attributes = new MultiMap();
            int n = Math.min(headers.length, tokens.length);
            for (int i = 0; i < n; i++) {
                String key = headers[i];
                String value = tokens[i];
                if (value.length() > 0) {
                    attributes.put(key, value);
                }
            }


            Mutation mut = new Mutation(sampleId, chr, start, end, type);
            mut.setAttributes(attributes);

            if (refAlleleColumn > 0) {
                mut.setRefAllele(tokens[refAlleleColumn].trim());
            }
            if (tumorAllele1Column > 0) {
                mut.setAltAllele1(tokens[tumorAllele1Column].trim());
            }
            if (tumorAllele2Column > 0) {
                mut.setAltAllele2(tokens[tumorAllele2Column].trim());
            }

            return mut;
        }
    }

/*
TCGA
-----
Hugo_Symbol
Entrez_Gene_Id
Center
NCBI_Build
Chromosome
Start_position
End_position
Strand
Variant_Classification
Variant_Type
Reference_Allele
Tumor_Seq_Allele1
Tumor_Seq_Allele2
dbSNP_RS
dbSNP_Val_Status
Tumor_Sample_Barcode
Matched_Norm_Sample_Barcode
Match_Norm_Seq_Allele1
Match_Norm_Seq_Allele2
Tumor_Validation_Allele1
Tumor_Validation_Alliele2
Match_Norm_Validation_Allele1
Match_Norm_Validation_Allele2
Verification_Status
Validation_Status
Mutation_Status	Sequencing_Phase

BROAD
-----
Hugo_Symbol
Entrez_Gene_Id
Chromosome
Start_position
End_position
Variant_Classification
Variant_Type
Reference_Allele
Tumor_Seq_Allele1
Tumor_Seq_Allele2
dbSNP_RS
dbSNP_Val_Status
Tumor_Sample_Barcode
Matched_Norm_Sample_Barcode
n_ref_count
t_ref_count
n_alt_count
t_alt_count

MAFLITE
-------
build
chr
start
end
ref_allele
alt_allele


tumor_barcode
normal_barcode
NCBI_Build
Strand
Center
source
status
phase
sequencer
Tumor_Validation_Allele1
Tumor_Validation_Allele2
Match_Norm_Validation_Allele1
Match_Norm_Validation_Allele2
Verification_Status
Validation_Status
Validation_Method
Score
BAM_file
Match_Norm_Seq_Allele1
Match_Norm_Seq_Allele2

 */

    private void setColumns(boolean isMAF, String[] tokens) {

        this.isMAF = isMAF;
        if (isMAF) {
            for (int i = 0; i < tokens.length; i++) {
                switch (tokens[i].toLowerCase()) {
                    case "chromosome":
                    case "chr":
                        chrColumn = i;
                        break;
                    case "start_position":
                    case "start":
                        startColumn = i;
                        break;
                    case "end_position":
                    case "end":
                        endColumn = i;
                        break;
                    case "tumor_sample_barcode":
                    case "tumor_barcode":
                        sampleColumn = i;          // Tumor_Sample_Barcode
                        break;
                    case "variant_classification":
                        typeColumn = i;
                        break;
                    case "reference_allele":
                    case "ref_allele":
                        refAlleleColumn = i;
                        break;
                    case "tumor_seq_allele1":
                    case "alt_allele":
                        tumorAllele1Column = i;
                        break;
                    case "tumor_seq_allele2":
                        tumorAllele2Column = i;
                        break;
                }
            }
        } else {
            chrColumn = 0;
            startColumn = 1;
            endColumn = 2;
            sampleColumn = 3;
            typeColumn = 4;
            refAlleleColumn = -1;
            tumorAllele1Column = -1;
            tumorAllele2Column = -1;
        }
    }

    /**
     * @param locator
     * @return
     * @throws IOException
     */
    public static boolean isMutationAnnotationFile(ResourceLocator locator) {

        String ext;
        {
            //Only want pathLC when checking extension, so we limit its scope
            String typeStringLC = locator.getTypeString().toLowerCase();
            if (typeStringLC.endsWith(".maf.annotated")) {
                // TCGA extension
                return true;
            }

            ext = ParsingUtils.getIGVExtension(typeStringLC);
        }

        if (ext.equals("mut")) {
            return true;
        } else  {

            BufferedReader reader = null;
            String path = locator.getPath();

            if(path == null) {
                return false;
            }

            try {
                reader = ParsingUtils.openBufferedReader(path);
                if (reader == null) {
                    return false;
                }

                String nextLine;
                while ((nextLine = reader.readLine()) != null && nextLine.startsWith("#")) {
                    if (nextLine.startsWith("#")) {
                        continue;
                    }
                }

                String[] tokens = Globals.tabPattern.split(nextLine);
                if (tokens.length > 5 && tokens[0].equalsIgnoreCase("Hugo_Symbol")) {
                    return true;
                } else if (isMafLiteFile(tokens)) {
                    return true;
                } else {
                    return false;
                }
            } catch (IOException e) {
                return false;
            } finally {
                if (reader != null) {
                    try {
                        reader.close();
                    } catch (IOException e) {
                        log.error("Error closing: " + path, e);
                    }
                }
            }
        }
    }


    public static boolean isMafLiteFile(String[] tokens) {
        HashSet<String> tokenSet = new HashSet(Arrays.asList(tokens));
        String[] requiredFields = {"build", "chr", "start", "end", "ref_allele", "alt_allele"};
        for (String f : requiredFields) {
            if (!tokenSet.contains(f)) {
                return false;
            }
        }
        return true;
    }

    @Override
    public boolean canDecode(String path) {
       return isMutationAnnotationFile(new ResourceLocator(path));
    }
}
