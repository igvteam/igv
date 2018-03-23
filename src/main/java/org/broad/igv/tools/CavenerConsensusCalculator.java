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

package org.broad.igv.tools;

import java.util.List;

/**
 * Rules taken from Cavener, Nucleic Acids Res. 15, 1353-1361, 1987.
 * Also used by CisRED (note that their literal meaning of rule 2 is slightly off), see http://www.cisred.org/content/methods/help/pfm
 * 1. If the frequency of a single nucleotide at a specific position is greater than 50% and greater than
 * twice the number of the second most frequent nucleotide it is assigned as the consensus nucleotide.
 *
 * 2. If the sum of the frequencies of two nucleotides is greater than 75% (but neither
 * meet the criteria for a single nucleotide assignment) they are assigned as co-consensus nucleotides.
 *
 * 3. If no single nucleotide or pair of nucleotides meet the criteria, we assign an 'n'
 *
 * @author jacob
 * @date 2013-Jun-24
 */
public class CavenerConsensusCalculator extends AbstractConsensusCalculator{


    @Override
    protected char calculateConsensusBase(List<BaseFraction> baseFractions) {

        BaseFraction most = baseFractions.get(0);
        BaseFraction secondMost = baseFractions.get(1);
        float highestFrac = most.fraction;

        if (highestFrac >= 0.5f && highestFrac >= 2.0 * secondMost.fraction) {
            return most.base;
        }else if(most.base == 'n' || secondMost.base == 'n'){
            return 'n';
        } else if(highestFrac + secondMost.fraction >= 0.75f){
            return getDegenerateCode(most.base, secondMost.base);
        }

        return 'n';
    }

}
