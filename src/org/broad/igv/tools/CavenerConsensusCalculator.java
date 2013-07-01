/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
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
