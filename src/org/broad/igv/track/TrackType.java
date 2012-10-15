/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broad.igv.track;

/**
 * @author jrobinso
 */
public enum TrackType {
    OTHER, COPY_NUMBER, GENE_EXPRESSION, CHIP, DNA_METHYLATION, TILING_ARRAY, PHASTCON,
    ALLELE_SPECIFIC_COPY_NUMBER, LOH, MUTATION, RNAI, POOLED_RNAI, CHIP_CHIP, CNV,
    ALLELE_FREQUENCY, COVERAGE, REPMASK, EXPR, AFFECTIVE, GENE, PLUGIN
}
