table BisulfiteSeq
"BED6 + 2 scores for bisulfite-seq CpG data"
    (
    string chrom;      "Reference chromosome or scaffold"
    uint   chromStart; "Start position in chromosome"
    uint   chromEnd;   "End position in chromosome"
    string name;       "Empty"
    uint   score;      "Score from 0-1000. PercentMeth*10
    char[1] strand;    "+ or - or . for unknown"
    float percentMeth;       "PercentMeth 0-100 (percent of C or T reads with C value)"
    uint numCorT;       "Number of C or T reads"
    )