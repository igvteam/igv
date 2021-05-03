Mm and Ml auxiliary tags
========================

The purpose of these test files is to test parsing of the Mm and Ml
tags.  These succint Mm and Ml tags are present in the .sam files,
with a more human readable expanded form in the .txt files.
Developers should check whether their implementation is able to
convert between the two forms.

The .sam files are SAM format, but the only fields used for this
test are the reverse-complementation flag (FLAG bit 0x10), the
sequence, and the Mm and Ml tags.

The .txt files uses one line per base, with a blank line separating
sequences.  Each line consists of two tab-separated fields
representing the top (original as-sequenced orientation) and bottom
strand.

Each field in the .txt files is a concatenation of the canonical base
and any modifications with their associated probability, expressed
here as a rounded-down percentage.  The modification and percentage
are joined together with no punctuation, but for ChEBI numeric codes
the code is bracketted to distinguish it from the numeric percentage.
For example "Cm80(76792)73" is canonical base C, modification m (80%)
and modification ChEBI 76792 (73%).


The parse_mm.pl perl script can convert .sam to .txt files, but is for
demonstrative purposes only.  Example usage to check files:

    for i in `echo *.sam|sed 's/\.sam//g'`
    do
        ./parse_mm.pl $i.sam > _
	cmp _ $i.txt
    done
