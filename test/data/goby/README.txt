April 27 2011:
Goby 1.9.5 no longer supports queryLengths stored in the aligment header.
This usage was deprecated in Goby 1.7. Versions of Goby 1.7-1.9.4 had a concat mode that transfered
information from the header to the alignment entries transparently.
Use  concatenate-alignments of Goby 1.7-1.94 when you need to migrate a 1.6- alignment to work with Goby 1.9.5+.

The alignment files in this directory were converted to work with Goby 1.9.5 with the concatenate-alignment
mode of Goby 1.9.1, as follows:

goby_1.9.1/goby 1g concatenate-alignments test/data/goby/DLTTEJH-Bullard-HBR-SRR037439 \
               -o test/data/goby/DLTTEJH-Bullard-HBR-SRR037439-1.9.5
