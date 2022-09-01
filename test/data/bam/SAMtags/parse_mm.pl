#!/usr/bin/perl -w

# Copyright (c) 2020 Genome Research Ltd.
# Author: James Bonfield <jkb@sanger.ac.uk>
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
#    1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
# 
#    2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
# 
#    3. Neither the names Genome Research Ltd and Wellcome Trust Sanger
# Institute nor the names of its contributors may be used to endorse or promote
# products derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH LTD OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

use strict;

# Complement a sequence
sub comp {
    my ($seq)=@_;
    $seq =~ tr/ACGTUMRWSYKVHDBN/TGCAAKYWSRMBDHVN/;
    return $seq;
}

# Reverse complement a sequence
sub rc {
    my ($seq)=@_;
    $seq = reverse($seq);
    return comp($seq);
}

my $nseq = 0;
while (<>) {
    chomp($_);                                 # trim newline
    next if /^\@/;

    my ($qname,$flag,$rname,$pos,$mapq,$cigar,$rnext,
	$pnext,$tlen,$seq,$qual,@aux) = split("\t",$_);
    my $dir = $flag & 0x10 ? "-" : "+";
    my ($mod_str)="@aux"=~m/M[mM]:Z:(\S+)/;
    my ($mod_prob)="@aux"=~m/M[lL]:B:C,(\S+)/;

    # All counting is from 5' end irrespective of BAM FLAG 0x10
    $seq = rc($seq) if ($dir eq "-");

    my @seq = split("", $seq);                 # array of seq bases
    my @seqt = @seq;                           # orientation as shown in SAM
    my @seqb = split("", comp($seq));          # plus a complemented copy

    $mod_str =~ s/;$//;                        # trim last ; to aid split()
    my @mods = split(";", $mod_str);
    my @probs = split(",", $mod_prob);
    my $pnum = 0;

    print "\n" if $nseq++ > 0;
    foreach (@mods) {
	my ($base, $strand, $types, $pos) = $_ =~ m/([A-Z])([-+])([^,]+),(.*)/;

	my $i = 0; # I^{th} bosition in sequence
	foreach my $delta (split(",", $pos)) {
	    # Skip $delta occurences of $base
	    do {
		$delta-- if ($base eq "N" || $base eq $seq[$i]);
		$i++;
	    } while ($delta >= 0);
	    $i--;

	    # TypePercent combo
	    my $type_perc = "";
	    foreach ($types =~ m/(\d+|.)/g) {
		if (/\d/) {
		    # Avoid ChEBI numbers running into quality values.
		    $type_perc .= "($_)" . int(($probs[$pnum++]+0.5)/256.0*100);
		} else {
		    # Integer qualities represent 1/256th of the total
		    # probability space, with each bin being equal size.
		    # We don't know where in that range our actual
		    # probability was, so we use the mid-point (+0.5 below).
		    $type_perc .= "$_" . int(($probs[$pnum++]+0.5)/256.0*100);
		}
	    }

	    # Add to top or bottom seq
	    if ($strand eq "+") {
		$seqt[$i] .= $type_perc;
	    } else {
		$seqb[$i] .= $type_perc;
	    }
	    $i++;
	}
    }

    for (my $i=0; $i<=$#seq; $i++) {
	print "$seqt[$i]\t$seqb[$i]\n";
    }
}
