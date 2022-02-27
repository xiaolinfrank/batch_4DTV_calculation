#!/usr/bin/perl
use strict;

##author: sun ming'an, sunma@genomics.org.cn
##modifier: fanwei, fanw@genomics.org.cn
##correction: LiJun, junli@genomics.org.cn
##Date: 2008-9-24

## Correction and batch function: Huang Peilin, huangplin@mail2.sysu.edu.cn / peilinfrank@foxmail.com
## Date: 2022-02-27

##4dtv (transversion rate on 4-fold degenerated sites) are calculated with HKY substitution models
##Reference: M. Hasegawa, H. Kishino, and T. Yano, J. Mol. Evol. 22 (2), 160 (1985)

die "perl $0 AXTfile > outfile\n" unless( @ARGV == 1);


my %codons=(
'CTT'=>'L', 'CTC'=>'L', 'CTA'=>'L', 'CTG'=>'L',
'GTT'=>'V', 'GTC'=>'V', 'GTA'=>'V', 'GTG'=>'V',
'TCT'=>'S', 'TCC'=>'S', 'TCA'=>'S', 'TCG'=>'S',
'CCT'=>'P', 'CCC'=>'P', 'CCA'=>'P', 'CCG'=>'P',
'ACT'=>'T', 'ACC'=>'T', 'ACA'=>'T', 'ACG'=>'T',
'GCT'=>'A', 'GCC'=>'A', 'GCA'=>'A', 'GCG'=>'A',
'CGT'=>'R', 'CGC'=>'R', 'CGA'=>'R', 'CGG'=>'R',
'GGT'=>'G', 'GGC'=>'G', 'GGA'=>'G', 'GGG'=>'G'
);

my %transversion = (
        "A" => "TC",
        "C" => "AG",
        "G" => "TC",
        "T" => "AG",
);

my $axtFile = shift;

open(AXT,"$axtFile")||die"Cannot open $axtFile\n";
#$/ = "\n\n";
#my @seqs = <AXT>;
#$/ ="\n";
#close AXT;

print "tag\t4dtv_corrected\t4dtv_raw\tcondon_4d\tcodon_4dt\n";
#foreach my $line ( @seqs ){
#        chomp $line;
while(<AXT>){
#        if( $line =~ /^(\S+)\n(\S+)\n(\S+)$/ ){
	if( ~/(\S+)\n(\S+)\n(\S+)$/ ){
                my $tag = $_;
		chomp($tag);
                my $seq1 = <AXT>;
                my $seq2 = <AXT>;
                my ($corrected_4dtv, $raw_4dtv, $condon_4d, $codon_4dt) = &calculate_4dtv($seq1, $seq2);
                print "$tag\t$corrected_4dtv\t$raw_4dtv\t$condon_4d\t$codon_4dt\n";
        }
}

close AXT;

sub calculate_4dtv {
        my($str1, $str2) = @_;

        my ($condon_4d, $codon_4dt) = (0,0);
		my ($V,$a,$b,$d,$tt,$t) = (0,0,0,0,0,0);
		my %fre=();
        for( my $i = 0; $i < length($str1); $i += 3){
                my $codon1 = substr($str1, $i, 3);
                my $codon2 = substr($str2, $i, 3);
                my $base1= uc(substr($str1, $i+2, 1));
                my $base2= uc(substr($str2, $i+2, 1));

                if( exists $codons{$codon1} && exists $codons{$codon2} && $codons{$codon1} eq $codons{$codon2} ){
					$fre{$base1}++;
					$fre{$base2}++;
                    $condon_4d++;
                    $codon_4dt++ if(is_transversion($codon1,$codon2));
                }
        }

		if($condon_4d > 0){
			$V=$codon_4dt / $condon_4d; ##this is raw 4dtv value
			##correction the raw 4dtv values by HKY substitution model
			$fre{"Y"}=$fre{"T"}+$fre{"C"};
			$fre{"R"}=$fre{"A"}+$fre{"G"};
			foreach (keys %fre){
				$fre{$_}=0.5*$fre{$_}/$condon_4d;
			}
      $t=2*($fre{T}*$fre{C}*$fre{R}+$fre{A}*$fre{G}*$fre{Y});
      if($t!=0){
        $tt=$V*($fre{T}*$fre{C}*$fre{R}/$fre{Y}+$fre{A}*$fre{G}*$fre{Y}/$fre{R})/$t;
      }else{
        $tt=1;
      }


			if($fre{Y}!=0 && $fre{R}!=0 && $fre{A}!=0 && $fre{C}!=0 && $fre{G}!=0 && $fre{T}!=0 && $tt!=1){
				$a=-1*log(1-$V*($fre{T}*$fre{C}*$fre{R}/$fre{Y}+$fre{A}*$fre{G}*$fre{Y}/$fre{R})/(2*($fre{T}*$fre{C}*$fre{R}+$fre{A}*$fre{G}*$fre{Y})));
				if (1-$V/(2*$fre{Y}*$fre{R}) > 0) {
					$b=-1*log(1-$V/(2*$fre{Y}*$fre{R}));
					$d=2*$a*($fre{T}*$fre{C}/$fre{Y}+$fre{A}*$fre{G}/$fre{R})-2*$b*($fre{T}*$fre{C}*$fre{R}/$fre{Y}+$fre{A}*$fre{G}*$fre{Y}/$fre{R}-$fre{Y}*$fre{R});
				}else{
					$d = "NA";
				}
			}else{
				$d = "NA";
			}


		}else{
			$V="NA";
			$d="NA";
		}

        return ($d,$V,$condon_4d, $codon_4dt);

}


sub is_transversion{
        my ($codon1,$codon2) = @_;
        my $is_transversion = 0;
        my $base1 = substr($codon1,2,1);
        my $base2 = substr($codon2,2,1);
        $is_transversion = 1 if (exists $transversion{$base1} && $transversion{$base1} =~ /$base2/);
        return $is_transversion;
}
