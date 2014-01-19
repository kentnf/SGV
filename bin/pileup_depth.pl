#!/usr/bin/perl -w 
use strict; 
# usage: pileup_depth.pl inputfile outputfile
# ��pileup����ȡÿ��chr��ƽ��depth��Ϣ


if (@ARGV < 1)
{
  print "usage: pileup_depth.pl inputfile outputfile\n";
  exit(0);
}

our $input = $ARGV[0]; #ȥ���ļ���׺����
our $output = $ARGV[1];

my $last_chr = ''; 
my %depth;#�����ϣ��ӳ��chr������ƽ��depth
my $mean_depth;
my $total_depth = 0; 
my $total_positions = 1; 
open(IN, "$input");
while (<IN>) {
	my @cols = split;
	if ($cols[0] ne $last_chr) { #�����ǰchr��ǰһ����ͬ
                if ($last_chr){      #���Ҳ��ǵ�һ��chr
			#$depth{$last_chr} = 1.0*$total_depth/$total_positions;#�洢chr��Ӧ��ƽ��depth
			$mean_depth = 1.0*$total_depth/$total_positions;
			push(@{$depth{$last_chr}}, ($mean_depth,$total_positions)); 
		}			
		$total_depth = $cols[3];#���³�ʼ��
		$total_positions = 1;   #���³�ʼ��
	}
	else{
		$total_depth = $total_depth + $cols[3];#�ۼӵ�4�������Ϣ
		$total_positions = $total_positions + 1;#�ۼ�λ����Ϣ
	}
        $last_chr = $cols[0];	
}
#$depth{$last_chr} = 1.0*$total_depth/$total_positions;
$mean_depth = 1.0*$total_depth/$total_positions;
push(@{$depth{$last_chr}}, ($mean_depth,$total_positions));  
close(IN);

open(OUT, ">$output");
foreach my $chr (keys(%depth))
{
	print OUT $chr."\t".$depth{$chr}[0]."\t".$depth{$chr}[1]."\n";#���hit���ƣ�ƽ����Ⱥ�hit�����ǵ�λ������
}
close(OUT);


