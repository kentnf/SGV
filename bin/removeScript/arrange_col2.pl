#!/usr/bin/perl -w 
use strict;
# ���ȵ����д�����Ȼ���������ӱ�ͷ���
# ����column��Ҫ�������,�ȸ���decription�������Hit_ID����
# Contig_ID	Contig_length	Hit_ID	Hit_length	strand	Contig_start	Contig_end	Hit_start	Hit_end	identity	genus	 description	contig_sequence

if (@ARGV < 1)
{
  print "usage: arrange_col2.pl input > output\n";
  exit(0);
}

our $input = $ARGV[0]; #�����ļ�

my @all_data;
open(IN, "$input");
while (<IN>) {
	chomp; 
	my @ta = split(/\t/, $_);
    my $coverage= 1.0*$ta[7]/$ta[1];	
	push(@all_data, [@ta[0,1,7],$coverage,@ta[4,5,6,8,9]]);#ѡ����Ҫ���У�����������˳�� 
}
close(IN);
@all_data = sort { ($a->[7] cmp $b->[7]) || ($b->[2] <=> $a->[2])} @all_data; #����Genus��hit_covered(bp)����

foreach my $each (@all_data){
	print join("\t", @$each)."\n";#���������׼���	
}
