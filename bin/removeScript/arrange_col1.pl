#!/usr/bin/perl -w 
use strict;
# ���ȵ����д�����Ȼ���������ӱ�ͷ���
# ����column��Ҫ�������,�ȸ���decription�������Hit_ID����
# Contig_ID	Contig_length	Hit_ID	Hit_length	strand	Contig_start	Contig_end	Hit_start	Hit_end	identity	genus	 description	contig_sequence

if (@ARGV < 1)
{
  print "usage: arrange_col1.pl input > output\n";
  exit(0);
}

our $input = $ARGV[0]; #�����ļ�

my @all_data;
open(IN, "$input");
while (<IN>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	push(@all_data, [@ta[0,19,1,2,3,17,18,9,10,11,12,13,6,8]]);#ѡ����Ҫ���У�����������˳�� 
}
close(IN);
@all_data = sort { ($a->[5] cmp $b->[5]) || ($a->[3] cmp $b->[3])} @all_data; #����Genus��hit����

print "Contig_ID\tContig_Seq\tContig_Len\tHit_ID\tHit_Len\tGenus\tDescription\tContig_start\tContig_end\tHit_start\tHit_end\tHsp_identity\tE_value\tHsp_strand\n";#���������׼���
foreach my $each (@all_data){
	print join("\t", @$each)."\n";#���������׼���	
}
