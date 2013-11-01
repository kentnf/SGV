#!/usr/bin/perl -w 
use strict;
# ���ȵ����д�����Ȼ���������ӱ�ͷ���
# ����column��Ҫ�������,�ȸ���decription�������Hit_ID����
# Contig_ID	Contig_length	Hit_ID	Hit_length	strand	Contig_start	Contig_end	Hit_start	Hit_end	identity	genus	 description	contig_sequence

if (@ARGV < 2)
{
  print "usage: arrange_col.pl input output\n";
  exit(0);
}

our $input = $ARGV[0]; #�������Ϊindex���ļ�
our $output = $ARGV[1]; #����Ҫ�õ�������ļ���������������ȡ������

my @all_data;
open(IN, "$input");
while (<IN>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	push(@all_data, [@ta[0,13,1,2,3,11,12,6,7,8,9,10,4,5]]);#contig������hit���Ⱥ����� 
}
close(IN);
@all_data = sort { ($a->[6] cmp $b->[6]) || ($a->[4] cmp $b->[4])} @all_data; #����
open(OUT, ">$output");
	print OUT "Contig_ID\tContig_Seq\tContig_Len\tHit_ID\tHit_Len\tGenus\tDescription\tContig_start\tContig_end\tHit_start\tHit_end\tHsp_identity\tE_value\tHsp_strand\n";#���������׼���
	foreach my $each (@all_data){
		print OUT join("\t", @$each)."\n";#���������׼���	
	}
close(OUT);
