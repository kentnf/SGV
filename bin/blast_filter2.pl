#!/usr/bin/perl -w 
use strict; 
# usage: ./blast_fileter2.pl inputfile outputfile
# ����blast_parse_table2.pl��������
# query_name	query_length	hit_name	hit_length	hsp_length	identity	evalue	score	strand	query_start	query_end	hit_start	hit_end

# ֻ����ÿ��<query,hit>�� evalue��͵�(����һ��)

if (@ARGV < 2)
{
  print "usage: hit_cov.pl input output1 query_coverage\n";
  exit(0);
}

our $input = $ARGV[0]; #ȥ���ļ���׺����
our $output = $ARGV[1];

open(IN, "$input");
open(OUT, ">$output");

my $last_query="";
my $last_hit="";
while (<IN>) {
	my @cols = split(/\t/, $_); 
	if($cols[0] ne $last_query || $cols[2] ne $last_hit){#�����ǰ������һ����һ��
		print OUT $_;
	}
	$last_query=$cols[0];
	$last_hit=$cols[2];
}
close(IN);
close(OUT);


