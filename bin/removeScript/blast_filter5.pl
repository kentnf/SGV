#!/usr/bin/perl -w 
use strict; 

# ����blast_parse_table2.pl��������
# query_name	query_length	hit_name	hit_length	hsp_length	identity	evalue	score	strand	query_start	query_end	hit_start	hit_end
# query��hit����Ҫ��hsp����һ������
# ��ͬ��<query,hit>�У�ֻ����evalue��͵�(����һ��)
# һ��query��Ӧ��hits�� ֻ���������evalue�Ǹ�hsp��identity����drop_off(5%)���ڵ�
if (@ARGV < 2)
{
  print "usage: blast_filter5.pl input output coverage identity_dropoff\n";
  exit(0);
}

our $input = $ARGV[0]; #ȥ���ļ���׺����
our $output = $ARGV[1];
our $coverage=$ARGV[2];
our $identity_dropoff = $ARGV[3];

open(IN, "$input");
open(OUT, ">$output");

my $last_query="";
my $last_hit="";
my $current_query;
my $current_hit;

my $high_identity;#����ÿ��query��ߵ�identity
my $current_identity;
while (<IN>) {
	my @cols = split(/\t/, $_);
	if(3.0*$cols[4]/$cols[1] >= $coverage || 3.0*$cols[4]/$cols[3] >= $coverage){#query��hit����Ҫ��hsp����һ������
		$current_query=$cols[0];
		$current_hit=$cols[2];	
		$current_identity=$cols[5];
		if($current_query ne $last_query){#�������һ���µ�query
			$high_identity=$cols[5];#��¼��ߵ�identity
			print OUT $_;#query�ĵ�һ��hit�϶�Ҫ����Ϊ��evalue��͵�
		}else{
			if ($current_hit ne $last_hit && ($current_identity >= $high_identity-$identity_dropoff))#��ͬ��<query,hit>�У�ֻ����evalue��͵�(����һ��)
			{
				print OUT $_;	
			}	
		}
		$last_query=$cols[0];
		$last_hit=$cols[2];
	}
}
close(IN);
close(OUT);


