#!/usr/bin/perl -w 
use strict; 
# ����blast_parse_table2.pl��������
# query_name	query_length	hit_name	hit_length	hsp_length	identity	evalue	score	strand	query_start	query_end	hit_start	hit_end
# ����0��1��9��10�Ǳ����
# ��blast�����table�е�query��hit���ǵ�ȫ�����ȼ�������Ȼ�����query�ܳ��ȣ��õ�һ��ratio
# ���е�hit����Ҫ�����һ����identity
# ��query�����ǳ���һ��ratio�ļ�¼������Ȼ���һ��fasta�ļ�����ȡ��¼�в���������Ӧ���У������output1
# ʣ�µ����е����������ouput2
# ���hit���ǲ��������������ratio���ж����õ�contig�ǲ��ǲ���
if (@ARGV < 2)
{
  print "usage: query_filter2.pl inputfile1 inputfile2 output1\n";
  exit(0);
}

our $input1 = $ARGV[0]; #�������Ϊindex���ļ�
our $input2 = $ARGV[1]; #�����blast table
our $output = $ARGV[2]; #����Ҫ�õ�������ļ���������������ȡ������

open(IN1, "$input1");
my %hit_index; 
my %query_exists;#�ⲿ�����ݶ�Ҫ����
while (<IN1>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	defined $hit_index{$ta[0]} or $hit_index{$ta[0]} = 1; 
}
close(IN1);

open(IN2, "$input2");
open(OUT, ">$output");
while(<IN2>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	if(defined $hit_index{$ta[2]}) {#����������hit�������ļ�1�д���
		print OUT join("\t", @ta[0,1,2,3,6,8,9,10,11,12,13])."\n";#��Ӧ��contig��¼�����
		defined $query_exists{$ta[0]} or $query_exists{$ta[0]} = 1;
	}else{#����һ�£�ÿ��contigֻ���һ�Σ�ȥ���ظ���¼��
		if(not defined $query_exists{$ta[0]}){
			print OUT join("\t", @ta[0,1,2,3,6,8,9,10,11,12,13])."\n";#�����
			$query_exists{$ta[0]} = 1;
		}
	}  
}
close(IN2);
close(OUT);



