#!/usr/bin/perl -w 
use strict; 
# 来自blast_parse_table2.pl程序的输出
# query_name	query_length	hit_name	hit_length	hsp_length	identity	evalue	score	strand	query_start	query_end	hit_start	hit_end
# 其中0、1、9、10是必须的
# 将blast结果的table中的query被hit覆盖的全部长度加起来，然后除以query总长度，得到一个ratio
# 所有的hit必须要求大于一定的identity
# 把query被覆盖超过一定ratio的记录下来，然后从一个fasta文件中提取记录中不包括的相应序列，输出到output1
# 剩下的序列的名称输出到ouput2
# 如果hit都是病毒，可以用这个ratio来判断所得到contig是不是病毒
if (@ARGV < 2)
{
  print "usage: query_filter2.pl inputfile1 inputfile2 output1\n";
  exit(0);
}

our $input1 = $ARGV[0]; #这个是作为index的文件
our $input2 = $ARGV[1]; #这个是blast table
our $output = $ARGV[2]; #这是要得到的输出文件，包括了所有提取的序列

open(IN1, "$input1");
my %hit_index; 
my %query_exists;#这部分数据都要处理
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
	if(defined $hit_index{$ta[2]}) {#如果包括这个hit在输入文件1中存在
		print OUT join("\t", @ta[0,1,2,3,6,8,9,10,11,12,13])."\n";#对应的contig记录就输出
		defined $query_exists{$ta[0]} or $query_exists{$ta[0]} = 1;
	}else{#否则看一下，每个contig只输出一次（去除重复记录）
		if(not defined $query_exists{$ta[0]}){
			print OUT join("\t", @ta[0,1,2,3,6,8,9,10,11,12,13])."\n";#就输出
			$query_exists{$ta[0]} = 1;
		}
	}  
}
close(IN2);
close(OUT);



