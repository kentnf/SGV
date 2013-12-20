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
my %hit_index; #�����ڵ�һ���ļ��г��ֵ�hit����Ϊreference

while (<IN1>) {
	chomp; 
	my @cols = split(/\t/, $_); 
	defined $hit_index{$cols[0]} or $hit_index{$cols[0]} = 1; #�����г��ֵ�hit��hash����
}
close(IN1);

open(IN2, "$input2");
open(OUT, ">$output");
my $last_query="";
my $current_query;

my $high_evalue;#����ÿ��query��ߵ�identity����ʵ���ǵ�һ�γ��ֵ�(��Ϊ������
my $current_evalue;

my $high_identity;#����ÿ��query��ߵ�identity
my $current_identity;

my $high_record;#����ǰһ��query�����hsp���м�¼
my $if_has_ref=0;#��ʾǰһ��query��ǰqueryû��reference

while (<IN2>) {
	my @cols = split(/\t/, $_);
	$current_query=$cols[0];
	$current_identity=$cols[5];
	$current_evalue=$cols[6];	
	if($current_query ne $last_query){#ÿ�γ���һ���µ�query
		if($if_has_ref==0 && $last_query ne ""){#�ȿ�ǰһ��query������hsp����ref�������û��ref
			print OUT $high_record;#�Ͱ�ǰһ��query������hsp�����hsp���
		}
		#�������ϸ�query����ʼ������ǰquery�ĵ�һ��hsp
		$if_has_ref=0;   #��ʾ��ǰquery��û���ҵ�ref
		if (defined($hit_index{$cols[2]})){#�����һ��hsp����ref
			$high_identity=$cols[5];#���������ߵ�identity
			$high_evalue=$cols[6];	#���������ߵ�evalue
			print OUT $_;           #���������¼
			$if_has_ref=1;	        #��ʾ��ǰquery�ҵ���ref		
		}else{
			$high_record=$_; #�����һ��û��ref������Ϊ��ߵģ�������ε�����hsp��û��ref�������´�������			
		}
	}else{#���ڵ�ǰquery������hsp
		if($if_has_ref==1){#�����ǰquery�Ѿ��ҵ���ref
			if (defined($hit_index{$cols[2]}) && $current_evalue>=$high_evalue && $current_identity>=$high_identity)
			{
				print OUT $_;#ֻ�����refͬʱ����õ�hspͬ��identity��evalue�Ľ��(ע����>=,������п��ܸ���)			
			}
		}else{#�����ǰquery��û�ҵ�ref
			if (defined($hit_index{$cols[2]})){#һ���ҵ�ref
				$high_identity=$cols[5];#���������ߵ�
				$high_evalue=$cols[6];		
				print OUT $_;
				$if_has_ref=1;	#��ʾ��ǰ��query�ҵ���ref
			}
		}
	
	}
	$last_query=$cols[0];
}
close(IN2);
close(OUT);



