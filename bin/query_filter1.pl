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
  print "usage: query_filter1.pl inputfile1 inputfile2 output1 identity min_ratio > output2\n";
  exit(0);
}

our $input1 = $ARGV[0]; #
our $input2 = $ARGV[1]; #
our $output = $ARGV[2]; #����Ҫ�õ�������ļ���������������ȡ������
our $identity = $ARGV[3];#��Ҫ�ۻ���hsp��identity��Сֵ
our $min_ratio = $ARGV[4];#һ��query�ۼƱ�covered����С����

open(IN1, "$input1");
my %blk; 
my %query_len; 
my %query_filtered;
while (<IN1>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	if ($ta[5]>=$identity){#ֻ�������identityҪ��ļ�¼
		push(@{$blk{$ta[0]}}, [@ta[9,10]]); #����query��(query_start,query_end)֮���ӳ��
		defined $query_len{$ta[0]} or $query_len{$ta[0]} = $ta[1];#����query��query_length֮���ӳ�� 
	}
}
close(IN1);

#print join("\t", qw/query_name block_len block_start block_end covered_len/)."\n"; #���׼�������ĸ�������
#print OUT join("\t", qw/query_name query_len total_cov cov_rate%/)."\n";   #�����ļ�����ĸ�������

for my $tk (sort keys %blk) {#�ȸ���query����Ȼ����ȡ������Ҫfiltered����query����
	my @o; #�洢û��overlap��query�ϵ�block
	for my $ar (sort { $a->[0]<=>$b->[0] || $a->[1]<=>$b->[1];} @{$blk{$tk}}) {
		if (scalar(@o) == 0) {
			push(@o, [@$ar]); 
		}elsif ($o[-1][1] >= $ar->[0]-1) {#��һ��query���ص�������ϲ�
			$o[-1][1] < $ar->[1] and $o[-1][1] = $ar->[1]; 
		}else{
			push(@o, [@$ar]); 
		}
	}

	my $total_cov = 0; 
	for my $ar (@o) {
		#print join("\t", $tk, $query_len{$tk}, $ar->[0], $ar->[1], $ar->[1]-$ar->[0]+1)."\n"; #���׼�������ĸ���
		#query���ƣ�query���ȣ�queryһ��block����㣬queryһ��block���յ㣬queryһ��block�ĳ���
		$total_cov += ($ar->[1]-$ar->[0]+1); #һ��query�����з��ص�block�ĳ���֮��
	}
	my $ratio=int($total_cov/$query_len{$tk}*10000+0.5)/100;
	if($ratio>=$min_ratio){#���з���������query,����Ҫ������
		defined $query_filtered{$tk} or $query_filtered{$tk} = 1;
	}
	#�����ļ�����ĸ���
	#query���ƣ�query���ȣ�query�����ǵ��ܳ��ȣ�query�����ǵİٷֱ�
}

#��input2�аѲ�������%query_filtered�е�������ȡ�����������OUT
open(IN2, $input2);
open(OUT, ">$output");
my $flag = "off";
while(<IN2>) {
	if($_ =~ m/^>/) {
		my $head = $_;
		chomp($head);
		$head=~s/>//;

		if(defined $query_filtered{$head}) {#����������name
			print $head."\t";#���������׼���
			$flag = "on";#ͬʱ�ı��־����ʾ�����������Ҫ������OUT1���
		}
		else {#������������name
			print OUT $_;#�������OUT
			$flag = "off";#ͬʱ�ı��־����ʾ�����������Ҫ������OUT2���
		}
	}
	else {
		if($flag eq "on") {#��ʾΪ"on"
			print $_;#�����������Ҫ�������׼������
		}
		else {
			print OUT $_;#���򣬺����������Ҫ������OUT���
		}
	}
}
close(IN2);
close(OUT);



