#!/usr/bin/perl -w 
use strict; 
# usage: ./hit_cov1.pl inputfile outputfile identity query_cov > output2
# ����blast_parse_table2.pl��������
# query_name	query_length	hit_name	hit_length	hsp_length	identity	evalue	score	strand	query_start	query_end	hit_start	hit_end
# ����2��3��4��5��11��12��(0��ʼ)�Ǳ����
# ��blast�����table�е�ÿ��hit��query���ǵ�hspȫ���������ĳ��ȣ�Ȼ����Ը�hit�ܳ��ȣ��õ�һ��ratio(���)
# �����Եõ����������������hit��query
# ����ÿ��hit��ֻ���������ƶȸ���identity��query������query_cov��hsp�����Կ���
if (@ARGV < 2)
{
  print "usage: hit_cov1.pl input output1 identity query_cov > output2\n";
  exit(0);
}

our $input = $ARGV[0]; #ȥ���ļ���׺����
our $output = $ARGV[1];
our $identity = $ARGV[2];
our $query_cov = $ARGV[3];

open(IN, "$input");
open(OUT, ">$output");

my %blk; 
my %hit_len; #hit�ĳ���
while (<IN>) {
	chomp; 
	my @ta = split(/\t/, $_);
	my $query_covered= 1.0*$ta[4]/$ta[1]; #hsp length����query length	
	if ($ta[5]>=$identity && $query_covered>=$query_cov){#ֻ�������identityҪ���hsp
		push(@{$blk{$ta[2]}}, [@ta[11,12]]); #����hit��(hit_start,hit_end)֮���ӳ��
		$blk{$ta[2]}[-1][2]{$ta[0]}=1; 
		defined $hit_len{$ta[2]} or $hit_len{$ta[2]} = $ta[3];#����hit��hit_length֮���ӳ�� 
	}
}
close(IN);

#print OUT join("\t", qw/hit_name hit_len total_cov cov_rate% query_names/)."\n";   #�����ļ�����ĸ�������
print join("\t", qw/hit_name block_len block_start block_end covered_len query_names/)."\n"; #���׼�������ĸ�������

for my $tk (sort keys %blk) {#�ȸ���hit����
	my @o; #�洢û��overlap��hit�ϵ�block
	for my $ar (sort {$a->[0]<=>$b->[0] || $a->[1]<=>$b->[1] || $a->[2] cmp $b->[2];} @{$blk{$tk}}) {
		if (scalar(@o) == 0) {
			push(@o, [@$ar]); 
		}elsif ($o[-1][1] >= $ar->[0]-1) {#��һ��hit���ص�������ϲ�
			for my $qname (keys %{$ar->[2]}) {
				$o[-1][2]{$qname}=1; 
			}
			$o[-1][1] < $ar->[1] and $o[-1][1] = $ar->[1]; 
		}else{
			push(@o, [@$ar]); 
		}
	}

	my $total_cov = 0;#hit�����ǵ��ܳ��� 
	my %aa; 
	for my $ar (@o) {
		my @query_names = sort keys %{$ar->[2]}; 
		print join("\t", $tk, $hit_len{$tk}, $ar->[0], $ar->[1], $ar->[1]-$ar->[0]+1, join(",", @query_names))."\n"; #���׼����������
		#hit���ƣ�hit���ȣ�hitһ��block����㣬hitһ��block���յ㣬hitһ��block�ĳ���
		$total_cov += ($ar->[1]-$ar->[0]+1); #һ��hit�����з��ص�block�ĳ���֮��
		@aa{@query_names} = (1) x scalar(@query_names); #������contigs�����ƴ���hash��
	}
	print OUT join("\t", $tk, $hit_len{$tk}, $total_cov, $total_cov/$hit_len{$tk}*1.0, join(",", sort keys %aa), scalar(keys %aa))."\n"; #�����ļ��������
	#hit���ƣ�hit���ȣ�hit������bp����hit�����ǵİٷֱȣ�����contigs���ƣ�contigs����
}
close(OUT);


