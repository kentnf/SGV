#!/usr/bin/perl -w 
#
use strict; 

-t and !@ARGV and die "perl $0 2.sam\n"; 

our $max_distance = 2; #Ҫ��ѡ�е�����hit�ı༭���붼���ܴ���$max_distance
my @lls; 
my $pre_query_name; 
my $nmv; #ÿ��query����һ����С�༭���룬������������hits�ű���
$pre_query_name = ''; 
$nmv = -1; 
my ($query_col, $opt_col) = (0, 11); #sam�ļ�һ���еĵ�1�к͵�12�У��ֱ��Ӧquery��option�Ŀ�ʼ��
while (<>) {
	chomp; s/\s+$//; 
	my @ta = split(/\t/, $_);#ta�������sam�ļ�һ���е�������
	my $query_name = $ta[$query_col];#���뵱ǰquery������ 

	if ($pre_query_name eq '') {#���ǰһ��queryΪ�գ���ֵ
		$pre_query_name = $query_name; 
	}elsif ($query_name ne $pre_query_name) {#�����ǰquery��ǰһ����ͬ
		$nmv < $max_distance and print STDOUT join("\n", @lls)."\n";#nmv����С��$max_distance���ű�����ע��bwa��ʱ���������Ԥ��$max_distance��hits
		$nmv=-1; 
		$pre_query_name = $query_name; 
	}

	for my $each_opt (@ta[$opt_col .. $#ta]) {#ɨ��11�е����һ�е�����
		if ($each_opt =~ /^NM:i:(\d+)$/) {#�ӵ�ǰ�����ҵ�'NM:i:X'ģʽ
			my $distance = $1; #ȡ��query��reference����֮���distance
			if ($nmv == -1 or $distance < $nmv) {#����ҵ���С�ľ���
				$nmv = $distance; #�͸�����С����
				@lls = ($_); 
			}elsif ($distance == $nmv) {
				push(@lls, $_); 
			}else{
				; 
			}
			last; #�ҵ��������NM:i:<N>ģʽ������������Ҫ��ÿ�е�NM:i:<N>ֻ����һ��
		}
	}
}

if ($nmv > -1) {
	$nmv < $max_distance and print STDOUT join("\n", @lls)."\n"; # for MAX NM:i:? 
	$nmv = -1; 
	$pre_query_name = ''; 
}


