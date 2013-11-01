#!/usr/bin/perl -w 
use strict; 

our $max_Xratio = 80/100; 
our $min_BaseNum = 40; 
#��fasta�ļ��е������������㣬[XxNn]�������С��һ��ratio�����߷�[XxNn]�����������һ������ֵ�����������output
 
if (@ARGV < 2)
{
  print "usage: trm_XNseq1.pl input1 input2 max_Xratio min_BaseNum > output\n";
  exit(0);
}

our $input1 = $ARGV[0];
our $input2 = $ARGV[1];
$max_Xratio = $ARGV[2];
$min_BaseNum = $ARGV[3];

my %seq_required;
my ($pre_name, $seq) = ('', ''); 
open(IN1, "$input1");
while (<IN1>) {
	if (/^>(\S+)/) {
		my $current_name = $1; 
		if ($pre_name ne '') {
			my $is_out = 1; #��ʾ������Ӧ��������������������
			$seq =~ s/\s//g; 
			my $xnnum = ($seq =~ tr/XxNn/XxNn/); #�õ���4�ֳ������[XxNn]�ļ������
			my $seqlen = length($seq); 
			if ($xnnum/$seqlen >= $max_Xratio) {#[XxNn]�����������һ��ratio
				$is_out = 0; #��ʾ�����в�Ӧ�����
			}elsif ($seqlen-$xnnum < $min_BaseNum) {#����[ATCG]�����������һ������
				$is_out = 0; #��ʾ�����в�Ӧ�����
			}else{
				#$seq =~ s/(.{50})/$1\n/g; chomp($seq); 
				#print ">$pre_name\n$seq\n"; #��������ľ��������׼���
				defined $seq_required{$pre_name} or $seq_required{$pre_name} = 1;#����query��query_length֮���ӳ��
			}
			#$is_out == 0 and warn "[Record] [$pre_name] dropped.\n"; #������
		}
		$pre_name=$current_name; $seq = ''; 
	}else{
		$seq .= $_; 
	}
}
#��Ҫ���Ǵ���ʣ�µ�
if ($pre_name ne '') {
	my $is_out = 1; 
	$seq =~ s/\s//g; 
	my $xnnum = ($seq =~ tr/XxNn/XxNn/);
	my $seqlen = length($seq);
	if ($xnnum/$seqlen >= $max_Xratio) {
		$is_out = 0;
	}elsif ($seqlen-$xnnum < $min_BaseNum) {
		$is_out = 0;
	}else{
		#$seq =~ s/(.{50})/$1\n/g; chomp($seq); 
		#print ">$pre_name\n$seq\n"; #��������ľ��������׼���
		defined $seq_required{$pre_name} or $seq_required{$pre_name} = 1;#����query��query_length֮���ӳ��
	}
	#$is_out == 0 and warn "[Record] [$pre_name] dropped.\n"; 
}
close(IN1);

#��input2�а���%seq_required�е�������ȡ����
open(IN2, $input2);

my $flag = "off";
while(<IN2>) {
	if($_ =~ m/^>/) {
		my $head = $_;
		chomp($head);
		$head=~s/>//;

		if(defined $seq_required{$head}) {#����������name
			print $_;#�����
			$flag = "on";#ͬʱ�ı��־����ʾ�����������Ҫ������OUT1���
		}
		else {#������������name
			#print OUT $_;#�������OUT
			$flag = "off";#ͬʱ�ı��־����ʾ�����������Ҫ������OUT2���
		}
	}
	else {
		if($flag eq "on") {#��ʾΪ"on"
			print $_;#�����������Ҫ�������
		}
	}
}
close(IN2);
