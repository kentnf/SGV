#!/usr/bin/perl -w 
#findRedundancy(),��ȡһ��query��Ӧ������hit������hsp
use strict; 
use IO::File; 
use Getopt::Long; #���ģ����Խ�����������
use Cwd;
my $usage = <<_EOUSAGE_;

#########################################################################################
# hit_filter2.pl --input <FILE> --diff_ratio <FLOAT> --output <FILE>
#
# Required[1]:
#  --input              A name of an input file containing sequences in fasta format
#
# options[1]:
#   
##########################################################################################

_EOUSAGE_
	;
#################
##   ȫ�ֱ���  ##
#################	
our $input;            #��Ҫ������ļ�
our $diff_ratio;     
our $output;            

################################
##   ��������Ŀ¼���ļ���·�� ##
################################
our $WORKING_DIR=cwd();#����Ŀ¼���ǵ�ǰĿ¼
our $BIN_DIR=$WORKING_DIR."/bin";#���п�ִ���ļ����ڵ�Ŀ¼

##################
## ����������� ##
##################
&GetOptions( 'input=s' => \$input, 
	'diff_ratio=f' => \$diff_ratio,	
	'output=s' => \$output 	
			 );
unless ($input) {#�������ͨ������õ�
die $usage;}
			 
#################
##  ������ʼ ##
#################
			 
my (@all_data, $name, $sequence); 
$name = ''; $sequence = ''; 

open(IN, "$input");
while (<IN>) {
	chomp; 
	my @each_line = split(/\t/, $_);
        my @contigs= split(/,/, $each_line[4]);#��5��������contig������	
	push(@all_data, [scalar(@contigs), $each_line[3], $_]);#contig������covered���ʺ���������
}
close(IN);

my @inset;  #���ڴ洢���������У��������ٴδ洢@all_data�е����У�Ӧ�øĽ�Ϊֻ����index
my @restset = ''; #���ڴ洢�������У���Щ���ݵı����Ϊ��У����
@all_data = sort { -1*($a->[0] <=> $b->[0]) || -1*($a->[1] <=> $b->[1])} @all_data; #@all_data�����ݰ���contig������Ȼ��hit���Ƚ�������
=head;#�����һ��
for my $tr (@all_data) {
 print $tr->[2]."\n";
}
=cut;
my $contig_count=1; 
for my $tr (@all_data) {
	if (scalar(@inset)  == 0) {#��һ�ΰ���������ȷ���inset����
		push(@inset, $tr->[2]); 
	}else{
		my @aa = split(/\t/, $tr->[2]);
		#print $tr->[2]."\n";
		my $return_string = &ifRedundant(\@inset, \$aa[4]);
		if ($return_string eq "n") {#���query�Ƿ����࣬��������
			push(@inset, $tr->[2]); 
		}else{#���query�����࣬�ӵ�
			push(@restset, $tr->[2]); 
		}
	}
}

open(OUT1, ">$output");#�����з������������
for my $tr (@inset) {
	print OUT1 $tr."\n";
}
close(OUT1);

open(OUT2, ">restset");#�����������������
for my $tr (@restset) {
	print OUT2 $tr."\n";
}
close(OUT2);

#######################
##     �ӳ���ʼ    ##
#######################
sub ifRedundant {
	my ($inset, $query) = @_; 
	my @query_contigs= split(/,/, $$query);
	my $ratio;
	for my $tr (@$inset) {
		my @aa = split(/\t/, $tr);
		my @contigs= split(/,/, $aa[4]);
		my %contigs;
		for my $each_contig (@contigs) {#���Ȱѱ������е�
			$contigs{$each_contig}=1;
		}
		my $total_contigs=0;#ע��ÿ�α����ʼ��
		my $diff_contigs=0; #ע��ÿ�α����ʼ��
		for my $each_contig (@query_contigs) {
			$total_contigs++;
			if(not defined $contigs{$each_contig}){
				$diff_contigs++;			
			}
		}
		$ratio=$diff_contigs*1.0/$total_contigs;
		if($ratio<=$diff_ratio){return "r";} #һ������������contig������Ϊ�����࣬���̷���
	}
	return "n";#���inset�е��������ݱȽ�һ���û�з������Ƶģ�����Ϊ��������ļ�¼
}
sub process_cmd {
	my ($cmd) = @_;	
	print "CMD: $cmd\n";
	my $ret = system($cmd);	
	if ($ret) {
		die "Error, cmd: $cmd died with ret $ret";
	}
	return($ret);
}