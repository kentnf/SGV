#!/usr/bin/perl
#use strict;   #�������붨��ſ���ʹ��
#use warnings; #��δ��ʼ�������Ǹ�������
use Getopt::Long; #���ģ����Խ�����������

#��������
my $usage = <<_EOUSAGE_;

#########################################################################################
# files_combine.pl --filelist <FILE> --suffix <String> --startNumber <INT> --endNumber <INT>
# Required:
#  --filelist       a fastq file name without suffix for processing 
#  --suffix
#  --number
##########################################################################################

_EOUSAGE_
	;
#����������������в���	
our $filelist;        #������Ҫ���������fastq�ļ������б�
our $suffix;        #
our $startNumber = 1;  
our $endNumber; 

#�������в����������ֵ
&GetOptions( 'filelist=s' => \$filelist,
			'suffix=s' => \$suffix,
			'startNumber=i' => \$startNumber,
			'endNumber=i' => \$endNumber
			 );

#������ʼ
open(IN1,$filelist) || die "Can't open the $filelist file\n";
my $files1="";
while(<IN1>){
	chomp;
	my $files="";
	
	for(my $i = $startNumber; $i<= $endNumber; $i++){
		$files=$files.$_.".$suffix$i ";#������Ҫ�ϲ����ļ���������һ��
	}
	
	my $finalFile=$_.".$suffix";
	print "cat $files > $finalFile\n";
	my $result = `cat $files> $finalFile`;		
}
close(IN1);

