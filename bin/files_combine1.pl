#!/usr/bin/perl

use strict;   
use warnings; 
use Getopt::Long; 
use IO::File;

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
#定义变量接受命令行参数	
our $filelist;        #包括需要处理的所有fastq文件名称列表
our $suffix;        #
our $startNumber = 1;  
our $endNumber; 

&GetOptions( 
	'filelist=s'	=> \$filelist,
	'suffix=s'	=> \$suffix,
	'startNumber=i' => \$startNumber,
	'endNumber=i'	=> \$endNumber
);

open(IN1,$filelist) || die "Can't open the $filelist file\n";
my $files1="";
while(<IN1>)
{
	chomp;
	my $file_base = $_;

	my $files="";
	
	for( my $i=$startNumber; $i<=$endNumber; $i++)
	{
		$files.=" ".$file_base.".".$suffix.$i." ";
	}
	
	my $finalFile = $_.".".$suffix;
	print "cat $files > $finalFile\n";
	my $result = `cat $files > $finalFile`;
}
close(IN1);

