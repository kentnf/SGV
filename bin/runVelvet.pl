#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case bundling);
use Bio::SeqIO;
use Cwd;
use FindBin;

my $usage = <<_EOUSAGE_;

#########################################################################################
# runVelvet.pl --velvet_dir [String] --parameters <File> --input_suffix <String>  --file_type <String> --output_suffix <String> 
# Required:
#  --parameters		A txt file containing a list of input file names with related parameters
#  --input_suffix	The suffix name of all the input files
#
# Options:
#  --velvet_dir		The directory to install the program velvet [workingdir/bin]
#  --file_type		The format of input files(fastq or fasta)  [fastq]
###########################################################################################


_EOUSAGE_
;

#################
# global vars   #
#################
our $velvet_dir;	# velvet的安装目录
our $parameters;	# 包括所有待处理样本和相应参数的文本文件名称（无后缀）
our $input_suffix = '';	# 输入文件的后缀名称
our $file_type="fastq";	# 输入文件的类型，当前只支持fastq和fasta格式
our $output_suffix; 	# 输出文件的后缀名称

################################
# set folder and file path     #
################################
my $WORKING_DIR = cwd();			# current folder
$velvet_dir = ${FindBin::RealBin};		# velvet path 
my $tf = $WORKING_DIR."/temp";			# temp folder

##################
# get input para #
##################
GetOptions( 
	'velvet_dir=s'		=> \$velvet_dir,
	'parameters=s' 		=> \$parameters,
	'input_suffix=s' 	=> \$input_suffix,
	'file_type=s' 		=> \$file_type,
	'output_suffix=s' 	=> \$output_suffix
);
	     
die $usage unless $parameters;

# main
main: {
        my $i=0;
	my $resultDir;

	open(IN, "$parameters");
	while (<IN>) {
		$i=$i+1;
		chomp;
		my @a = split(/\t/, $_);				# array: sample, Hash_len, Coverage_cutoff
		print "#processing sample $i: $a[0]\n";
		runVelvet($a[0],$a[1],$a[2]);				# assemblied reads into contigs using Velvet

		# move the assemblied contigs to ouput file
		$resultDir = $a[0]."_".$a[1]."_".$a[2];
		my $contigName = $a[0].".".$output_suffix;		# output file name
		system("mv $resultDir/contigs.fa $contigName");		# move assemblied contigs to output file
		#system("rm -r $resultDir");				# remove assemblied output folder

		my $count = `grep -c \'>\' $contigName `;		#  stat the number of contigs
		print $a[0]."\tAssemblied Contig No.:".$count;
	}
	close(IN);
	print "###############################\n";
	print "All the samples have been processed by $0\n";
}

# subroutine
sub runVelvet {
	my ($sample, $hash_length, $cov_cutoff) = @_;
	my $outputDir = $sample."_".$hash_length."_".$cov_cutoff;
	my $file;
	if ($input_suffix)	{ $file = "$sample.$input_suffix"; }
	else			{ $file = $sample; }
	process_cmd($velvet_dir."/velveth $outputDir $hash_length -$file_type $file > $tf/velvet.log");
	process_cmd($velvet_dir."/velvetg $outputDir -cov_cutoff $cov_cutoff -min_contig_lgth 30 > $tf/velvet.log");	
}

sub process_cmd {
	my ($cmd) = @_;	
	print "CMD: $cmd\n";
	my $ret = system($cmd);	
	if ($ret) {
		print "Error, cmd: $cmd died with ret $ret";
	}
	return($ret);
}
