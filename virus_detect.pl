#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;
use FindBin;
my $usage = <<_EOUSAGE_;

#########################################################################################
# virus_detect.pl --file_type [String] --reference [FILE] --coverage [Float]
#                 --host_removal --host_reference [FILE] --objective_type [String]
#                 --max_dist[INT] --max_open [INT] --max_extension [INT] --len_seed [INT] --dist_seed [INT] --thread_num [INT]
#                 --strand_specific --min_overlap [INT] --max_end_clip [INT] --cpu_num [INT] --mis_penalty [INT] --gap_cost [INT] --gap_extension [INT]
#		  input_file1 input_file2 ... input_fileN
#  
# Options(6):
#  --file_type The format of input files(fastq or fasta)  [fastq]
#  --reference The name of a fasta file containing all of the virus reference sequences  [vrl_genbank.fasta]
#  --coverage  The reference sequence with at least a pencetage of its length covered by reads is detected as a aligned virus contig  [0.3]
#  --host_removal         Remove host sequences before assembly or not  [Not selected]
#  --host_reference       If using host_removal, the name of a host reference file should be provided  [Null]
#  --objective_type       The objective type of optimization(maxLen, n50 or avgLen)  [maxLen]
# 
# BWA-related options(6):
#  --max_dist      Maximum edit distance [1]  
#  --max_open      Maximum number of gap opens [1]  
#  --max_extension Maximum number of gap extensions [1]  
#  --len_seed      Take the first INT subsequence as seed [15] 
#  --dist_seed     Maximum edit distance in the seed [1]  
#  --thread_num    Number of threads (multi-threading mode) [8] 
# 
# Megablast-related options(7): for virus detection
#  --strand_specific Only for sequences assembled from strand specific RNA-seq  [Not selected]
#  --min_overlap     The minimum overlap length between two contigs to be combined [30]
#  --max_end_clip    The maximum length of end clips [6]
#  --cpu_num         Number of processors to use [8] 
#  --mis_penalty     Penalty for a nucleotide mismatch [-1]
#  --gap_cost        Cost to open a gap [2] 
#  --gap_extension   Cost to extend a gap [1] 
#
# Megablast-related options(): for virus identification
#  --word_size	     [11] 
#  --exp_value	     [1e-5]
#  --identity_percen [80] 
#  --mis_penalty_b   Penalty for a nucleotide mismatch [-1]
#  --gap_cost_b      Cost to open a gap [2] 
#  --gap_extension_b Cost to extend a gap [1]
#
###########################################################################################

_EOUSAGE_

	;

################################
##  set file folder path      ##
################################
our $WORKING_DIR  = cwd();				# set current folder as working folder
our $DATABASE_DIR = ${FindBin::RealBin}."/databases";	# set database folder
our $BIN_DIR	  = ${FindBin::RealBin}."/bin";		# set programs folder 

###############################
##   global vars	     ##
###############################
our $file_list;					#�������д�����������ļ������ݣ����б��ļ����޺�׺��
our $file_type= "fastq";			#�����ļ������ͣ���ǰֻ֧��fastq��fasta��ʽ
our $reference= "vrl_genbank.fasta";		#����ȫ���ο����е��ļ����ƣ�FASTA��ʽ��
our $coverage=0.3;  				#ÿ���ο����������reads���ǵĲ���ռȫ����������ֵ
our $host_removal;         			#�Ƿ���Ҫremove host
our $host_reference;       			#���Ҫremove host�������ṩһ���ļ�����host��ȫ���ο�����
our $objective_type='maxLen';			#�Ż�Ŀ��ֵ�����ͣ�ֻ��n50��maxLen��avgLen����

# paras for BWA

our $max_dist = 1;  				#bwa��������༭���� 
our $max_open = 1;  				#bwa��������gap����
our $max_extension = 1; 			#bwa��������gap����,-1��ʾ������gap
our $len_seed = 15; 				#bwa�е�����������
our $dist_seed = 1; 				#bwa��������������༭����
our $thread_num = 8; 				#bwa������õ��߳����� 

# paras for megablast detection

our $strand_specific;  				#ר������strand specificת¼��
our $min_overlap = 30; 				#hsp�ϲ�ʱ����̵�overlap
our $max_end_clip = 6; 				#hsp�ϲ�ʱ�������������Сclip
our $cpu_num = 8;          			#megablastʹ�õ�cpu��Ŀ
our $mis_penalty = -1;     			#megablast�У��Դ���ķ��֣������Ǹ�����
our $gap_cost = 2;         			#megablast�У���gap open�ķ��֣�������������
our $gap_extension = 1;    			#megablast�У���gap open�ķ��֣�������������

# paras for blast && identification 

our $diff_ratio= 0.25;
our $word_size = 11;
our $exp_value = 1e-5;				#
our $identity_percen = 25;			# tblastx�Ե������������ȶ�ʱhsp����Сͬһ��
our $mis_penalty_b = -1;			# megablast�У��Դ���ķ��֣������Ǹ�����
our $gap_cost_b = 2;				# megablast�У���gap open�ķ��֣�������������
our $gap_extension_b = 1;			# megablast�У���gap open�ķ��֣�������������

our $filter_query = "F";			# Ĭ�ϲ���Ҫȥ�������У��������Ҫ�û��趨
our $hits_return = 500;				# megablast���ص�hit��Ŀ���������Ҫ�û��趨

# disabled parameters
our $input_suffix='clean'; 			# �����ļ���׺��clean��unmapped��

# get input paras #
GetOptions(
	'file_type=s'	=> 	\$file_type,
	'reference=s'	=> 	\$reference,
	'coverage=f'	=> 	\$coverage,
	'host_removal!' => 	\$host_removal,
	'host_reference=s' => 	\$host_reference,
	'objective_type=s' => 	\$objective_type,
		
	'max_dist=i' => 	\$max_dist,
	'max_open=i' => 	\$max_open,
	'max_extension=i' => 	\$max_extension,
	'len_seed=i' => 	\$len_seed,
	'dist_seed=i' => 	\$dist_seed,			 
	'thread_num=i' => 	\$thread_num,
		
	'strand_specific!' => 	\$strand_specific,
	'min_overlap=i' => 	\$min_overlap,
	'max_end_clip=i' => 	\$max_end_clip,
	'cpu_num=i' => 		\$cpu_num ,
	'mis_penalty=i' => 	\$mis_penalty,
	'gap_cost=i' => 	\$gap_cost,
	'gap_extension=i' => 	\$gap_extension,


	'diff_ratio=s' => 	\$diff_ratio,
	'word_size=i' =>  	\$word_size,
	'exp_value-s' =>  	\$exp_value,
	'identity_percen=s' => 	\$identity_percen,	# tblastx�Ե������������ȶ�ʱhsp����Сͬһ��
	'mis_penalty_b=i' => 	\$mis_penalty_b,
	'gap_cost_b=i' => 	\$gap_cost_b,
	'gap_extension_b=i' => 	\$gap_extension_b
);

#die $usage unless $file_list;
my @filelist;
if ( !@ARGV )
{
	die "\nPlease input files:\n\n$usage\n";
}
else
{
	foreach my $m (@ARGV)
	{
		if (-s $m) {
			push(@filelist, $m);
		} else {
			print "Error, can not parse input file $m\n";
		}
	}
}

$file_list = "file_list_tp";
open(LS, ">".$file_list) || die $!;
if (scalar(@filelist) > 0) {
	foreach my $f (@filelist) {
		print LS $f."\n";
	}
} else {
	die "\nPlease input files:\n\n$usage\n";
}
close(LS);

# main
main: {
	# detect known virus
	# align small RNA reads to know plant virus sequence 
	print "alignend small RNA to know plant virus database: $reference\n";
	system("mkdir aligned");

	my $cmd_align = "$BIN_DIR/bwa-alignAndCorrect.pl --file_list $file_list --reference $DATABASE_DIR/$reference --coverage $coverage --max_dist $max_dist --max_open $max_open --max_extension $max_extension --len_seed $len_seed --dist_seed $dist_seed --thread_num $thread_num";
	print $cmd_align."\n";
	system($cmd_align) && die "Error in command $cmd_align\n";  

	if($strand_specific){
		system("$BIN_DIR/removeRedundancy_batch.pl --file_list $file_list --file_type $file_type --input_suffix contigs1.fa --contig_type aligned --contig_prefix KNOWN --min_overlap $min_overlap --max_end_clip $max_end_clip --cpu_num $cpu_num --mis_penalty $mis_penalty --gap_cost $gap_cost --gap_extension $gap_extension");
	}
	else{
		system("$BIN_DIR/removeRedundancy_batch.pl --file_list $file_list --file_type $file_type --input_suffix contigs1.fa --contig_type aligned --contig_prefix KNOWN --strand_specific --min_overlap $min_overlap --max_end_clip $max_end_clip --cpu_num $cpu_num --mis_penalty $mis_penalty --gap_cost $gap_cost --gap_extension $gap_extension");
	}
	
	#detect unknown virus
	print "\nassembly small RNA to contigs\n\n";
	system("mkdir assembled");

	if($host_removal){	
		system("$BIN_DIR/bwa_remove.pl --file_list $file_list --reference $DATABASE_DIR/$host_reference --max_dist $max_dist --max_open $max_open --max_extension $max_extension --len_seed $len_seed --dist_seed $dist_seed --thread_num $thread_num");
		$input_suffix='.unmapped'; #�Ժ�����ļ���׺������clean��Ϊunmapped
		#����Ҫע�⣬��sam�ļ�����ȡ�����ݶ���fastq�ļ�������ԭʼ�ļ���fastq����fasta
		system("$BIN_DIR/Velvet_Optimiser.pl --file_list $file_list --input_suffix $input_suffix  --file_type fastq --objective_type $objective_type --hash_end 19 --coverage_end 25");
		system("$BIN_DIR/runVelvet.pl --parameters optimization.result --input_suffix $input_suffix --file_type fastq --output_suffix contigs1.fa");
		system("rm *.unmapped");#�õ�contig1.fa���Ͳ���Ҫunmapped�ļ���
	}	
	else{
		system("$BIN_DIR/Velvet_Optimiser.pl --file_list $file_list --file_type $file_type --objective_type $objective_type --hash_end 19 --coverage_end 25");
		system("$BIN_DIR/runVelvet.pl --parameters optimization.result --file_type $file_type --output_suffix contigs1.fa");
	}

	print "\nremove redundancy for assemblied unknown small RNA\n\n";
	if($strand_specific){
	
		system("$BIN_DIR/removeRedundancy_batch.pl --file_list $file_list --file_type $file_type --input_suffix contigs1.fa --contig_type assembled --contig_prefix NOVEL --min_overlap $min_overlap --max_end_clip $max_end_clip --cpu_num $cpu_num --mis_penalty $mis_penalty --gap_cost $gap_cost --gap_extension $gap_extension");
	}
	else{
		system("$BIN_DIR/removeRedundancy_batch.pl --file_list $file_list --file_type $file_type --input_suffix contigs1.fa --contig_type assembled --contig_prefix NOVEL --strand_specific --min_overlap $min_overlap --max_end_clip $max_end_clip --cpu_num $cpu_num --mis_penalty $mis_penalty --gap_cost $gap_cost --gap_extension $gap_extension");
	}

	# combine the known and unknown virus
	print "combine the known and unkown contigs\n";
	system("mkdir combined");
	system("$BIN_DIR/files_combine2.pl --filelist $file_list --folder1 aligned --folder2 assembled");

	print "remove redundancy\n";
	system("$BIN_DIR/removeRedundancy_batch.pl --file_list $file_list --file_type $file_type --input_suffix contigs1.fa --contig_type combined --contig_prefix CONTIG --strand_specific --min_overlap $min_overlap --max_end_clip $max_end_clip --cpu_num $cpu_num --mis_penalty $mis_penalty --gap_cost $gap_cost --gap_extension $gap_extension");
	
	print "finished\n";

	unlink("optimization.result");
	unlink("velvet_assembly.result");
	#system("touch virus_detect.run.finished");	# create sign for the end of virus detection

	# identify the virus
	my $cmd_identify = "$BIN_DIR/virus_identify.pl ";
       	$cmd_identify .= "--file_list $file_list --file_type $file_type --reference $reference --contig_type combined ";
	$cmd_identify .= "--diff_ratio $diff_ratio --word_size $word_size --exp_value $exp_value --identity_percen $identity_percen ";
	$cmd_identify .= "--cpu_num $cpu_num --mis_penalty $mis_penalty_b --gap_cost $gap_cost_b --gap_extension $gap_extension_b ";
	system($cmd_identify) && die "Error in command: $cmd_identify\n";

	# delete list file
	unlink($file_list);
}


