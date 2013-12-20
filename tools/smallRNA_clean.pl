#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;
my $usage = <<_EOUSAGE_;

#########################################################################################
# smallRNA_clean.pl --file_list <FILE>  --suffix [String] --adapter <String> --adapterLen [String]
# 		    --run_R --n_Cutoff [INT] --read_Length [INT] --read_PerYield [FLOAT] 
#                   --rRNA_removal --rRNA_reference [FILE] --virus_removal --virus_reference [FILE]
#                   --max_dist[INT] --max_open [INT] --max_extension [INT] --len_seed [INT] --dist_seed [INT] --thread_num [INT]
#                                  
# Required(2):
#  --file_list  The name of a txt file containing a list of input file names without any suffix 
#  --suffix     The suffix of input files
#
# Options(7):
#  --run_R            run R script to remove low quality, ambiguous and adapter bases [Not selected]
#  --rRNA_removal     Remove rRNA sequences [Not selected]
#  --rRNA_reference   a fasta file containing the rRNA reference sequences
#  --virus_removal    Remove virus sequences [Not selected]
#  --virus_reference  a fasta file containing the virus reference sequences
#  --adapter          the reverse complimentary of 3' end adapter(5'->3')
#  --adapterLen       the adapter offset from the 5' end [11]
#
# R script-related options(3): 
#  --n_Cutoff         after low quality bases trimmed, the reads containing "N" with a number above this parameter will be removed [1]  
#  --read_Length      reads with the length shorter than this parameter will be removed [15] 
#  --read_PerYield    Number of reads in each iteration to limit the memory usage [5e5] 
#
# BWA-related options(6):
#  --max_dist      Maximum edit distance [1]  
#  --max_open      Maximum number of gap opens [1]  
#  --max_extension Maximum number of gap extensions [1] 
#  --len_seed      Take the first INT subsequence as seed for the BWA alignment [15] 
#  --dist_seed     Maximum edit distance in the seed [1]  
#  --thread_num    Number of threads (multi-threading mode) [8] 
###########################################################################################

_EOUSAGE_

	;

#################
##   ȫ�ֱ���  ##
#################
our $file_list;
our $suffix="fastq";;          #�����ļ���fastq��ʽ���ĺ�׺��
our $run_R;
our $rRNA_removal;
our $virus_removal;
our $rRNA_reference="rRNA.fasta";
our $virus_reference= "virus_genbank186.fasta";#����ȫ���ο����е��ļ����ƣ�FASTA��ʽ��

our $adapter;         #adapter����,���е����ݶ���ͬһ��adapter
our $adapterLen=11;  #ԭʼadapter���е�ǰadapter_offset��bp���������adapterȫ������ƥ��
our $n_Cutoff=1;
our $read_Length=15; 
our $read_PerYield=5e5;#5Mreads*4*100=2G�ֽ�
 
our $max_dist = 1;  #bwa��������༭���� 
our $max_open = 1;  #bwa��������gap����
our $max_extension = 1; #bwa��������gap����,-1��ʾ������gap
our $len_seed=15; #bwa�е�����������,100bp������50��50bp������36
our $dist_seed = 1; #bwa��������������༭����
our $thread_num = 8; #bwa������õ��߳����� 

################################
##   ��������Ŀ¼���ļ���·�� ##
################################
our $WORKING_DIR=cwd();#����Ŀ¼���ǵ�ǰĿ¼
our $DATABASE_DIR=$WORKING_DIR."/databases";#�������ݿ��ļ����ڵ�Ŀ¼
our $BIN_DIR=$WORKING_DIR."/bin";#���п�ִ���ļ����ڵ�Ŀ¼
our $Tools=$WORKING_DIR."/tools";#������Ϊ���ߵĿ�ִ���ļ����ڵ�Ŀ¼

##################
## ����������� ##
##################
&GetOptions( 'file_list=s' => \$file_list,#�������д�����������ļ����ƣ��޺�׺��
	'suffix=s' => \$suffix,
	'adapter=s' => \$adapter,
	'adapterLen=i' => \$adapterLen,
	'run_R!' => \$run_R,	
	'n_Cutoff=i' => \$n_Cutoff,
	'read_Length=i' => \$read_Length,
	'read_PerYield=f' => \$read_PerYield,	
	'rRNA_removal!' => \$rRNA_removal,
	'rRNA_reference=s' => \$rRNA_reference,
	'virus_removal!' => \$virus_removal,	
	'virus_reference=s' => \$virus_reference,
	'max_dist=i' => \$max_dist,
	'max_open=i' => \$max_open,
	'max_extension=i' => \$max_extension,
	'len_seed=i' => \$len_seed,
	'dist_seed=i' => \$dist_seed,			 
	'thread_num=i' => \$thread_num
			 );

unless ($file_list) {
	die $usage;
}

#################
##  ������ʼ ##
################
main: {
    #�ⲿ����clipper���̣��������Ҫ����ע�͵�
	system("$Tools/fastq_clipper.pl --filelist $file_list --suffix $suffix --adapter $adapter --adapterLen $adapterLen --distance 1");
	system("$Tools/fastq_clipper.pl --filelist $file_list --suffix unmatched1 --adapter $adapter --adapterLen $adapterLen --distance 2");
	system("$BIN_DIR/files_combine1.pl --filelist $file_list --suffix trimmed --startNumber 1 --endNumber 2");
	system("rm *.trimmed1");
	system("rm *.trimmed2");
	system("rm *.null1");
	system("rm *.null2");
	system("rm *.unmatched1");
	system("rm *.unmatched2");
    #��ȥ������1��"N"���ϵ�reads�����ȥ����reads��<15bp��
	if($run_R){
		system("Rscript $BIN_DIR/run_sRNA_clean.R filelist=$file_list nCutoff=$n_Cutoff readLength=$read_Length RdPerYield=$read_PerYield");
 	}       
	#�����Ҫ����ȥ��rRNA��Ⱦ
	if($rRNA_removal){	
		system("$BIN_DIR/bwa_remove.pl --file_list $file_list --reference $DATABASE_DIR/$rRNA_reference --max_dist $max_dist --max_open $max_open --max_extension $max_extension --len_seed $len_seed --dist_seed $dist_seed --thread_num $thread_num");
		system("$BIN_DIR/files_name_change.pl --file_list $file_list --suffix1 unmapped --suffix2 clean");
	}
	#�����Ҫ����ȥ��virus��Ⱦ
	if($virus_removal){	
		system("$BIN_DIR/bwa_remove.pl --file_list $file_list --reference $DATABASE_DIR/$virus_reference --max_dist $max_dist --max_open $max_open --max_extension $max_extension --len_seed $len_seed --dist_seed $dist_seed --thread_num $thread_num");
		system("$BIN_DIR/files_name_change.pl --file_list $file_list --suffix1 unmapped --suffix2 clean");
	}
}