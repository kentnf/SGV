#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Cwd;
use FindBin;
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
# global vars	#	
#################
our $file_list;
our $suffix="fastq";		# input file must have fastq suffix
our $run_R;
our $rRNA_removal;
our $virus_removal;
our $rRNA_reference="rRNA.fasta";
our $virus_reference= "virus_genbank186.fasta";	#virus sequence

our $adapter;			# small RNA adapter sequence
our $adapterLen=11;		# the leftmost length of adapter sequence
our $n_Cutoff=1;
our $read_Length=15; 
our $read_PerYield=5e5;		#5Mreads*4*100=2G bt
 
our $max_dist = 1;		#bwa max editable distance
our $max_open = 1;		#bwa max allowed gap count
our $max_extension = 1; 	#bwa max allowed gap length, -1 means not allowd long gap
our $len_seed=15; 		#bwa 中的种子区长度,100bp测序用50，50bp测序用36
our $dist_seed = 1; 		#bwa 种子区允许的最大编辑距离
our $thread_num = 8; 		#bwa cpu number

################################
# set folder and path	       #
################################
our $WORKING_DIR=cwd();				# current folder
our $DATABASE_DIR=$WORKING_DIR."/databases";	# database
our $BIN_DIR=$WORKING_DIR."/bin";		# bin script
our $Tools=$WORKING_DIR."/tools";		# tool script




##################
# input parameter#
##################
&GetOptions( 'file_list=s' => \$file_list,
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
# main program  #
#################
main: {

	my $test = ${FindBin::RealBin};
	print $test."\n";

	system("$Tools/fastq_clipper.pl --filelist $file_list --suffix $suffix --adapter $adapter --adapterLen $adapterLen --distance 1");
	system("$Tools/fastq_clipper.pl --filelist $file_list --suffix unmatched1 --adapter $adapter --adapterLen $adapterLen --distance 2");
	system("$BIN_DIR/files_combine1.pl --filelist $file_list --suffix trimmed --startNumber 1 --endNumber 2");
	#system("rm *.trimmed1");
	#system("rm *.trimmed2");
	#system("rm *.null1");
	#system("rm *.null2");
	#system("rm *.unmatched1");
	#system("rm *.unmatched2");

	# remove reads with more than 1 "N"
	# the remove reads short than 15 bp
	if($run_R){
		my $R_cmd = "Rscript $BIN_DIR/run_sRNA_clean.R filelist=$file_list nCutoff=$n_Cutoff readLength=$read_Length RdPerYield=$read_PerYield";
		system($R_cmd) && die "Error in command $R_cmd\n";
	}
      
	# remove rRNA reads
	if($rRNA_removal){
		my $BWA_remove_cmd = "$BIN_DIR/bwa_remove.pl --file_list $file_list --reference $DATABASE_DIR/$rRNA_reference --max_dist $max_dist --max_open $max_open --max_extension $max_extension --len_seed $len_seed --dist_seed $dist_seed --thread_num $thread_num";
		system($BWA_remove_cmd) && die "Error in command: $BWA_remove_cmd\n";

		my $file_change_cmd = "$BIN_DIR/files_name_change.pl --file_list $file_list --suffix1 unmapped --suffix2 clean";
		system($file_change_cmd) && die "Error in command : $file_change_cmd\n";
	}

	# remove virus reads
	if($virus_removal){	
		system("$BIN_DIR/bwa_remove.pl --file_list $file_list --reference $DATABASE_DIR/$virus_reference --max_dist $max_dist --max_open $max_open --max_extension $max_extension --len_seed $len_seed --dist_seed $dist_seed --thread_num $thread_num");
		system("$BIN_DIR/files_name_change.pl --file_list $file_list --suffix1 unmapped --suffix2 clean");
	}
}
