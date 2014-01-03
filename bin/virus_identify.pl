#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;
use IO::File;
use Bio::SeqIO;
use FindBin;
use File::Basename;

my $usage = <<_EOUSAGE_;

#########################################################################################
# virus_itentify.pl --file_list <FILE> --contig_type <String> --file_type <String> --reference [FILE]
#                   --diff_ratio --word_size [INT] --exp_value <Float> --identity_percen <Float>
#                   --cpu_num  [INT] --mis_penalty [INT] --gap_cost[INT] --gap_extension [INT]
#                                  
# Required(2):
#  --file_list The name of a txt file containing a list of input file names without any suffix
#  --contig_type The type of contig files(aligned��assembled or combined) 
#
# Options(3):
#  --file_type The format of files containing reads are fasta or fastq  [fastq]
#  --reference The name of a fasta file containing all of the virus reference sequences  [vrl_genbank.fasta] 
#  --diff_ratio The hits with distance less than 0.2 will be combined into one  [0.2] 
# 
# blast-related options(7):
#  --word_size [11] 
#  --exp_value [1e-5]
#  --identity_percen [80] 
#  --cpu_num         Number of processors to use [8] 
#  --mis_penalty     Penalty for a nucleotide mismatch [-1]
#  --gap_cost        Cost to open a gap [2] 
#  --gap_extension   Cost to extend a gap [1] 
#
###########################################################################################

_EOUSAGE_

	;

################################
## set folder and file path   ##
################################
our $WORKING_DIR = cwd();				# current folder : working folder
our $BIN_DIR = ${FindBin::RealBin};			# bin folder
our $result_dir = $WORKING_DIR."/result";		# result folder
my  $tf = $WORKING_DIR."/temp";				# temp folder

our $DATABASE_DIR= $WORKING_DIR."/databases";		# database folder
our $seq_info = $DATABASE_DIR."/vrl_genbank.info";	# virus sequence info
###############################
##  global vars		     ##
###############################
our $file_list;			# temp list of input file
our $file_type= "fastq";
our $contig_type;		# aligned assembled conbined
our $reference= "";		# virus reference, fasta format

our $diff_ratio= 0.25;
our $word_size = 11;
our $cpu_num = 8;		# megablast: thread number
our $mis_penalty = -1;		# megablast: penalty for mismatch
our $gap_cost = 2;		# megablast: penalty for gap open
our $gap_extension = 1;		# megablast: penalty for gap extension
our $exp_value = 1e-5;		#
our $identity_percen = 25;	# tblastx�Ե������������ȶ�ʱhsp����Сͬһ��

our $filter_query = "F";	# megablast: F - disable remove simple sequence
our $hits_return = 500;		# megablast: hit number
our $input_suffix='';

########################
# get input parameters #
########################
GetOptions( 
	'file_list=s'	=> \$file_list,
	'contig_type=s' => \$contig_type,
	'file_type=s' 	=> \$file_type,
	'reference=s' 	=> \$reference,
	'diff_ratio=f' 	=> \$diff_ratio,		
	'word_size=i' 	=> \$word_size,
	'exp_value=f' 	=> \$exp_value,
	'identity_percen=f' => \$identity_percen,
	'cpu_num=i' 	=> \$cpu_num,
	'mis_penalty=i' => \$mis_penalty,
	'gap_cost=i' 	=> \$gap_cost,
	'gap_extension=i' => \$gap_extension
	 );

die $usage unless $file_list;	# required parameter

#################
#    main       #
#################
main: {

my $sample;
my $i=0;
my $format = "-q"; # default input file is fastq format

my $in =IO::File->new($file_list) || die "Can't open the list file $file_list $!\n";
while(<$in>)
{
	chomp;
	$i=$i+1; $sample=$_; 
	my $sample_base = basename($sample);
	print "#processing sample $i : $sample_base using $0\n";

	# create folder for each sample
	my $sample_dir = $result_dir."_".$sample_base;
	system("mkdir $sample_dir") unless -e $sample_dir;

	# copy aligned, assemblied, or combined contigs to each folder, named as contig_sequences
	my $contigs = $sample.".".$contig_type;
	process_cmd("cp $contigs $sample_dir/contig_sequences.fa");

	# format the reference plant virus sequence, blast it against contig sequences	
	system("$BIN_DIR/formatdb -i $DATABASE_DIR/$reference -p F") unless (-e "$DATABASE_DIR/$reference.nhr");
	my $blast_program = $BIN_DIR."/megablast";
	my $blast_param = "-i $contigs -d $DATABASE_DIR/$reference -o $sample.blastn.paired -F $filter_query -a $cpu_num -W $word_size -q $mis_penalty -G $gap_cost -E $gap_extension -b $hits_return -e $exp_value";				
	process_cmd($blast_program." ".$blast_param);

	# parse blast results to table format (NOT USE SEARCH::IO)
	# output table format:
	# query_name \t query_length \t hit_name \t hit_length \t hsp_length \t identity \t evalue \t score \t strand \t 
	# query_start \t query_end \t hit_start \t hit_end \t identity2 \t aligned_query \t aligned_hit \t aligned_string\n
	process_cmd("$BIN_DIR/blast_parse_table22.pl $sample.blastn.paired $sample.blastn.table1");

	# get paired query and hit according to e value	
	# the coverage of hsp should be more than 0.75 for query or hit
	# keep the best pair of <query,hit> according to evalue
	# һ��query��Ӧ��hits�� ֻ���������evalue�Ǹ�hsp��identity����drop_off(5%)���ڵ�
	#process_cmd("$BIN_DIR/blast_filter3.pl $sample.blastn.table1 $sample.blastn.table 0.75 5");
	# *** do not understand why using the last parameter:1 for this function	
	blast_filter("$sample.blastn.table1", "$sample.blastn.table", 0.75, 5, 1);

	# remove temp file
	#system("rm $sample.blastn.paired");
	#system("rm $sample.blastn.table1");
	
	# separate known and novel contig
	# 1. get the all the hsps with > 60 identity for one query
	# 2. combine all the hsps for one query
	# 3. get ratio of combined hsp and query
	# 4. separate know if ratio > 50 and novel contigs  
	process_cmd("$BIN_DIR/query_filter1.pl $sample.blastn.table $contigs $sample.novel.contigs 60 50 > $sample.known.contigs");

	my $file_size= -s "$sample.known.contigs";	#���ݴ��ļ���С�ǲ���0���������洦������
	if($file_size==0){				#����ļ���С��0���ͽ�������ѭ��
		#system("rm $sample.known.contigs");	#�����������ܹ�align����֪�������contigs��table��ʽ��
		#system("rm $sample.novel.contigs");	#���������в���align����֪�������contigs��fasta��ʽ��
		#system("rm $sample.blastn.table");
		system("touch $sample_dir/no_virus_detected");#��������ļ�����ʾ��������û�м�⵽����
		next;
	}

	# get the blast result (tab delimit format ) for known contigs
	process_cmd("$BIN_DIR/uniqComb.pl $sample.blastn.table -index $sample.known.contigs -col 0 -newCol 0 -exist > $sample.known.table");

	# get coverage for each hit  
	# process_cmd("$BIN_DIR/hit_cov1.pl $sample.known.table $sample.known.cov 60 0.5 > $sample.known.block");
	hit_cov("$sample.known.table", "$sample.known.cov", "$sample.known.block", 60, 0.5, 1);

	# get hit virus sequence ID
	process_cmd("cut -f1 $sample.known.cov > $tf/hit_virus.list");

	# generate hit virus sequence and remainding sequence
	process_cmd("$BIN_DIR/extractFromFasta.pl -i $DATABASE_DIR/$reference --type list --query $tf/hit_virus.list --output1 $tf/hit_virus.fa --output2 $tf/remainding.fa");

	# align reads to hit virus to get average depth
	my $sample_reads= $sample;
	process_cmd("bowtie-build --quiet $tf/hit_virus.fa $tf/virus");

	if($file_type eq "fasta"){$format="-f"};	
	process_cmd("bowtie --quiet $tf/virus -v 1 -p 8 -a --best --strata $format $sample_reads -S --sam-nohead $sample.sam");
	process_cmd("$BIN_DIR/samtools faidx $tf/hit_virus.fa");	
	process_cmd("$BIN_DIR/samtools view -bt $tf/hit_virus.fa.fai $sample.sam > $sample.bam");
	process_cmd("$BIN_DIR/samtools sort $sample.bam $sample.sorted");
	process_cmd("$BIN_DIR/samtools mpileup -f $tf/hit_virus.fa $sample.sorted.bam > $sample.pileup");	
	process_cmd("$BIN_DIR/pileup_depth.pl $sample.pileup $sample.known.depth");

	# extract the column 2, and 3 form sample.know.depth, 
	# add them to file sample.known.cov, generate new file 1.tem	
	process_cmd("$BIN_DIR/ColLink.pl $sample.known.cov -keyC1 0 -Col1 1,2 -keyC2 0 -add -f1 $sample.known.depth > $tf/1.tem");

	# extract the column 3, and 4 from seq_info
	# add them to file 1.tem, generate new file sample.known.identified1
	process_cmd("$BIN_DIR/ColLink.pl $tf/1.tem -keyC1 0 -Col1 2,3 -keyC2 0 -add -f1 $seq_info > $sample.known.identified1");	

	# combine hit accoridng to contig ? how it works
	process_cmd("$BIN_DIR/hit_filter2.pl --input $sample.known.identified1 --diff_ratio $diff_ratio --output $sample.known.identified");

	##��contig��Ӧ��hit���ˣ����߶�Ӧ��known.identified��hit�ϣ����߶�Ӧ��e value��ߵ�hit����	
	process_cmd("$BIN_DIR/query_filter2.pl $sample.known.identified $sample.known.table $sample.contigs.table");	

	#���$seq_info�ļ��еĵ�3��4�У�genus��spicies��
	process_cmd("$BIN_DIR/ColLink.pl $sample.contigs.table -keyC1 0 -Col1 2,3 -keyC2 2 -add -f1 $seq_info > $sample.contigs.info");

	#���$sample.known.contigs�ļ��еĵ�2�У�contig���У� 
	process_cmd("$BIN_DIR/ColLink.pl $sample.contigs.info -keyC1 0 -Col1 1 -keyC2 0 -add -f1 $sample.known.contigs > $sample.contigs.table1");

	#���¶Ա�񲼾֣��õ�known contigs�ıȶ��ļ�
	arrange_col1("$sample.contigs.table1", "$sample_dir/$sample_base.known.xls");

	#ȡ��hit name����3�У�,����������sam�ļ�
	process_cmd("cut -f1 $sample.known.identified > $tf/known.references.list");

	#ֻҪ.known.identified�е�hit��Ӧ��
	process_cmd("$BIN_DIR/uniqComb.pl $sample.contigs.table1 -index $tf/known.references.list -col 0 -newCol 2 -exist > $sample.known.table1");

	# #��ʾsam��	
	process_cmd("$BIN_DIR/extractFromFasta.pl -i $DATABASE_DIR/$reference --type list --query $tf/known.references.list --output1 $sample_dir/known.references.fa --output2 $tf/remainding.fa");

	#������Ӧhit��sam�ļ�
	process_cmd("$BIN_DIR/blastTable2sam.pl $sample.known.table1 > $sample_dir/$sample_base.known.sam");

	#����һ����Ŀ¼�����ڷ�����html�ļ�	
	system("mkdir $sample_dir/known_references");

	#��������$sample.known.cov
	# process_cmd("$BIN_DIR/hit_cov1.pl $sample.known.table1 $sample.known.cov 60 0.5 > $sample.known.block");
	hit_cov("$sample.known.table1", "$sample.known.cov", "$sample.known.block", 60, 0.5, 1);

	#��������$sample.known.identified1
	process_cmd("$BIN_DIR/ColLink.pl $sample.known.cov -keyC1 0 -Col1 6,7,8,9 -keyC2 0 -add -f1 $sample.known.identified > $sample.known.identified1");

	#��������$sample.known.identified����ʽ�б仯��
	arrange_col2("$sample.known.identified1",  "$sample.known.identified");
	process_cmd("$BIN_DIR/plot_results.pl $sample.known.identified $sample.known.table1 $sample_dir known");

	#system("rm $sample.blastn.table");
	#system("rm $sample.known.table");
	#system("rm $sample.known.cov");
	#system("rm $sample.known.block");#����ļ�ֻ����У��
	#system("rm $sample.sam");
	#system("rm $sample.bam");
	#system("rm $sample.sorted.bam");
	#system("rm $sample.pileup");
	#system("rm $sample.known.depth");
	#system("rm $sample.known.identified1");
	#system("rm $sample.known.identified");
	#system("rm $sample.contigs.table");
	#system("rm $sample.contigs.info");
	#system("rm $sample.contigs.table1");
	#system("rm known.references.list");
	#system("rm $sample.known.table1");
	#system("rm *.ebwt");

	# check novel virus of contig type is not aligned 
	if($contig_type ne "aligned")
	{
		# compare noval contigs against virus database using tblastx
		$blast_program = $BIN_DIR."/blastall -p tblastx";
		$blast_param = "-i $sample.novel.contigs -d $DATABASE_DIR/$reference -o $sample.novel.paired -F $filter_query -a $cpu_num -e $exp_value";				
		process_cmd($blast_program." ".$blast_param);
		process_cmd("$BIN_DIR/blast_parse_table44.pl $sample.novel.paired $sample.novel1.table");


		# get best paired query and hit according to e value	
		# the coverage of hsp should be more than 0 for query or hit
		# keep the best pair of <query,hit> according to evalue
		# һ��query��Ӧ��hits�� ֻ���������evalue�Ǹ�hsp��identity����drop_off(5%)���ڵ�
		#process_cmd("$BIN_DIR/blast_filter5.pl $sample.novel1.table $sample.novel.table 0 5");
		blast_filter("$sample.novel1.table", "$sample.novel.table", 0, 5, 3);
		#system("rm $sample.novel.paired");
		#system("rm $sample.novel1.table");

		# check the blast result
		# of the blast do not have any hit (file_size is 0), exit it
		my $file_size= -s "$sample.novel.table";
		if($file_size==0){
			#system("rm $sample.novel.table");
			#system("rm $sample.novel.contigs");			
			system("touch $sample_dir/no_novel_virus_detected"); # create this sign for no virus identified
			system("cut -f1 $sample.known.contigs > $tf/known.contigs.list");
			process_cmd("$BIN_DIR/extractFromFasta.pl -i $sample_dir/contig_sequences.fa --type list --query $tf/known.contigs.list --output1 $tf/known.contigs.fa --output2 $sample_dir/unknown.contigs.fa");
			#system("rm known.contigs.list");
			#system("rm $sample.known.contigs");
			#system("rm known.contigs.fa");
			next;
		}
		
		#����ע���hsp��Ҫ��Ҫ��
		# process_cmd("$BIN_DIR/hit_cov1.pl $sample.novel.table $sample.novel.cov $identity_percen 0 > $sample.novel.block");
		hit_cov1("$sample.novel.table", "$sample.novel.cov", "$sample.novel.block", $identity_percen, 0, 1);
		
		process_cmd("cut -f1 $sample.novel.cov > $tf/hit_virus.list");
		process_cmd("$BIN_DIR/extractFromFasta.pl -i $DATABASE_DIR/$reference --type list --query $tf/hit_virus.list --output1 $tf/hit_virus.fa --output2 $tf/remainding.fa");
		
		#�ٵõ�average depth��Ϣ
		my $sample_reads= $sample.$input_suffix;#read�ļ�����Ҫalign����contigs��
		process_cmd("bowtie-build --quiet $tf/hit_virus.fa $tf/virus");				
		process_cmd("bowtie --quiet $tf/virus -v 1 -p 8 -a --best --strata $format $sample_reads -S --sam-nohead $sample.sam");
		process_cmd("$BIN_DIR/samtools faidx $tf/hit_virus.fa");	
		process_cmd("$BIN_DIR/samtools view -bt $tf/hit_virus.fa.fai $sample.sam > $sample.bam");
		process_cmd("$BIN_DIR/samtools sort $sample.bam $sample.sorted");
		process_cmd("$BIN_DIR/samtools mpileup -f $tf/hit_virus.fa $sample.sorted.bam > $sample.pileup");	
		process_cmd("$BIN_DIR/pileup_depth.pl $sample.pileup $sample.novel.depth");


		#��depth�ļ���ȡ��2��3�мӵ��ļ�$sample.novel.cov�ĺ���
		process_cmd("$BIN_DIR/ColLink.pl $sample.novel.cov -keyC1 0 -Col1 1,2 -keyC2 0 -add -f1 $sample.novel.depth > $tf/1.tem");

		#��$seq_info��ȡ��3��4�мӵ��ļ�1.tem�ĺ���
		process_cmd("$BIN_DIR/ColLink.pl $tf/1.tem -keyC1 0 -Col1 2,3 -keyC2 0 -add -f1 $seq_info > $sample.novel.identified1");			
		#���ǵõ��������ļ�
		process_cmd("$BIN_DIR/hit_filter2.pl --input $sample.novel.identified1 --diff_ratio $diff_ratio --output $sample.novel.identified");

		#�ѵ�һ�ļ��д��ڵ�hit��Ӧ��contig��������
		process_cmd("$BIN_DIR/query_filter2.pl $sample.novel.identified $sample.novel.table $sample.contigs.table");

		#���hit��ע����Ϣ	
		process_cmd("$BIN_DIR/ColLink.pl $sample.contigs.table -keyC1 0 -Col1 2,3 -keyC2 2 -add -f1 $seq_info > $sample.contigs.info");

		# convert fasta sequence to tab delimit format
		# column1: seqID
		# column2: sequence
		fasta2tab("$sample.novel.contigs", "$sample.novel.contigs1");

		#���contig��������Ϣ
		process_cmd("$BIN_DIR/ColLink.pl $sample.contigs.info -keyC1 0 -Col1 1 -keyC2 0 -add -f1 $sample.novel.contigs1 > $sample.contigs.table1");		
		# re-generate table 
		#process_cmd("$BIN_DIR/arrange_col1.pl $sample.contigs.table1 > $sample_dir/$sample.novel.xls");
		arrange_col1("$sample.contigs.table1", "$sample_dir/$sample.novel.xls");	

		# get hit name (3rd col) for generate sam file	
		process_cmd("cut -f1 $sample.novel.identified > $tf/novel.references.list");

		#ֻҪ.known.identified�е�hit��Ӧ��
		process_cmd("$BIN_DIR/uniqComb.pl $sample.contigs.table1 -index $tf/novel.references.list -col 0 -newCol 2 -exist > $sample.novel.table1");

		#����һ����Ŀ¼�����ڷ�����html�ļ�	
		system("mkdir $sample_dir/novel_references");

		#��������$sample.novel.cov
		# process_cmd("$BIN_DIR/hit_cov2.pl $sample.novel.table1 $sample.novel.cov $identity_percen 0 > $sample.novel.block");
		hit_cov("$sample.novel.table1", "$sample.novel.cov", "$sample.novel.block", $identity_percen, 0, 3);


		#��������$sample.novel.identified1
		process_cmd("$BIN_DIR/ColLink.pl $sample.novel.cov -keyC1 0 -Col1 6,7,8,9 -keyC2 0 -add -f1 $sample.novel.identified > $sample.novel.identified1");

		# re-generate sample.novel.identified (format changed)
		#process_cmd("$BIN_DIR/arrange_col2.pl $sample.novel.identified1 > $sample.novel.identified");#��������$sample.novel.identified����ʽ�б仯��
		arrange_col2("$sample.novel.identified1", "$sample.novel.identified");
		process_cmd("$BIN_DIR/plot_results.pl $sample.novel.identified $sample.novel.table1 $sample_dir novel");
		
		#������ȡ�Ȳ���knownҲ����novel��contigs
		system("cut -f1 $sample.known.contigs > $tf/known.contigs.list");
		system("cut -f1 $sample.novel.table | sort | uniq > $tf/novel.contigs.list");
		system("cat $tf/known.contigs.list $tf/novel.contigs.list > $tf/final.contigs.list");#������hit������
		process_cmd("$BIN_DIR/extractFromFasta.pl -i $sample_dir/contig_sequences.fa --type list --query $tf/final.contigs.list --output1 $tf/known.contigs.fa --output2 $sample_dir/unknown.contigs.fa");
		
		#system("rm $sample.novel.table");
		#system("rm $sample.novel.cov");
		#system("rm $sample.novel.block");#����ļ�ֻ����У��
		#system("rm $sample.sam");
		#system("rm $sample.bam");
		#system("rm $sample.sorted.bam");
		#system("rm $sample.pileup");
		#system("rm $sample.novel.depth");
		#system("rm $sample.novel.identified1");
		#system("rm $sample.novel.identified");
		#system("rm $sample.contigs.table");#���ǹ�known������ͬ���ļ�
		#system("rm $sample.contigs.info");#���ǹ�known������ͬ���ļ�
		#system("rm $sample.contigs.table1");#���ǹ�known������ͬ���ļ�
		#system("rm novel.references.list");
		#system("rm $sample.novel.table1");
		#system("rm *.ebwt");
		#system("rm $sample.known.contigs");#�������ɾ��
		#system("rm $sample.novel.contigs");
		#system("rm $sample.novel.contigs1");	
		#system("rm known.contigs.list");
		#system("rm novel.contigs.list");
		#system("rm final.contigs.list");
		#system("rm known.contigs.fa");
	}
}
$in->close;
print "###############################\n";
print "All the samples have been processed by $0\n";
}

# put folder new folder

#################################################################
# kentnf: subroutine						#
#################################################################
sub process_cmd 
{
	my ($cmd) = @_;	
	print "CMD: $cmd\n";
	my $ret = system($cmd);	
	if ($ret) {
		print "Error, cmd: $cmd died with ret $ret";
	}
	return($ret);
}


=head2
 blast_filter : 
 blast_filter(input, output, coverage, identity_dropoff, ???param??? );
=cut
sub blast_filter
{
	my ($input, $output, $coverage, $identity_dropoff, $param) = @_;

	my $in = IO::File->new($input) || die $!;
	my $out = IO::File->new(">".$output) || die $!;

	my $last_query="";
	my $last_hit="";
	my $current_query;
	my $current_hit;
	my $high_identity;
	my $current_identity;

	while (<$in>) 
	{
		my @cols = split(/\t/, $_);

		if( $param*$cols[4]/$cols[1] >= $coverage || $param*$cols[4]/$cols[3] >= $coverage)	#query��hit����Ҫ��hsp����һ������
		{ 

			$current_query=$cols[0];
			$current_hit=$cols[2];	
			$current_identity=$cols[5];

			if($current_query ne $last_query){#�������һ���µ�query
				$high_identity=$cols[5];#��¼��ߵ�identity
				print $out $_;#query�ĵ�һ��hit�϶�Ҫ����Ϊ��evalue��͵�
			}else{
				if ($current_hit ne $last_hit && ($current_identity >= $high_identity-$identity_dropoff))#��ͬ��<query,hit>�У�ֻ����evalue��͵�(����һ��)
				{
					print $out $_;	
				}	
			}
			$last_query=$cols[0];
			$last_hit=$cols[2];
		}
	}
	
	$in->close;
	$out->close;
}

=head2
 hit_cov : 
 usage: hit_cov(input, output1, output2, identify, query_cov, output2, ???param???);

 # input file is the output of blast_parse_table2.pl
 # format of input file:
 # query_name \t query_length \t hit_name \t hit_length \t hsp_length \t identity \t evalue \t score \t	
 # strand \t query_start \t query_end \t hit_start \t hit_end ... \n
 # index 2��3��4��5��11��12 are required
 # for each pair of [query, hit] should meet the requirement of cutoff_identity and query_cov
 # ratio = all query coverage / hit_length
 # hit: query1, query2, ... queryN

=cut
sub hit_cov
{
	my ($input, $output1, $output2, $cutoff_identity, $query_cov, $param) = @_;

	# put hit info to hash hash
	# %bkl
	# key: hit_name, 
	# value: arrays of block and query
	# 	 for each element in array is another array include three element [hit_start, hit_end, hash_of_query]
	#        
	#        % hash_of_query
	#        key: queryID
	#        value: 1
	#
 	# %hit_len
	# key: hit_name
	# value: hit_length
	my %blk;
	my %hit_len;

	my $in = IO::File->new($input) || die $!;
	my $out1 = IO::File->new(">".$output1) || die $!;
	my $out2 = IO::File->new(">".$output2) || die $!;
	while (<$in>) 
	{
		chomp; 
		my @ta = split(/\t/, $_);
		die "Error in line $_" if scalar @ta < 13;
		my ($query_name, $query_length,	$hit_name, $hit_length,	$hsp_length, $identity,	$evalue,
		$score, $strand, $query_start, $query_end, $hit_start, $hit_end) = (@ta[0..12]);

		# covered = param * hsp length / query length
		my $query_covered= $param * $hsp_length / $query_length; 
		if ( $identity >= $cutoff_identity && $query_covered >= $query_cov)
		{
			#����hit��(hit_start,hit_end)֮���ӳ��
			push(@{$blk{$hit_name}}, [$hit_start, $hit_end]); 
			$blk{$ta[2]}[-1][2]{$ta[0]}=1; 

			defined $hit_len{$hit_name} or $hit_len{$hit_name} = $hit_length; 
		}
	}
	$in->close;

	#print OUT join("\t", qw/hit_name hit_len total_cov cov_rate% query_names/)."\n";
	print $out2 join("\t", qw/hit_name block_len block_start block_end covered_len query_names/)."\n"; 

	foreach my $tk (sort keys %blk) # sort by hit name 
	{
		my @o; #�洢û��overlap��hit�ϵ�block

				# sort by hit start, hit end, and query name
		foreach my $ar (sort {$a->[0]<=>$b->[0] || $a->[1]<=>$b->[1] || $a->[2] cmp $b->[2];} @{$blk{$tk}}) 
		{
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

		my $total_cov = 0; # hit totoal coverage 
		my %aa; 
		for my $ar (@o) {
			my @query_names = sort keys %{$ar->[2]};
			# output to file
			print $out2 join("\t", $tk, $hit_len{$tk}, $ar->[0], $ar->[1], $ar->[1]-$ar->[0]+1, join(",", @query_names))."\n";

			# hit name, hit len, hitһ��block����㣬hitһ��block���յ㣬hitһ��block�ĳ���
			$total_cov += ($ar->[1]-$ar->[0]+1); #һ��hit�����з��ص�block�ĳ���֮��
			@aa{@query_names} = (1) x scalar(@query_names); #������contigs�����ƴ���hash��
		}
		# output
		print $out1 join("\t", $tk, $hit_len{$tk}, $total_cov, $total_cov/$hit_len{$tk}*1.0, join(",", sort keys %aa), scalar(keys %aa))."\n"; 
		#hit���ƣ�hit���ȣ�hit������bp����hit�����ǵİٷֱȣ�����contigs���ƣ�contigs����
	}
	$out1->close;
	$out2->close;
}

=head2
 query_filter1 : 
 query_filter1(inputfile1, inputfile2, output1, output2, identity, min_ratio);
 # inputFile1 is the output of blast_parse_table2.pl
 # format of inputFile1
 # query_name \t query_length \t hit_name \t hit_length \t hsp_length \t identity \t evalue \t score \t strand \t 
 # query_start \t query_end \t hit_start \t hit_end \n
 # the index 0, 1, 9, 10 are required
 # all hit must meet the cutoff_identity
 # ratio = all hsp of hit / query_length
 # ��query�����ǳ���һ��ratio�ļ�¼������Ȼ���һ��fasta�ļ�����ȡ��¼�в���������Ӧ���У������output1
 # ʣ�µ����е����������ouput2
 # ���hit���ǲ��������������ratio���ж����õ�contig�ǲ��ǲ���
=cut
=head
sub query_filter1
{
	my ($input1, $input2, $output1, $output2, $identity, $min_ratio) = @_;

	my %blk; 
	my %query_len; 
	my %query_filtered;
	
	my $in1 = IO::File->new($input1) || die $!;
	while(<$in1>)
	{	
		chomp; 
		my @ta = split(/\t/, $_); 
		if ($ta[5]>=$cutoff_identity){			# only use record meet the cutoff_identity
			push(@{$blk{$ta[0]}}, [@ta[9,10]]); 	# create map between queryName and (query_start,query_end)
			defined $query_len{$ta[0]} or $query_len{$ta[0]} = $ta[1];
		}
	}
	$in1->close;

	# print join("\t", qw/query_name block_len block_start block_end covered_len/)."\n"; #���׼�������ĸ�������
	# print OUT join("\t", qw/query_name query_len total_cov cov_rate%/)."\n";   #�����ļ�����ĸ�������

	foreach my $tk (sort keys %blk) # sort by query name��Ȼ����ȡ������Ҫfiltered����query����
	{
		my @o; #�洢û��overlap��query�ϵ�block
		for my $ar (sort { $a->[0]<=>$b->[0] || $a->[1]<=>$b->[1];} @{$blk{$tk}}) 
		{
			if (scalar(@o) == 0) {
				push(@o, [@$ar]); 
			}elsif ($o[-1][1] >= $ar->[0]-1) {#��һ��query���ص�������ϲ�
				$o[-1][1] < $ar->[1] and $o[-1][1] = $ar->[1]; 
			}else{
				push(@o, [@$ar]); 
			}
		}

		my $total_cov = 0; 
		for my $ar (@o) {
			#print join("\t", $tk, $query_len{$tk}, $ar->[0], $ar->[1], $ar->[1]-$ar->[0]+1)."\n"; #���׼�������ĸ���
			#query���ƣ�query���ȣ�queryһ��block����㣬queryһ��block���յ㣬queryһ��block�ĳ���
			$total_cov += ($ar->[1]-$ar->[0]+1); #һ��query�����з��ص�block�ĳ���֮��
		}
		my $ratio=int($total_cov/$query_len{$tk}*10000+0.5)/100;
		if($ratio >= $min_ratio){#���з���������query,����Ҫ������
			defined $query_filtered{$tk} or $query_filtered{$tk} = 1;
		}
		#�����ļ�����ĸ���
		#query���ƣ�query���ȣ�query�����ǵ��ܳ��ȣ�query�����ǵİٷֱ�
	}

	#��input2�аѲ�������%query_filtered�е�������ȡ�����������OUT

	open(IN2, $input2);
	open(OUT, ">$output");
	my $flag = "off";
	while(<IN2>) {
		if($_ =~ m/^>/) {
			my $head = $_;
			chomp($head);
			$head=~s/>//;

			if(defined $query_filtered{$head}) {#����������name
				print $head."\t";#���������׼���
				$flag = "on";#ͬʱ�ı��־����ʾ�����������Ҫ������OUT1���
			}
			else {#������������name
				print OUT $_;#�������OUT
				$flag = "off";#ͬʱ�ı��־����ʾ�����������Ҫ������OUT2���
			}
		} else {
			if($flag eq "on") {#��ʾΪ"on"
				print $_;#�����������Ҫ�������׼������
			} else {
				print OUT $_;#���򣬺����������Ҫ������OUT���
			}
		}
	}
	close(IN2);
	close(OUT);
}
=cut

=head2
 fasta2tab: convert fasta sequence to tab delimit file
=cut
sub fasta2tab
{
	my ($input, $output) = @_;

	my $out = IO::File->new(">".$output) || die $!;
	my $in = Bio::SeqIO->new(-format=>'fasta', -file=>$input);
	while(my $inseq = $in->next_seq)
	{
		print $out $inseq->id."\t".$inseq->seq."\n";
	}
	$out->close;
}

=head2
 arrange_col1 : 

# ���ȵ����д�����Ȼ���������ӱ�ͷ���
# ����column��Ҫ�������,�ȸ���decription�������Hit_ID����
# Contig_ID	Contig_length	Hit_ID	Hit_length	strand	Contig_start	Contig_end	Hit_start	Hit_end	identity	genus	 description	contig_sequence
=cut
sub arrange_col1 
{
	my ($input, $output) = @_;

	my @all_data;
	my $in = IO::File->new($input) || die "Can not open input file $input $!\n";
	while(<$in>) {
		chomp;
		my @a = split(/\t/, $_);
		push(@all_data, [@a[0,19,1,2,3,17,18,9,10,11,12,13,6,8]]);	# re-order the data, then put them to array 	
	}
	$in->close;
	
	@all_data = sort { ($a->[5] cmp $b->[5]) || ($a->[3] cmp $b->[3])} @all_data; # sort according to Genus and hit

	my $out = IO::File->new(">".$output) || die "Can not open output file $output $!\n";	
	print $out "Contig_ID\tContig_Seq\tContig_Len\tHit_ID\tHit_Len\tGenus\tDescription\tContig_start\tContig_end\tHit_start\tHit_end\tHsp_identity\tE_value\tHsp_strand\n";
	foreach my $data ( @all_data ) {
		print $out join("\t", @$data)."\n";
	}
	$out->close;
}
=head2
 arrange_col2 : 
=cut
sub arrange_col2
{
	my ($input, $output) = @_;

	my @all_data;
	my $in = IO::File->new($input) || die "Can not open input file $input $!\n";
	while(<$in>) {
		chomp;
		my @a = split(/\t/, $_);
		my $coverage= 1.0*$a[7]/$a[1];					# get coverage
		push(@all_data, [@a[0,1,7],$coverage,@a[4,5,6,8,9]]);		# re-order the data, then put them to array
	}
	$in->close;
	
	@all_data = sort { ($a->[7] cmp $b->[7]) || ($b->[2] cmp $a->[2])} @all_data; # sort according to Genus and hit_covered(bp)

	my $out = IO::File->new(">".$output) || die "Can not open output file $output $!\n";	
	foreach my $data ( @all_data ) {
		print $out join("\t", @$data)."\n";
	}
	$out->close;	
}



