#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;
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
#  --word_size [15] 
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
##   ��������Ŀ¼���ļ���·�� ##
################################
our $WORKING_DIR=cwd();#����Ŀ¼���ǵ�ǰĿ¼
our $DATABASE_DIR=$WORKING_DIR."/databases";#�������ݿ��ļ����ڵ�Ŀ¼
our $BIN_DIR=$WORKING_DIR."/bin";#���п�ִ���ļ����ڵ�Ŀ¼
our $seq_info= $DATABASE_DIR."/vrl_genbank.info";
our $result_dir= $WORKING_DIR."/result";

###############################
##   ȫ�ֱ�������Ĭ�����ã�  ##
###############################
our $file_list;#�������д�����������ļ������ݣ����б��ļ����޺�׺��
our $contig_type;#����contig������
our $reference= "vrl_genbank.fasta";#����ȫ���ο����е��ļ����ƣ�FASTA��ʽ��
our $file_type= "fastq";

our $diff_ratio= 0.2;
our $word_size = 15;
our $cpu_num = 8;          #megablastʹ�õ�cpu��Ŀ
our $mis_penalty = -1;     #megablast�У��Դ���ķ��֣������Ǹ�����
our $gap_cost = 2;         #megablast�У���gap open�ķ��֣�������������
our $gap_extension = 1;    #megablast�У���gap open�ķ��֣�������������
our $exp_value = 1e-5;    #
our $identity_percen = 25;    #tblastx�Ե������������ȶ�ʱhsp����Сͬһ��

our $filter_query = "F";     #Ĭ�ϲ���Ҫȥ�������У��������Ҫ�û��趨
our $hits_return = 500;   #megablast���ص�hit��Ŀ���������Ҫ�û��趨
our $input_suffix='clean';           #�����ļ���׺��clean��unmapped��
#################
## �����������##
#################
&GetOptions( 'file_list=s' => \$file_list,
	'contig_type=s' => \$contig_type,
	'file_type=s' => \$file_type,
	'reference=s' => \$reference,
	'diff_ratio=f' => \$diff_ratio,		
	'word_size=i' => \$word_size,
	'exp_value=f' => \$exp_value,
	'identity_percen=f' => \$identity_percen,
	'cpu_num=i' => \$cpu_num,
	'mis_penalty=i' => \$mis_penalty,
	'gap_cost=i' => \$gap_cost,
	'gap_extension=i' => \$gap_extension
	 );

unless ($file_list) {#������Ҫ1������
die $usage;
}

#################
##  ������ʼ ##
#################
main: {
my $sample;
my $i=0;
open(IN,$file_list) || die "Can't open the file $file_list\n";
while(<IN>){
	$i=$i+1;
	chomp;
	$sample=$_; #ÿ��ѭ������һ�У��������붼�Ǵ���������ļ��������޺�׺����
	print "#processing sample $i by $0: $sample\n";
	my $sample_dir=$result_dir."_$sample";
	system("mkdir $sample_dir");#ÿ����������һ���Լ��Ľ���ļ���
	#�ȵõ�coverage_ratio��Ϣ
	system("$BIN_DIR/formatdb -i $DATABASE_DIR/$reference -p F") unless (-e "$DATABASE_DIR/$reference.nhr");
	my $blast_program;
	my $blast_param;
	
	#�ȿ���û��known��virus
	$blast_program = $BIN_DIR."/megablast";
	print
	$blast_param = "-i ./$contig_type/$sample.$contig_type.fa -d $DATABASE_DIR/$reference -o $sample.blastn.paired -F $filter_query -a $cpu_num -W $word_size -q $mis_penalty -G $gap_cost -E $gap_extension -b $hits_return -e $exp_value";				
	&process_cmd($blast_program." ".$blast_param);
	&process_cmd("$BIN_DIR/blast_parse_table22.pl $sample.blastn.paired $sample.blastn.table");
	system("rm $sample.blastn.paired");#�����꣬����ɾ��
	#known��novel contig����
	&process_cmd("$BIN_DIR/query_filter1.pl $sample.blastn.table ./$contig_type/$sample.$contig_type.fa $sample.novel.contigs 80 25 > $sample.known.contigs");
	my $file_size= -s "$sample.known.contigs";#����pileup�ļ���С�ǲ���0���������洦������
	if($file_size==0){#����ļ���С��0���ͽ�������ѭ��
		system("rm $sample.known.contigs");
		system("rm $sample.blastn.table");
		next;
	}
	&process_cmd("$BIN_DIR/uniqComb.pl $sample.blastn.table -index $sample.known.contigs -col 0 -newCol 0 -exist | cut -f1-14 > $sample.known.table");
	&process_cmd("$BIN_DIR/hit_cov1.pl $sample.known.table $sample.known.cov 60 0.5 > $sample.known.block");
	&process_cmd("cut -f1 $sample.known.cov > hit_virus.list");
	&process_cmd("$BIN_DIR/extractFromFasta1.pl -i $DATABASE_DIR/$reference --type list --query hit_virus.list --output1 hit_virus.fa --output2 remainding.fa");
	system("rm $sample.known.block");#����ļ�ֻ����У��
	
	#�ٵõ�average depth��Ϣ
	my $sample_reads= $sample.".$input_suffix";#read�ļ�����Ҫalign����contigs��
	&process_cmd("bowtie-build hit_virus.fa virus");				
	if($file_type eq "fasta")
	{
		&process_cmd("bowtie virus -v 1 -p 8 -a --best --strata -f $sample_reads -S --sam-nohead $sample.sam");
	}
	else
	{
		&process_cmd("bowtie virus -v 1 -p 8 -a --best --strata $sample_reads -S --sam-nohead $sample.sam");	
	}		

	&process_cmd("$BIN_DIR/samtools faidx hit_virus.fa");	
	&process_cmd("$BIN_DIR/samtools view -bt hit_virus.fa.fai $sample.sam > $sample.bam");
	&process_cmd("$BIN_DIR/samtools sort $sample.bam $sample.sorted");
	&process_cmd("$BIN_DIR/samtools mpileup -f hit_virus.fa $sample.sorted.bam > $sample.pileup");	
	&process_cmd("$BIN_DIR/pileup_depth.pl $sample.pileup $sample.known.depth");
	&process_cmd("$BIN_DIR/ColLink.pl $sample.known.cov -keyC1 0 -Col1 1 -keyC2 0 -add -f1 $sample.known.depth > 1.tem");
	&process_cmd("$BIN_DIR/ColLink.pl 1.tem -keyC1 0 -Col1 2,3 -keyC2 0 -add -f1 $seq_info > 2.tem");		
	&process_cmd("sort -k6,6 -n -r 2.tem > $sample.known.identified1");
	&process_cmd("$BIN_DIR/hit_filter2.pl --input $sample.known.identified1 --diff_ratio $diff_ratio --output $sample.known.identified");#������������contig���ظ�������ϲ�hit
	&process_cmd("$BIN_DIR/query_filter2.pl $sample.known.identified $sample.known.table $sample.contigs.table");#�ѵ�һ�ļ��д��ڵ�hit��Ӧ��contig�������ݣ���ѡ����Ҫ����
	&process_cmd("$BIN_DIR/ColLink.pl $sample.contigs.table -keyC1 0 -Col1 2,3 -keyC2 2 -add -f1 $seq_info > $sample.contigs.info");
	&process_cmd("$BIN_DIR/ColLink.pl $sample.contigs.info -keyC1 0 -Col1 1 -keyC2 0 -add -f1 $sample.known.contigs > $sample.known.table1");
	&process_cmd("$BIN_DIR/arrange_col.pl $sample.known.table1 $sample_dir/$sample.known.xls");#���¶Ա�񲼾�
	system("mkdir $sample_dir/known_references");#����һ����Ŀ¼�����ڷ�����html�ļ�
	&process_cmd("$BIN_DIR/plot_results.pl $sample.known.identified $sample.blastn.table result_$sample known");

	system("rm $sample.sam");
	system("rm $sample.bam");
	system("rm $sample.sorted.bam");
	system("rm $sample.pileup");
	system("rm $sample.known.cov");
	system("rm $sample.known.depth");
	system("rm *.ebwt");
	system("rm $sample.blastn.table");
	system("rm $sample.known.contigs");
	system("rm $sample.known.table");
	system("rm $sample.known.identified1");
	system("rm $sample.known.identified");
	system("rm $sample.contigs.table");
	system("rm $sample.contigs.info");
	system("rm $sample.known.table1");
	
	if($contig_type ne "aligned")#�����Ҫ��һ�����novel�Ĳ���
	{
		#���Ȱ�novel��contigs��tblastx�Աȵ�������
		$blast_program = $BIN_DIR."/blastall -p tblastx";
		$blast_param = "-i $sample.novel.contigs -d $DATABASE_DIR/$reference -o $sample.novel.paired -F $filter_query -a $cpu_num -e $exp_value";				
		&process_cmd($blast_program." ".$blast_param);
		&process_cmd("$BIN_DIR/blast_parse_table44.pl $sample.novel.paired $sample.novel1.table");		
		system("rm $sample.novel.paired");
		my $file_size= -s "$sample.novel1.table";#���ݴ��ļ���С�ǲ���0���������洦������
		if($file_size==0){#����ļ���С��0���ͽ�������ѭ��
			system("rm $sample.novel.contigs");
			system("rm $sample.novel1.table");
			next;
		}
		&process_cmd("$BIN_DIR/blast_filter2.pl $sample.novel1.table $sample.novel.table");#��ÿ��query��hit�ԣ�ֻѡȡeֵ��ߵĽ��
		#����ע���hsp��Ҫ��Ҫ��
		&process_cmd("$BIN_DIR/hit_cov1.pl $sample.novel.table $sample.novel.cov $identity_percen 0.1 > $sample.novel.block");
		&process_cmd("cut -f1 $sample.novel.cov > hit_virus.list");
		&process_cmd("$BIN_DIR/extractFromFasta1.pl -i $DATABASE_DIR/$reference --type list --query hit_virus.list --output1 hit_virus.fa --output2 remainding.fa");
		system("rm $sample.novel.block");
		
		#�ٵõ�average depth��Ϣ
		my $sample_reads= $sample.".$input_suffix";#read�ļ�����Ҫalign����contigs��
		&process_cmd("bowtie-build hit_virus.fa virus");				
		if($file_type eq "fasta")
		{
			&process_cmd("bowtie virus -v 1 -p 8 -a --best --strata -f $sample_reads -S --sam-nohead $sample.sam");
		}
		else
		{
			&process_cmd("bowtie virus -v 1 -p 8 -a --best --strata $sample_reads -S --sam-nohead $sample.sam");	
		}		

		&process_cmd("$BIN_DIR/samtools faidx hit_virus.fa");	
		&process_cmd("$BIN_DIR/samtools view -bt hit_virus.fa.fai $sample.sam > $sample.bam");
		&process_cmd("$BIN_DIR/samtools sort $sample.bam $sample.sorted");
		&process_cmd("$BIN_DIR/samtools mpileup -f hit_virus.fa $sample.sorted.bam > $sample.pileup");	
		&process_cmd("$BIN_DIR/pileup_depth.pl $sample.pileup $sample.novel.depth");
		&process_cmd("$BIN_DIR/ColLink.pl $sample.novel.cov -keyC1 0 -Col1 1 -keyC2 0 -add -f1 $sample.novel.depth > 1.tem");
		&process_cmd("$BIN_DIR/ColLink.pl 1.tem -keyC1 0 -Col1 2,3 -keyC2 0 -add -f1 $seq_info > 2.tem");		
		&process_cmd("sort -k6,6 -n -r 2.tem > $sample.novel.identified1");
		&process_cmd("$BIN_DIR/hit_filter2.pl --input $sample.novel.identified1 --diff_ratio $diff_ratio --output $sample.novel.identified");#���ǵõ��������ļ�
		&process_cmd("$BIN_DIR/query_filter2.pl $sample.novel.identified $sample.novel.table $sample.contigs.table");#�ѵ�һ�ļ��д��ڵ�hit��Ӧ��contig��������
		&process_cmd("$BIN_DIR/ColLink.pl $sample.contigs.table -keyC1 0 -Col1 2,3 -keyC2 2 -add -f1 $seq_info > $sample.contigs.info");#���hit��ע����Ϣ
		&process_cmd("$BIN_DIR/fasta2tab.pl $sample.novel.contigs $sample.novel.contigs1");#��ʽת��
		&process_cmd("$BIN_DIR/ColLink.pl $sample.contigs.info -keyC1 0 -Col1 1 -keyC2 0 -add -f1 $sample.novel.contigs1 > $sample.novel.table1");#���contig��������Ϣ
		&process_cmd("$BIN_DIR/arrange_col.pl $sample.novel.table1 $sample_dir/$sample.novel.xls");#���¶Ա�񲼾�
		system("mkdir $sample_dir/novel_references");#����һ����Ŀ¼�����ڷ�����html�ļ�
		&process_cmd("$BIN_DIR/plot_results.pl $sample.novel.identified $sample.novel.table result_$sample novel");
		
		&process_cmd("cp ./$contig_type/$sample.$contig_type.fa $sample_dir/contig_sequences.fa");#���Ҫ������copy��ȥ
		system("rm $sample.sam");
		system("rm $sample.bam");
		system("rm $sample.sorted.bam");
		system("rm $sample.pileup");
		system("rm $sample.novel.cov");
		system("rm $sample.novel.depth");
		system("rm *.ebwt");	
		system("rm $sample.novel.contigs");
		system("rm $sample.novel.contigs1");
		system("rm $sample.novel.table");
		system("rm $sample.novel.identified1");
		system("rm $sample.novel.identified");
		system("rm $sample.contigs.table");#���ǹ�known������ͬ���ļ�
		system("rm $sample.contigs.info");#���ǹ�known������ͬ���ļ�
		system("rm $sample.novel.table1");
		system("rm $sample.novel1.table");		
	}
}
close(IN);
system("rm remainding.fa");
system("rm hit_virus.list");
system("rm hit_virus.fa");
system("rm hit_virus.fa.fai");
system("rm *.tem");
system("rm *.log");
print "###############################\n";
print "All the samples have been processed by $0\n";
system("touch virus_identify.$contig_type.finished");#��������ļ�����ʾ������־

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