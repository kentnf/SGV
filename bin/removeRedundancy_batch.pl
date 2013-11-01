#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;
#��ͨ�����ѭ��ȥ���࣬Ȼ������base correction
my $usage = <<_EOUSAGE_;

#########################################################################################################
# removeRedundancy_batch.pl --file_list <FILE> --input_suffix <String> --contig_type <String> --contig_prefix <String> 
#                 --strand_specific --min_overlap [INT] --max_end_clip [INT] --cpu_num [INT] --mis_penalty [INT] --gap_cost [INT] --gap_extension [INT]
#
# Required(4):
#  --file_list       A txt file containing a list of input file names without any suffix
#  --input_suffix    A file name suffix of input data 
#  --contig_type     The contig type for redundancy removal (aligned or assembled)
#  --contig_prefix   A prefix of each contig name in the data (fasta format)
#
# Megablast-related options(7):
#  --strand_specific Only for sequences assembled from strand specific RNA-seq [Not selected]
#  --min_overlap     The minimum overlap length between two contigs to be combined [30]
#  --max_end_clip    The maximum length of end clips [4]
#  --cpu_num         Number of processors to use [8] 
#  --mis_penalty     Penalty for a nucleotide mismatch [-1]
#  --gap_cost        Cost to open a gap [2] 
#  --gap_extension   Cost to extend a gap [1]  
#########################################################################################################

_EOUSAGE_
	;
	
#################
##   ȫ�ֱ���  ##
#################
our $file_list; #�������д�������������ı��ļ����ƣ��޺�׺��
our $input_suffix;#�����ļ��ĺ�׺����
our $contig_type; #contig���ͣ�����������ļ�����aligned����assembled
our $contig_prefix;#fasta�ļ�����ʱ��ÿ����¼�õ���ǰ׺

our $strand_specific;  #ר������strand specificת¼��
our $min_overlap = 30; #hsp�ϲ�ʱ����̵�overlap
our $max_end_clip = 4; #hsp�ϲ�ʱ�������������Сclip
our $cpu_num = 8;          #megablastʹ�õ�cpu��Ŀ
our $mis_penalty = -1;     #megablast�У��Դ���ķ��֣������Ǹ�����
our $gap_cost = 2;         #megablast�У���gap open�ķ��֣�������������
our $gap_extension = 1;    #megablast�У���gap open�ķ��֣�������������
################################
##   ��������Ŀ¼���ļ���·�� ##
################################
our $WORKING_DIR=cwd();#����Ŀ¼���ǵ�ǰĿ¼
our $DATABASE_DIR=$WORKING_DIR."/databases";#�������ݿ��ļ����ڵ�Ŀ¼
our $BIN_DIR=$WORKING_DIR."/bin";#���п�ִ���ļ����ڵ�Ŀ¼
our $CONTIG_DIR;#ָ���Ĳ���contig���ڵ�Ŀ¼

#################
## �����������##
#################
&GetOptions( 'file_list=s' => \$file_list,
	'input_suffix=s' => \$input_suffix,
	'contig_type=s' => \$contig_type,
	'contig_prefix=s' => \$contig_prefix,
	'strand_specific!' => \$strand_specific,
	'min_overlap=i' => \$min_overlap,
	'max_end_clip=i' => \$max_end_clip,
	'cpu_num=i' => \$cpu_num,
	'mis_penalty=i' => \$mis_penalty,
	'gap_cost=i' => \$gap_cost,
	'gap_extension=i' => \$gap_extension			 			 
			 );
			 
unless ($file_list&&$input_suffix&&$contig_type&&$contig_prefix) {#��4������ͨ������õ�
die $usage;}

$CONTIG_DIR=$WORKING_DIR."/$contig_type";
#################
##  ������ʼ ##
#################
open(IN1,$file_list) || die "Can't open the file $file_list\n";
my $sample;
my $j=0;
while(<IN1>){
	$j=$j+1;
	chomp;
	$sample=$_; #ÿ��ѭ������һ�У��������붼�Ǵ���������ļ��������޺�׺����
	print "#processing sample $j by $0: $sample\n";
		
	my $file_size= -s "$CONTIG_DIR/$sample.$input_suffix";#���ݸ������ļ��ǲ���0���������洦������
	if($file_size==0){#����ļ���С��0�����������ָ��ֱ�ӽ�����һ��ѭ��
		system ("mv $CONTIG_DIR/$sample.$input_suffix $CONTIG_DIR/$sample.$contig_type.fa");
		next;
	}
	#�����ƶ�����ǰĿ¼��ȥ����
	system ("$BIN_DIR/dust $CONTIG_DIR/$sample.$input_suffix > $sample.masked");#����mask���ظ�����
	system ("$BIN_DIR/trim_XNseq1.pl $sample.masked $CONTIG_DIR/$sample.$input_suffix 0.8 40 > $sample.$input_suffix");
	system ("rm $sample.masked");
	system ("rm $CONTIG_DIR/$sample.$input_suffix");
	my $contig_num1 =  `grep \'>\' $sample.$input_suffix | wc -l `;#���ȵõ�ȥ����ǰ����������
	chomp($contig_num1);
	my $contig_num2 = 0;#�˱�������ÿ��ȥ����������������	
	my @a= split(/\./,$input_suffix);
	my ($i)= $a[0] =~ /(\d+)/;#��ȡcontigs��ʼ�ı��,Ĭ����1
	
	
	while($contig_num2!=$contig_num1)#�����contigs���������ھ��������������㣬ֱ������Ϊֹ
	{
		process_cmd ("$BIN_DIR/removeRedundancy.pl --input $sample.contigs$i.fa --min_overlap $min_overlap --max_end_clip $max_end_clip --cpu_num $cpu_num --mis_penalty $mis_penalty --gap_cost $gap_cost --gap_extension $gap_extension");		
		process_cmd("rm $sample.contigs$i.fa");#ɾ�����ļ�
		$i++;
		$contig_num1=$contig_num2;#����contig_num1
		#����contig_num2
		$contig_num2 =  `grep \'>\' inset | wc -l `;#У�������ɵ�contigs
		chomp($contig_num2);	
              	system ("mv inset $sample.contigs$i.fa");	
	}
	#ȥ������ϣ���ʼ����base correction
	my $sample_reference= $sample.".contigs$i.fa";#ȥ������Ľ���ļ���������contigs������Ϊ�ο�������У��
	my $sample_reads= $sample.".clean";#read�ļ�����Ҫalign����contigs��	
	#aligment -> sam -> bam -> sorted bam -> pileup
	&process_cmd("$BIN_DIR/bowtie-build -f $sample_reference $sample") unless (-e "$sample.1.amb");
	&process_cmd("$BIN_DIR/samtools faidx $sample_reference") unless (-e "$sample_reference.fai");
	&process_cmd("$BIN_DIR/bowtie $sample -v 1 -p $cpu_num -q $sample.clean -S -a --best $sample.sam") unless (-s "$sample.sam");
	&process_cmd("$BIN_DIR/samtools view -bt $sample_reference.fai $sample.sam > $sample.bam") unless (-s "$sample.bam");
	&process_cmd("$BIN_DIR/samtools sort $sample.bam $sample.sorted") unless (-s "$sample.sorted.bam");
	&process_cmd("$BIN_DIR/samtools mpileup -f $sample_reference $sample.sorted.bam > $sample.pileup") unless (-s "$sample.pileup");	

	$file_size= -s "$sample.pileup";#����pileup�ļ���С�ǲ���0���������洦������
	if($file_size==0){#����ļ���С��0����Ҫ����һ�����ļ���Ȼ������ѭ��
		system ("touch $CONTIG_DIR/$sample.$contig_type.fa");
		next;
	}
	#����������洦��
	$i++;
	&process_cmd("java -cp $BIN_DIR extractConsensus $sample 1 40 $i");
	system("$BIN_DIR/renameFasta.pl --inputfile $sample.contigs$i.fa --outputfile contigs --prefix $contig_prefix");#ÿ�����е�����Ҫͳһ������ȥ���ظ�
	system ("mv contigs $CONTIG_DIR/$sample.$contig_type.fa");#ȥ����õ��������ļ�ת�ص�$CONTIG_DIRĿ¼
	#ɾ������ѭ���в������м��ļ�
	system("rm $sample.sam");
	system("rm $sample.bam");
	system("rm $sample.sorted.bam");
	system("rm $sample.pileup");#ÿ�α���ɾ���������´β��ܼ���
	system("rm $sample_reference");
	system("rm $sample_reference.fai");
	system("rm *.ebwt");
	system("rm $sample.contigs$i.fa");	
}
close(IN1);
print "###############################\n";
print "All the samples have been processed by $0\n";
system("touch removeRedundancy_batch.$contig_type.finished");#��������ļ�����ʾ������־
system("rm tem.*");
system("rm tem");
system("rm *.log");
system("rm restset");
system("rm query");


#################
##    �ӳ���   ##
#################
sub process_cmd {
	my ($cmd) = @_;	
	print "CMD: $cmd\n";
	my $ret = system($cmd);	#�ɹ��ͷ���0�������ʧ��
	if ($ret) {
		die "Error, cmd: $cmd died with ret $ret";
	}
	return($ret);
}