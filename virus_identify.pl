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
#  --contig_type The type of contig files(aligned、assembled or combined) 
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
##   设置所有目录和文件的路径 ##
################################
our $WORKING_DIR=cwd();#工作目录就是当前目录
our $DATABASE_DIR=$WORKING_DIR."/databases";#所有数据库文件所在的目录
our $BIN_DIR=$WORKING_DIR."/bin";#所有可执行文件所在的目录
our $seq_info= $DATABASE_DIR."/vrl_genbank.info";
our $result_dir= $WORKING_DIR."/result";

###############################
##   全局变量（含默认设置）  ##
###############################
our $file_list;#包括所有待处理的输入文件（数据）的列表文件（无后缀）
our $contig_type;#输入contig的类型
our $reference= "vrl_genbank.fasta";#包括全部参考序列的文件名称（FASTA格式）
our $file_type= "fastq";

our $diff_ratio= 0.2;
our $word_size = 15;
our $cpu_num = 8;          #megablast使用的cpu数目
our $mis_penalty = -1;     #megablast中，对错配的罚分，必须是负整数
our $gap_cost = 2;         #megablast中，对gap open的罚分，必须是正整数
our $gap_extension = 1;    #megablast中，对gap open的罚分，必须是正整数
our $exp_value = 1e-5;    #
our $identity_percen = 25;    #tblastx以蛋白质序列来比对时hsp的最小同一性

our $filter_query = "F";     #默认不需要去除简单序列，这个不需要用户设定
our $hits_return = 500;   #megablast返回的hit数目，这个不需要用户设定
our $input_suffix='clean';           #数据文件后缀（clean或unmapped）
#################
## 输入参数处理##
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

unless ($file_list) {#至少需要1个参数
die $usage;
}

#################
##  主程序开始 ##
#################
main: {
my $sample;
my $i=0;
open(IN,$file_list) || die "Can't open the file $file_list\n";
while(<IN>){
	$i=$i+1;
	chomp;
	$sample=$_; #每次循环读入一行，后续代码都是处理该样本文件（名称无后缀）。
	print "#processing sample $i by $0: $sample\n";
	my $sample_dir=$result_dir."_$sample";
	system("mkdir $sample_dir");#每个样本都有一个自己的结果文件夹
	#先得到coverage_ratio信息
	system("$BIN_DIR/formatdb -i $DATABASE_DIR/$reference -p F") unless (-e "$DATABASE_DIR/$reference.nhr");
	my $blast_program;
	my $blast_param;
	
	#先看有没有known的virus
	$blast_program = $BIN_DIR."/megablast";
	print
	$blast_param = "-i ./$contig_type/$sample.$contig_type.fa -d $DATABASE_DIR/$reference -o $sample.blastn.paired -F $filter_query -a $cpu_num -W $word_size -q $mis_penalty -G $gap_cost -E $gap_extension -b $hits_return -e $exp_value";				
	&process_cmd($blast_program." ".$blast_param);
	&process_cmd("$BIN_DIR/blast_parse_table22.pl $sample.blastn.paired $sample.blastn.table");
	system("rm $sample.blastn.paired");#解析完，立刻删除
	#known与novel contig分离
	&process_cmd("$BIN_DIR/query_filter1.pl $sample.blastn.table ./$contig_type/$sample.$contig_type.fa $sample.novel.contigs 80 25 > $sample.known.contigs");
	my $file_size= -s "$sample.known.contigs";#根据pileup文件大小是不是0，进入下面处理流程
	if($file_size==0){#如果文件大小是0，就结束本次循环
		system("rm $sample.known.contigs");
		system("rm $sample.blastn.table");
		next;
	}
	&process_cmd("$BIN_DIR/uniqComb.pl $sample.blastn.table -index $sample.known.contigs -col 0 -newCol 0 -exist | cut -f1-14 > $sample.known.table");
	&process_cmd("$BIN_DIR/hit_cov1.pl $sample.known.table $sample.known.cov 60 0.5 > $sample.known.block");
	&process_cmd("cut -f1 $sample.known.cov > hit_virus.list");
	&process_cmd("$BIN_DIR/extractFromFasta1.pl -i $DATABASE_DIR/$reference --type list --query hit_virus.list --output1 hit_virus.fa --output2 remainding.fa");
	system("rm $sample.known.block");#这个文件只用于校对
	
	#再得到average depth信息
	my $sample_reads= $sample.".$input_suffix";#read文件，需要align到新contigs上
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
	&process_cmd("$BIN_DIR/hit_filter2.pl --input $sample.known.identified1 --diff_ratio $diff_ratio --output $sample.known.identified");#根据所包括的contig的重复情况，合并hit
	&process_cmd("$BIN_DIR/query_filter2.pl $sample.known.identified $sample.known.table $sample.contigs.table");#把第一文件中存在的hit对应的contig更新内容，并选择需要的列
	&process_cmd("$BIN_DIR/ColLink.pl $sample.contigs.table -keyC1 0 -Col1 2,3 -keyC2 2 -add -f1 $seq_info > $sample.contigs.info");
	&process_cmd("$BIN_DIR/ColLink.pl $sample.contigs.info -keyC1 0 -Col1 1 -keyC2 0 -add -f1 $sample.known.contigs > $sample.known.table1");
	&process_cmd("$BIN_DIR/arrange_col.pl $sample.known.table1 $sample_dir/$sample.known.xls");#重新对表格布局
	system("mkdir $sample_dir/known_references");#产生一个子目录，用于放最后的html文件
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
	
	if($contig_type ne "aligned")#如果需要进一步检查novel的病毒
	{
		#首先把novel的contigs用tblastx对比到病毒库
		$blast_program = $BIN_DIR."/blastall -p tblastx";
		$blast_param = "-i $sample.novel.contigs -d $DATABASE_DIR/$reference -o $sample.novel.paired -F $filter_query -a $cpu_num -e $exp_value";				
		&process_cmd($blast_program." ".$blast_param);
		&process_cmd("$BIN_DIR/blast_parse_table44.pl $sample.novel.paired $sample.novel1.table");		
		system("rm $sample.novel.paired");
		my $file_size= -s "$sample.novel1.table";#根据此文件大小是不是0，进入下面处理流程
		if($file_size==0){#如果文件大小是0，就结束本次循环
			system("rm $sample.novel.contigs");
			system("rm $sample.novel1.table");
			next;
		}
		&process_cmd("$BIN_DIR/blast_filter2.pl $sample.novel1.table $sample.novel.table");#对每个query和hit对，只选取e值最高的结果
		#这里注意对hsp的要求都要低
		&process_cmd("$BIN_DIR/hit_cov1.pl $sample.novel.table $sample.novel.cov $identity_percen 0.1 > $sample.novel.block");
		&process_cmd("cut -f1 $sample.novel.cov > hit_virus.list");
		&process_cmd("$BIN_DIR/extractFromFasta1.pl -i $DATABASE_DIR/$reference --type list --query hit_virus.list --output1 hit_virus.fa --output2 remainding.fa");
		system("rm $sample.novel.block");
		
		#再得到average depth信息
		my $sample_reads= $sample.".$input_suffix";#read文件，需要align到新contigs上
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
		&process_cmd("$BIN_DIR/hit_filter2.pl --input $sample.novel.identified1 --diff_ratio $diff_ratio --output $sample.novel.identified");#这是得到的最终文件
		&process_cmd("$BIN_DIR/query_filter2.pl $sample.novel.identified $sample.novel.table $sample.contigs.table");#把第一文件中存在的hit对应的contig保留下来
		&process_cmd("$BIN_DIR/ColLink.pl $sample.contigs.table -keyC1 0 -Col1 2,3 -keyC2 2 -add -f1 $seq_info > $sample.contigs.info");#添加hit的注释信息
		&process_cmd("$BIN_DIR/fasta2tab.pl $sample.novel.contigs $sample.novel.contigs1");#格式转换
		&process_cmd("$BIN_DIR/ColLink.pl $sample.contigs.info -keyC1 0 -Col1 1 -keyC2 0 -add -f1 $sample.novel.contigs1 > $sample.novel.table1");#添加contig的序列信息
		&process_cmd("$BIN_DIR/arrange_col.pl $sample.novel.table1 $sample_dir/$sample.novel.xls");#重新对表格布局
		system("mkdir $sample_dir/novel_references");#产生一个子目录，用于放最后的html文件
		&process_cmd("$BIN_DIR/plot_results.pl $sample.novel.identified $sample.novel.table result_$sample novel");
		
		&process_cmd("cp ./$contig_type/$sample.$contig_type.fa $sample_dir/contig_sequences.fa");#最后要把序列copy过去
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
		system("rm $sample.contigs.table");#覆盖过known产生的同名文件
		system("rm $sample.contigs.info");#覆盖过known产生的同名文件
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
system("touch virus_identify.$contig_type.finished");#建立这个文件，表示结束标志

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