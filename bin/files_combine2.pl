#!/usr/bin/perl
#use strict;   #�������붨��ſ���ʹ��
#use warnings; #��δ��ʼ�������Ǹ�������
use Getopt::Long; #���ģ����Խ�����������

#��������
my $usage = <<_EOUSAGE_;

#########################################################################################
# files_combine2.pl --filelist <FILE> --folder1 <String> --folder2 <String>
# Required:
#  --filelist       a fastq file name without suffix for processing 

##########################################################################################

_EOUSAGE_
	;
#����������������в���	
our $filelist;        #������Ҫ���������fastq�ļ������б�
our $folder1;  
our $folder2; 
our $output_folder="combined";

#�������в����������ֵ
&GetOptions( 'filelist=s' => \$filelist,
	     'folder1=s' => \$folder1,
	     'folder2=s' => \$folder2
			 );

#������ʼ
open(IN1,$filelist) || die "Can't open the $filelist file\n";
while(<IN1>){
	chomp;
	$sample=$_; #ÿ��ѭ������һ�У��������붼�Ǵ���������ļ��������޺�׺����
	$file1="./".$folder1."/".$sample.".".$folder1.".fa";
	$file2="./".$folder2."/".$sample.".".$folder2.".fa";
	my $result_file = "./".$output_folder."/".$sample.".contigs1.fa";
	&process_cmd("cat $file1 $file2 > $result_file");
}
close(IN1);
####
sub process_cmd {
	my ($cmd) = @_;	
	print "CMD: $cmd\n";
	my $ret = system($cmd);	#�ɹ��ͷ���0�������ʧ��
	if ($ret) {
		die "Error, cmd: $cmd died with ret $ret";
	}
	return($ret);
}
