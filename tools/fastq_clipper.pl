#!/usr/bin/perl
#use strict;   #�������붨��ſ���ʹ��
#use warnings; #��δ��ʼ�������Ǹ�������
use Getopt::Long; #���ģ����Խ�����������

#��������
my $usage = <<_EOUSAGE_;

#########################################################################################
# fastq_clipper.pl --filelist <FILE> --suffix <String> --adapter <String> --adapterLen INT --distance <INT> --levenshtein
# Required:
#  --filelist       a file containing all fastq file name for processing without suffix
#  --adapter        the reverse complimentary of 3' end adapter(5'->3')
#  --adapterLen     the adapter offset from the 5' end
#  --distance       the largest edit distance can be allowed between a adapter and reads
#  --levenshtein       using levenshtein distance, default is hamming distance
# Options(1):
#  --suffix	    the suffix of the input files
##########################################################################################

_EOUSAGE_
	;
#����������������в���	
our $filelist;        #������Ҫ���������fastq�ļ������б�
our $suffix="fastq";          #�����ļ���fastq��ʽ���ĺ�׺��
our $adapter;         #adapter����,���е����ݶ���ͬһ��adapter
our $adapter_offset;  
#ԭʼadapter��
#�е�ǰadapter_offset��bp���������adapterȫ������ƥ��
our $distance;        #adapter��readƥ��ʱ����������༭���룬Ŀǰ�����1
our $levenshtein=0;#����趨����ֵ�Զ�תΪ1����ʾʹ��levenshtein�㷨��༭����


#�������в����������ֵ
&GetOptions( 'filelist=s' => \$filelist,
             'adapter=s' => \$adapter,
	         'suffix=s' => \$suffix,
	         'adapterLen=i' => \$adapter_offset,
             'distance=i' => \$distance,
	         'levenshtein!' => \$levenshtein
			 );
unless ($filelist&&$adapter) {
	print $distance."\n";
	die $usage;
}
print "distance is $distance\n";

#����ȫ�ֱ���
our $adapter_substr = substr($adapter, 0, $adapter_offset);#ȡ��ʵ������ƥ���adapter���֣��ӵ�1���ַ���ʼ�����õ�$adapter_offset���ַ�

#������ʼ
open(IN1,$filelist) || die "Can't open the sample list file\n";
while(<IN1>){
	chomp;	
	#print $levenshtein."\n";#����ר��
	trim_fastq($_);	
}
close(IN1);

print "all the files have been processed\n";
#���������

#�ӳ������
###########################################################
#	Process fastq file                                #
###########################################################

sub trim_fastq{
	my $basename =shift;        #�õ���1�������������ļ�����
	my $input_file =$basename.".$suffix";                  #�����ļ�������
	my $output_file1 =$basename.".trimmed".$distance;      #����ļ�1������
    my $output_file2 =$basename.".unmatched".$distance;    #����ļ�2������
	my $output_file3 =$basename.".null".$distance;      #����ļ�3������
	
	my $match_position;     #@match_result�ĵ�2��Ԫ�أ�barcodeƥ�䵽read�ϵĵ�λ�ã�
	my $id1;                #fastq�ļ���ÿ��read�ĵ�1��
	my $DNAseq;             #fastq�ļ���ÿ��read�ĵ�2��
	my $id2;                #fastq�ļ���ÿ��read�ĵ�3��
	my $i = 0;              #���ڶ�fastq�ļ�ѭ���ڲ���i��ʾfastq�ļ��ĵ�i��
	
    open(IN2, $input_file) || die "Can't open the fastq file\n";
	open(OUT1, ">$output_file1") || die "can't create $output_file1 $!\n";#�ļ�������ļ���ӳ��
	open(OUT2, ">$output_file2") || die "can't create $output_file2 $!\n";#�ļ�������ļ���ӳ��
	open(OUT3, ">$output_file3") || die "can't create $output_file3 $!\n";#�ļ�������ļ���ӳ��
	while(<IN2>){
		chomp;
		$i ++;
		if($i%4 == 1){    #�����ǰ����ID��
			$id1 = $_;
		}
		if($i%4 == 2){    #�����ǰ����DNA����,ע��ֻ��������ܷ���ֵ
			$DNAseq = $_; #��ȡ��DNA����
			$match_position = adapter_match($DNAseq,$adapter_substr);#DNA������barcode����ƥ��
		}
		if($i%4 == 3){    #�����ǰ������һ��ID��
			$id2 = $_;
		}
		if($i%4 == 0){    #�����ǰ�����������У�4����������������ƥ��������ʼ���		  
			if ($match_position ==0){#��ʾ����DNA������Ⱦ�����浽".null"�ļ�
				print OUT3 $id1."\n".$DNAseq."\n".$id2."\n".$_."\n";#�����������
			}
			
			elsif($match_position ==length($DNAseq)){#û���ҵ�ƥ�䣬���浽".unmatched"�ļ�
				print OUT2 $id1."\n".$DNAseq."\n".$id2."\n".$_."\n";#�����������
			}	
			
			else{#trim��һ���֣����浽".trimmed"�ļ�
				my $trimmed_DNAseq = substr($DNAseq, 0, $match_position); #����DNA�����ϣ�barcodeƥ���֮ǰ�Ĳ���
		        my $trimmed_qual = substr($_, 0, $match_position);        #�������������ϣ�barcodeƥ���֮ǰ�Ĳ���
				print OUT1 $id1."\n".$trimmed_DNAseq."\n".$id2."\n".$trimmed_qual."\n";#���trimmed��ʣ�µ�����
			}		
		}#�������	
	}
	close(OUT1);
	close(OUT2);
	close(OUT3);
	close(IN2);
	my $trimmed_reads= count_reads($output_file1);            #ͳ�Ʋ���ʾfastq�ļ���reads����
	my $untrimmed_reads=count_reads($output_file2);           #ͳ�Ʋ���ʾtrim���ļ���reads��������������
	my $null_reads=count_reads($output_file3);
	print $basename."\t".$null_reads."\t".$trimmed_reads."\t".$untrimmed_reads."\n";#��ʾ�ļ�������Ϣ
	
}

sub adapter_match{
	my $read_seq = shift;    #�õ���1��������read����
	my $adapter_seq = shift; #�õ���2��������adapter����
	my $read_len = length($read_seq);	
	my $position = $read_len;	#����ֵ��ʼֵΪȫ������ʾû���ҵ�ƥ��
	my $min_offset=10;
 
	if($levenshtein==1){#�������levenshtein����
		for(my $i = 0; $i< $read_len-$adapter_offset+1; $i++){ 
			my $read_substr = substr($read_seq, $i, $adapter_offset);#������û����
			my $editDistance = levenshtein($read_substr,$adapter_seq);
			if ($editDistance <= $distance){
				$position=$i;  #�������λ�ã�����Ϊ����һ��trim
				last;          #�ҵ��ˣ��Ͳ��������ˣ�����ѭ��	
			}				
		}			
	}
	else{#�������hamming����
		for(my $i = 0; $i< $read_len-$min_offset+1; $i++){            
			my $read_substr = substr($read_seq, $i, $adapter_offset);#ע�����ʣ�µĳ��Ȳ���$adapter_offset���ж��ٽ�ȡ���ٶ�������
			my $editDistance = hamming($read_substr,$adapter_seq);#$read_substr�п��ܶ���$adapter_offset���̵ı�����ǰ��
			if ($editDistance <= $distance){				
				$position=$i;  #�������λ�ã�����Ϊ����һ��trim
				last;          #�ҵ��ˣ��Ͳ��������ˣ�����ѭ��	
			}			
		}
	}		
	return $position;
}
##############################################
#    ͳ�Ʋ���ʾfastq�ļ��е�read����         #
##############################################
sub count_reads{
        my $file = shift;
        my $count = `wc -l < $file`; #�õ��ļ�$file������
        die "wc failed: $?" if $?;
        chomp($count); 
        my $read = (int($count))/4;  #��������4���õ�read����
        return $read;
}
###############################################
#    ��2���ַ���֮���levenshtein����         #
###############################################
sub levenshtein
{
    my ($s1, $s2) = @_;
    my ($len1, $len2) = (length $s1, length $s2);
    return $len2 if ($len1 == 0);
    return $len1 if ($len2 == 0);
    my %mat;

    for (my $i = 0; $i <= $len1; ++$i)
    {
        for (my $j = 0; $j <= $len2; ++$j)
        {
            $mat{$i}{$j} = 0;
            $mat{0}{$j} = $j;
        }
        $mat{$i}{0} = $i;
    }

    my @ar1 = split(//, $s1);
    my @ar2 = split(//, $s2);

    for (my $i = 1; $i <= $len1; ++$i)
    {
        for (my $j = 1; $j <= $len2; ++$j)
        {
            my $cost = ($ar1[$i-1] eq $ar2[$j-1]) ? 0 : 1;
            $mat{$i}{$j} = min($mat{$i-1}{$j} + 1,
                                $mat{$i}{$j-1} + 1,
                                $mat{$i-1}{$j-1} + $cost);
        }
    }
    return $mat{$len1}{$len2};
}

sub min
{
    my @list = @_;
    my $min = $list[0];
    foreach my $i (@list)
    {
        $min = $i if ($i < $min);
    }
    return $min;
}
###########################################
#    ��2���ַ���֮���hamming����         #
###########################################
sub hamming($$) { length( $_[ 0 ] ) - ( ( $_[ 0 ] ^ $_[ 1 ] ) =~ tr[\0][\0] ) }
