#!/usr/bin/perl -w 
#findRedundancy(),��ȡһ��query��Ӧ������hit������hsp
use strict; 
use IO::File; 
use Getopt::Long; #���ģ����Խ�����������
use Cwd;
my $usage = <<_EOUSAGE_;

#########################################################################################
# removeRedundancy.pl --input <FILE> --strand_specific --min_overlap [INT] --max_end_clip [INT] 
#                      --cpu_num [INT] --mis_penalty [INT] --gap_cost [INT] --gap_extension [INT]
# Required[1]:
#  --input              A name of an input file containing sequences in fasta format
#
# Megablast-related options[7]:
#  --strand_specific    Only for sequences assembled from strand specific RNA-seq [Not selected]
#  --min_overlap        The minimum overlap length between two contigs to be combined [30]
#  --max_end_clip       The maximum length of end clips [4]
#  --cpu_num         Number of processors to use [8] 
#  --mis_penalty     Penalty for a nucleotide mismatch [-1]
#  --gap_cost        Cost to open a gap [2] 
#  --gap_extension   Cost to extend a gap [1]     
##########################################################################################

_EOUSAGE_
	;
#################
##   ȫ�ֱ���  ##
#################	
our $input;            #��Ҫ������ļ�
our $strand_specific;  #ר������strand specificת¼��
our $min_overlap = 30; #hsp�ϲ�ʱ����̵�overlap
our $max_end_clip = 4; #hsp�ϲ�ʱ�������������Сclip
our $min_identity;     #hsp�ϲ���Ҫ�������Сidentity����hsp���ȶ�����

our $filter_query = "F";     #Ĭ�ϲ���Ҫȥ��������
our $word_size = int($min_overlap/3);
our $cpu_num = 8;          #megablastʹ�õ�cpu��Ŀ
our $hits_return = 10;     #megablast���ص�hit��Ŀ���������Ҫ�û��趨
our $mis_penalty = -1;     #megablast�У��Դ���ķ��֣������Ǹ�����
our $gap_cost = 2;         #megablast�У���gap open�ķ��֣�������������
our $gap_extension = 1;    #megablast�У���gap open�ķ��֣�������������

################################
##   ��������Ŀ¼���ļ���·�� ##
################################
our $WORKING_DIR=cwd();#����Ŀ¼���ǵ�ǰĿ¼
our $BIN_DIR=$WORKING_DIR."/bin";#���п�ִ���ļ����ڵ�Ŀ¼

##################
## ����������� ##
##################
&GetOptions( 'input=s' => \$input, 
	'strand_specific!' => \$strand_specific,
	'min_overlap=i' => \$min_overlap,
	'max_end_clip=i' => \$max_end_clip,
	'cpu_num=i' => \$cpu_num ,
	'mis_penalty=i' => \$mis_penalty,
	'gap_cost=i' => \$gap_cost,
	'gap_extension=i' => \$gap_extension											 
			 );
unless ($input) {#�������ͨ������õ�
die $usage;}
			 
#################
##  ������ʼ ##
#################
			 
my (@all_data, $name, $sequence); 
$name = ''; $sequence = ''; 

open(IN, "$input");
while (<IN>) {
	chomp; 
	s/\s+$//; 
	if (/^>(\S+)/) {
		my $tk = $1; 
		if ($name ne '') {#����Ѿ��������������ƣ���ǰһ�εģ������������ų���һ�Σ�ǰһ��Ϊ�գ����
			$sequence =~ s/\s//g; #ȥ��sequence�еĿ��ַ�
			push(@all_data, [$name, length($sequence), $sequence]);
		}
		$name = $tk; #��ʼ����һ�������У����ȶ�����������
		$sequence = ''; #���������
	}else{
		$sequence .= $_; 
	}
}
close(IN);

if ($name ne '') { #���һ�ζ����fasta���ݣ���������
	$sequence =~ s/\s//g;
	push(@all_data, [$name, length($sequence), $sequence]); 
}
my %inset;  #���ڴ洢���������У��������ٴδ洢@all_data�е����У�Ӧ�øĽ�Ϊֻ����index
my $restset = ''; #���ڴ洢�������У���Щ���ݵı����Ϊ��У����
@all_data = sort { -1*($a->[1] <=> $b->[1]) } @all_data; #@all_data�����ݰ������г��Ƚ�������

my $contig_count=1; 
for my $tr (@all_data) {
	if (scalar(keys %inset)  == 0) {#��һ�ΰ���������ȷ���inset����
		$inset{$tr->[0]} = $tr->[2]; 
	}else{
		my $query = ">$tr->[0]\n$tr->[2]\n"; 
		my $return_string = &ifRedundant(\%inset, \$query);
		my @return_col = split(/\t/, $return_string);#��1����name����2����r��n�����Ǻϲ�������
		if ($return_col[1] eq "r") {#���query����ȫ���࣬�ӵ�
			$restset .= $query;
		}elsif($return_col[1] eq "n"){#���query�Ƿ����࣬��������
			$contig_count++;
			$inset{$tr->[0]} = $tr->[2]; 
		}
		else{#���query�ǲ������࣬�ϲ����滻
			my @names = split(/\:/, $return_col[0]);#��һ����hit name : query name		
			delete($inset{$names[0]}); #ԭ��hit name��Ӧ�ļ�¼ɾ��
			my $new_name =  $return_col[0];
			$new_name =~ s/\:/\+/;	
			$inset{$new_name} = $return_col[1];#����������һ����¼
			$restset .= $query;
		}
	}
}

open(OUT1, ">inset");#�����з������������
while (my ($k,$v) = each %inset) {
	print OUT1 ">$k\n$v\n"; 
}
close(OUT1); 

open(OUT2, ">restset");#�����������������
print OUT2 $restset; 
close(OUT2);

print "@@\t".$input."\t".$contig_count."\n";


#######################
##     �ӳ���ʼ    ##
#######################
sub ifRedundant {
	my ($inset, $query) = @_; 
	
	#�������б���ѯ���е��ļ�tem
	open(OUT, ">tem");
	while (my ($k,$v) = each %$inset) {
	print OUT ">$k\n$v\n"; 
	}
	close(OUT);
	
	#�����ѯ���е��ļ�query
	open(OUT, ">query");
	print OUT $$query; 
	close(OUT);
	
    #ִ��blast�������Ҫ�����������process_cmd()
	system("$BIN_DIR/formatdb -i tem -p F");#tem��inset�ļ��������������	
	my $blast_program = $BIN_DIR."/megablast";
	#my $blast_param = "-i query -d tem -o tem.paired -F F -a 8 -W $word_size -q -1 -G 2 -E 1 -b 5";
	my $blast_param = "-i query -d tem -o tem.paired -F $filter_query -a $cpu_num -W $word_size -q $mis_penalty -G $gap_cost -E $gap_extension -b $hits_return";
	if ($strand_specific)
	{
		$blast_param .= " -S 1";
	}
	system($blast_program." ".$blast_param);
	
	#��blast������ҵ�redundant��Ϣ
	my $result = findRedundancy($inset, $query);
#   if($$query =~ />NOVEL1\n/){#������
#	    print STDERR $result."good\n";
#		die "this is >NOVEL1";
#	}
	return $result;
}

sub findRedundancy
{
	my ($inset, $query) = @_;	
	my ($query_name, $query_length, $hit_name, $hit_length, $hsp_length, $identity, $evalue, $score, $strand, $query_start, $query_end, $hit_start, $hit_end, $query_to_end, $hit_to_end);
	my %hsp = ();  #�洢��ȡһ��query������hsp�����һ����
	my $hsp_count=0; 
	my $query_sequenc; 
	my $subject_sequenc; 	
	my $is_hsp = 1;

	my $bfh = IO::File->new("tem.paired") || "Can not open blast result file: tem.paired $!\n";
	while(<$bfh>)
	{
		chomp;
		if (/Query=\s(\S+)/ || eof)#һ������(��)Query Name�����ļ���ĩ����ֹֻ��һ��Query���������ǰ��һ��Query���е�hsp
		{
			#########################################
			#           ��¼ǰһ��hsp��¼           #
			#########################################
			if ($hsp_count > 0)#������ǵ�һ�γ��֣����Ѵ���hsp�����򱣴棨Ҳ���������ǰ��һ��hsp���
			{
				my $hsp_info =  $query_name."\t".$query_length."\t".$hit_name."\t".$hit_length."\t".
						$hsp_length."\t".$identity."\t".$evalue."\t".$score."\t".$strand."\t".
						$query_start."\t".$query_end."\t".$hit_start."\t".$hit_end;
				$hsp{$hsp_count} = $hsp_info;
			}

			#########################################
			#  ����������һ��query-hit������hsp	#
			#########################################
			if (scalar(keys(%hsp)) > 0)
			{
				foreach my $hsp_num (sort {$a<=>$b} keys %hsp)
				{
					my @one_hit = split(/\t/, $hsp{$hsp_num});
					unless(scalar(@one_hit) == 13) { die "Error! the parsed blast may lost some info:\n $hsp{$hsp_num} \n"; }
					
					my ($query_name, $query_length, $hit_name, $hit_length, $hsp_length, $identity, $evalue, $score, $strand,
					    $query_start, $query_end, $hit_start, $hit_end) =
					   ($one_hit[0], $one_hit[1], $one_hit[2], $one_hit[3], $one_hit[4], $one_hit[5], $one_hit[6], $one_hit[7],
					    $one_hit[8], $one_hit[9], $one_hit[10],$one_hit[11],$one_hit[12]);

					my $query_to_end = $query_length - $query_end;
					my $hit_to_end = $hit_length - $hit_end;
					my $hit_to_start = $hit_length - $hit_start;

					$identity =~ s/%//;
					if($hsp_length <= 40) {$min_identity = 96;}
					elsif($hsp_length > 40 && $hsp_length <= 50) {$min_identity = 97;}
					elsif($hsp_length > 50 && $hsp_length <= 100) {$min_identity = 98;}
					else{$min_identity = 99;}
                    
					#�����жϵ�˳��ǳ���Ҫ
					if ($identity < $min_identity)#hsp��identity���������ܺϲ��������ж�Ϊ�����࣬��Ҫ����һ��
					{
						next;#������֤�ˣ�����������Сidentity����ȥ�ж��������������������ȥ��
					}
					if ($query_start -1 <= $max_end_clip  && $query_to_end <= $max_end_clip)#query��subject�������ж�����
					{

						return $hit_name."\tr";
					}
					if($hsp_length >= $min_overlap)#���������(��overlap)��Ҫ�ϲ�
					{

					    my $combined_seq;
					    (my $query_seq = $$query) =~ s/^>[^\n]+\n//;#�滻��ǰ���query name
					        $query_seq =~ s/\s//g;
						my $hit_seq = $inset->{$hit_name}; 
						if ($strand==1)
						{
						        #������query��ǰ�����
							if($query_start -1 > $max_end_clip  && $query_to_end <= $max_end_clip && $hit_start <= $max_end_clip)
							{ 
						           my $query_string = substr($query_seq, 0, $query_start); 							   
							   my $hit_string = substr($hit_seq, $hit_start, $hit_to_start);
							   #print "type2\t".$hit_name."\t".$query_name."\t".$query_string."\t".$hit_string."\n";#����ϲ���Ϣ���˹�У�ԣ�������
							   $combined_seq = $hit_name.":".$query_name."\t".$query_string.$hit_string;
							   return $combined_seq;
							} 
							#������hit��ǰ�����
							if($query_start -1 <= $max_end_clip  && $query_to_end > $max_end_clip && $hit_to_end <= $max_end_clip)
							{ 
							   my $hit_string = substr($hit_seq, 0, $hit_end); 
							   my $query_string = substr($query_seq, $query_end, $query_to_end); 
							   #print "type3\t".$hit_name."\t".$query_name."\t".$hit_string."\t".$query_string."\n";#����ϲ���Ϣ���˹�У�ԣ�������
							   $combined_seq = $hit_name.":".$query_name."\t".$hit_string.$query_string;
							   return $combined_seq;
							} 
						}
						if ($strand==-1)
						{
							#������query��ǰ�����
							if($query_start -1 > $max_end_clip  && $query_to_end <= $max_end_clip && $hit_to_end <= $max_end_clip)
							{ 
							   my $query_string = substr($query_seq, 0, $query_start+1); 							   
							   my $hit_string = substr($hit_seq, 0, $hit_end-1);
							   rcSeq(\$hit_string, 'rc'); #�����еķ��򻥲�
							   #print "type4\t".$hit_name."\t".$query_name."\t".$query_string."\t".$hit_string."\n";#����ϲ���Ϣ���˹�У�ԣ�������
							   $combined_seq = $hit_name.":".$query_name."\t".$query_string.$hit_string;
							   return $combined_seq;
							} 
							#������hit��ǰ�����
							if($query_start -1 <= $max_end_clip  && $query_to_end > $max_end_clip && $hit_start-1 <= $max_end_clip)
							{ 
							   my $hit_string = substr($hit_seq, $hit_start, $hit_to_start); 
							   rcSeq(\$hit_string, 'rc'); #�����еķ��򻥲�
							   my $query_string = substr($query_seq, $query_end, $query_to_end); 
							   #print "type5\t".$hit_name."\t".$query_name."\t".$hit_string."\t".$query_string."\n";#����ϲ���Ϣ���˹�У�ԣ�������
							   $combined_seq = $hit_name.":".$query_name."\t".$hit_string.$query_string;
							   return $combined_seq;
							} 
						}
					}#�����������Ҫ�ϲ�
				}#ѭ��������û����������������жϣ��ж�Ϊ������
				return $hit_name."\tn";
			}
			
			#####################################
			#  ��ʼ��¼һ���µ�query	    #
			#####################################
			%hsp = ();$hsp_count = 0;
			$query_name = $1; $query_length = ""; $hit_name = ""; $hit_length = "";
		}
		elsif (/\s+\((\S+)\sletters\)/)#ȡQuery Length
		{
			$query_length = $1;
			$query_length =~ s/,//g;
		}
		
		elsif (/>(\S+)/)#һ������Hit Name
		{
			#########################################
			#           ��¼ǰһ��hsp��¼           #
			#########################################
			if ($hsp_count > 0 || eof)#���������Query���һ�γ��֣����Ѵ���hsp�����򱣴�ǰ��һ��hsp�����Ҳ���������
			{
				my $hsp_info =  $query_name."\t".$query_length."\t".$hit_name."\t".$hit_length."\t".
						$hsp_length."\t".$identity."\t".$evalue."\t".$score."\t".$strand."\t".
						$query_start."\t".$query_end."\t".$hit_start."\t".$hit_end;
				$hsp{$hsp_count} = $hsp_info;	
				$is_hsp = 0;
			}
			#################################
			#  ��ʼ��¼һ���µ�hit	        #
			#################################
		    $hit_name = $1; $hit_length = "";
		}
		elsif (/\s+Length = (\d+)/)
		{
				$hit_length = $1;
				$hit_length =~ s/,//g;
		}

		elsif (/Score =\s+(\S+) bits.+Expect(\(\d+\))? = (\S+)/)#һ������hsp
		{
			if ($hsp_count > 0 && $is_hsp == 1)
			{	
				my $hsp_info =  $query_name."\t".$query_length."\t".$hit_name."\t".$hit_length."\t".
								$hsp_length."\t".$identity."\t".$evalue."\t".$score."\t".$strand."\t".
								$query_start."\t".$query_end."\t".$hit_start."\t".$hit_end;
				$hsp{$hsp_count} = $hsp_info;
			}

			#################################
			#  ��ʼ��¼һ���µ�hsp		#
			#################################
			$is_hsp = 1;
			$hsp_count++; 
			$score = $1; $evalue = $3;
			$evalue = "1$evalue" if ($evalue =~ m/^e/);
			$query_start = 0; $query_end = 0; $hit_start = 0; $hit_end = 0;

		}
		elsif (/\s+Identities = (\d+)\/(\d+)\s+\((\S+)\)/ && $hsp_count >= 1)
		{
			$identity = $1/$2*100;
			$identity = sprintf("%.".(2)."f", $identity);
			if ( $1 > $2 ) { $hsp_length = $1; } else { $hsp_length = $2; }
		}

		elsif (/\s+Strand = (\S+) \/ (\S+)/ && $hsp_count >= 1) 
		{
			if ( $2 eq "Plus" ) { $strand = 1;} else { $strand = -1;}
		}

		elsif (/Query\:\s(\d+)\s+(\S+)\s(\d+)/ && $hsp_count >= 1)
		{
			if ($query_start == 0) { $query_start = $1; $query_start =~ s/,//g;}
			$query_end = $3;
			$query_end =~ s/,//g;
		}

		elsif (/Sbjct\:\s(\d+)\s+(\S+)\s(\d+)/ && $hsp_count >= 1) 
		{
			if ( $strand == -1 )#��Զ��֤$hit_start>=$hit_end
			{
				if ($hit_end == 0) { $hit_end = $1; $hit_end =~ s/,//g;  };
				$hit_start = $3;
				$hit_start =~ s/,//g;
			}
			else
			{
				if ($hit_start == 0) { $hit_start = $1; $hit_start =~ s/,//g; };
				$hit_end = $3;
				$hit_end =~ s/,//g;
			}
		}
		else
		{
			next;
		}
	}
	$bfh->close;
	return "null\tn";
}

sub rcSeq {
        my $seq_r = shift;
        my $tag = shift; defined $tag or $tag = 'rc'; # $tag = lc($tag);
        my ($Is_r, $Is_c) = (0)x2;
        $tag =~ /r/i and $Is_r = 1;
        $tag =~ /c/i and $Is_c = 1;
        #$tag eq 'rc' and ( ($Is_r,$Is_c) = (1)x2 );
        #$tag eq 'r' and $Is_r = 1;
        #$tag eq 'c' and $Is_c = 1;
        !$Is_r and !$Is_c and die "Wrong Input for function rcSeq! $!\n";
        $Is_r and $$seq_r = reverse ($$seq_r);
        # $Is_c and $$seq_r =~ tr/acgturyksbdhvnACGTURYKSBDHVN/tgcaayrmwvhdbnTGCAAYRMWVHDBN/;  # 2007-07-18 refer to NCBI;
        $Is_c and $$seq_r =~ tr/acgturykmbvdhACGTURYKMBVDH/tgcaayrmkvbhdTGCAAYRMKVBHD/; # edit on 2010-11-14;
        return 0;
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