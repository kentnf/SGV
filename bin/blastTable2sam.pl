#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;


#��blast�õ���paired����ļ�תΪsam�ļ�����
if (@ARGV ==0)
{
  print "usage: blast2sam.pl input > output\n";
  exit(0);
}

our $input = $ARGV[0]; 
my $hit_name;
my $hit_len;
my %name_len;
#��ͨ��һ��ѭ������ȡ���е�hit�����浽hash����
open(IN, "$input");
while (<IN>) {
	my @cols= split(/\t/,$_);
	$hit_name = $cols[2];
	$hit_len = $cols[3];
	defined $name_len{$hit_name} or $name_len{$hit_name} = $hit_len;	
}
close(IN);
#�����Щhit����Ϣ����Ϊsam�ļ���header
foreach my $key (keys(%name_len))
{
	print "\@SQ\tSN:".$key."\tLN:".$name_len{$key}."\n";#���hit���ƣ�ƽ�����
}

my ($qlen, $slen, $q, $s, $qbeg, $qend, @sam, @cigar, @cmaux);
@sam = (); @sam[4,6,7,8,10] = (255, '*', 0, 0, '*');#��Щ�ж��������ݶ��ǹ̶�ֵ
open(IN, "$input");
while (<IN>) {
	chomp;
	my @cols= split(/\t/,$_);
	$sam[0] = $cols[0]; #query������
	if($cols[8] eq "-1"){
		$sam[1] = 0x10;}
	else{$sam[1] = 0x00;}
	$sam[2] = $cols[2]; #subject������
	$qlen = $cols[1]; #query�ĳ���
	$qend = $cols[10]; #query�����һ�����
	$sam[3] = $cols[11]; #subject�Ŀ�ʼ�����Ҳ����query aligned��λ��
	@cigar = ();
	@cmaux = (0, 0);
	my @q_seqs= split(/,/,$cols[14]);#��������query�ķֶ�����
	my @s_seqs= split(/,/,$cols[15]);#��������query�ķֶ�����
	my $l = scalar(@q_seqs);#ȡ����size,���˶����һ��","����Ŀ��ַ�����������
	$sam[9] = '';#ÿ�α����ʼ��
	for (my $i = 0; $i < $l; ++$i)
	{
		$q_seqs[$i] =~ /Query\:\s(\d+)\s*(\S+)\s(\d+)/;
		$q = $2;
		my $x = $q;
		$x =~ s/-//g; 
		$sam[9] .= $x;#��������query����Ƭ��
		
		$s_seqs[$i] =~ /Sbjct\:\s(\d+)\s*(\S+)\s(\d+)/;
		$s = $2;
		&aln2cm(\@cigar, \$q, \$s, \@cmaux);		
	}	
	#��ȡ��������Ϣ
	my ($as, $ev) = (int($cols[7] + .499), $cols[6]);
	$ev = "1$ev" if ($ev =~ /^e/);
	@sam[11,12] = ("AS:i:$as", "EV:Z:$ev");
	&blast_print_sam(\@sam, \@cigar, \@cmaux, $qlen - $qend);
}
close(IN);


sub blast_print_sam {
  my ($sam, $cigar, $cmaux, $qrest) = @_;
  push(@$cigar, $cmaux->[1] . substr("MID", $cmaux->[0], 1));
  push(@$cigar, $qrest . 'H') if ($qrest);
  if ($sam->[1] & 0x10) {
	@$cigar = reverse(@$cigar);
	$sam->[9] = reverse($sam->[9]);
	$sam->[9] =~ tr/atgcrymkswATGCRYMKSW/tacgyrkmswTACGYRKMSW/;
  }
  $sam->[9] = '*' if (!$sam->[9]);
  $sam->[5] = join('', @$cigar);
  print join("\t", @$sam), "\n";
}

sub aln2cm {
  my ($cigar, $q, $s, $cmaux) = @_;
  my $l = length($$q);
  for (my $i = 0; $i < $l; ++$i) {
	my $op;
	# set $op
	if (substr($$q, $i, 1) eq '-') { $op = 2; }#���query�з���һ��'-'��Ӧ����"D"
	elsif (substr($$s, $i, 1) eq '-') { $op = 1; }#���subject�з���һ��'-'��Ӧ����"I"
	else { $op = 0; }#�������"M"
	# for CIGAR
	if ($cmaux->[0] == $op) {
	  ++$cmaux->[1];#��¼"MID"ÿ�����͵�����
	} else {
	  push(@$cigar, $cmaux->[1] . substr("MID", $cmaux->[0], 1));
	  $cmaux->[0] = $op; $cmaux->[1] = 1;
	}
  }
}
