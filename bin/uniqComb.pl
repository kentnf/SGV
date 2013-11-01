#!/usr/bin/perl

# 2006-12-1 14:21 
# 2008/-5/08 
# 2010/05/06


use strict;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"help!","index:s","col:s","newCol:s","exist!");

my $info = <<INFO;
perl $0 srchFile -index indxFile -col indxCol -newCol srchCol [-exist]
INFO

#srchFile����Ҫ��ȡ���ݵ�file
#-index����Ϊindex��file
#-col��index file���ļ��������Ƚϣ���","�ָ�
#-newCol��srchFile���ļ��������Ƚϣ���","�ָ�
#-exist:�����ͬ�Ļ��ǲ�ͬ��

!@ARGV and -t and die "$info"; 
$opts{help} and die "$info";
my @Cols = (0); 
# my $col= 0;
# defined $opts{col} and $col = $opts{col};
defined $opts{col} and @Cols = split(/,/, $opts{col}); 
my @newCols = (0); 
#my $newCol = 0;
#defined $opts{newCol} and $newCol = $opts{newCol};
defined $opts{newCol} and @newCols = split(/,/, $opts{newCol}); 

my %key;
open(INDEX,'<',$opts{index}) or die "index?\n";
while (<INDEX>) {
	chomp;
	# s/\s+$//; # added 20100506 
	my @temp = split(/\t/,$_);
	my $tkey = join("\t", @temp[@Cols]); 
	$key{$tkey} = 1;
}
close INDEX;

while (<>) {
	chomp;
	# s/\s+$//; # added 
	my @temp = split(/\t/,$_);
	my $tkey = join("\t", @temp[@newCols]); 
	if ($opts{exist}) {
		defined $key{$tkey} and print STDOUT "$_\n";
	}else{
		defined $key{$tkey} or print STDOUT "$_\n";
	}
}




