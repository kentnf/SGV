#!/usr/bin/perl
use strict;
use Bio::Graphics;
use Cwd;

if (@ARGV < 4)
{
  print "usage: plot_results.pl input1 input2 output_dir contig_type\n";
  exit(0);
}

our $input1 = $ARGV[0]; # ��һ��ѭ����������ļ�����ѡ���������html�ļ�(�ܵı����ļ�)������ָ�������ļ�������
our $input2 = $ARGV[1]; # ����ļ���������ÿ���ο������Լ���html�ļ��������ӻ���ͼ�Σ�
our $output = $ARGV[2]; # 
our $type = $ARGV[3]; #
our $WORKING_DIR=cwd();#����Ŀ¼���ǵ�ǰĿ¼
our $OUTPUT_DIR =  "$output/$type"."_references"; # ����ڹ���Ŀ¼�����Ŀ¼��ע�����Ŀ¼������ǰ����������������
our $output_file = "$output/$type.html"; # ����ڹ���Ŀ¼�����Ŀ¼�������ļ�

my %all_hits;
my %out;
my %index;
my $line_number=0;
open(IN1, "$input1");
open(OUT, ">$output_file" );#���һ�����
my $out_table = qq'
<style type="text/css">
td,th
{
  font-family : Arial, Sans-Serif, Helv, Helvetica, Verdana;
  font-size: 12px;
}
</style>
<table border=1 cellpadding=3 style=border-collapse:collapse bordercolor=#gray width=780 align=center>
  <tr bgcolor=#e2e8ec>
    <th>Reference</th>
    <th>Length</th>
    <th>Coverage(%)</th>
    <th>#contig</th>
    <th>Depth</th>
    <th>Genus</th>
    <th width="50%">Description</th>
  </tr>';

#�����������ĸ�ʽ�������html�ļ�
while (<IN1>) {
	chomp;
	$line_number++;	
	my @each_line = split(/\t/, $_); 
	$all_hits{$each_line[0]}=1; 
	$index{$each_line[0]}=$line_number;
	my $link="$type"."_references/".$each_line[0].".html";
	my $coverage= sprintf("%.3f", $each_line[3])*100;
	my $contigs= $each_line[5];
	my $depth= sprintf("%.1f", $each_line[6]);
	$out_table .= qq' 
   <tr>
        <td><a href="$link">$each_line[0]</a></td>
        <td>$each_line[1]</td>
        <td>$each_line[2]($coverage)</td>
        <td>$contigs</td>	
        <td>$depth</td>
        <td>$each_line[7]</td>	
	<td width="50%">$each_line[8]</td>	
   </tr>';	
}
close(IN1);
$out_table.= qq'</table>';
print OUT "<br>$out_table<br>";
close(OUT);

open(IN2, "$input2");
while (<IN2>) {
	chomp; 
	my @each_line = split(/\t/, $_); 
	if(defined $all_hits{$each_line[2]}){
		push(@{$out{$each_line[2]}},$_);
	}
}
close(IN2);


foreach my $hit (sort { $index{$a} <=> $index{$b} } keys %out){
	my @lines = @{$out{$hit}};#��ȡ�õ��������ڱ���ѭ����hit��������
	my $out_hit = $OUTPUT_DIR."/$hit.html";#ÿ��hit�����Լ���folder	
	open(OUT1, ">$out_hit") or print "you forgot to make $OUTPUT_DIR";#ÿ��hit����Ҫ����һ���ļ����
	print OUT1 qq'<div align=center>';#�����ǩ����֤��ͼƬ����
	draw_img(\$hit, \@lines, \*OUT1);#�����ļ���������Ӻ����л�ͼ
	print OUT1 qq'</div>';
	
	#�õ�ÿ��contig�ȶԵĽ��������һ��html�Ĵ�
	my $table_info = output_aligments(\@lines);
	print OUT1 "<br>$table_info<br>";
	close(OUT1);

	
=head	#�����������������
	foreach my $each_line (@lines){
		print $each_line."\n";	
	}
=cut
}

#################
##    �ӳ���   ##
#################
sub draw_img
{
	my ($hit_id, $align_info, $out) = @_;
	my @cols = split(/\t/, @$align_info[0]);#�����ڱ���hit�ĵ�һ�������ó���Ϊ����ȡhit�ĳ�����Ϣ
	my $full_length = $cols[3];

	my $panel = Bio::Graphics::Panel->new(
                                      -length => $full_length,
                                      -width  => 610
                                      -pad_left => 25,
                                      -pad_right => 25,
                                     );
        #��Ϊ�ο����еģ�������һ��feature ����
	my $full_length_seq = Bio::Graphics::Feature->new(
                                                -start => 1,
                                                -end => $full_length,
						-display_name => "$$hit_id",
                                                );
        #��������feature�����һ����ͷ��һ����track
	$panel->add_track( 
                        $full_length_seq,
                        -glyph => 'arrow',
                        -tick => 2,
                        -fgcolor => 'black',
                        -double => 1,
                 );
	#��������feature�������һ�����ε�һ����track
	$panel->add_track(
			$full_length_seq,
			-glyph => 'generic',
			-bgcolor => 'blue',
			-label => 1,
			-link => "http://www.ncbi.nlm.nih.gov/nuccore/$$hit_id"
		);
	#�ٽ�һ��track����ʾ���е�����feature��contigs��
	my $track = $panel->add_track(
                              -glyph => 'graded_segments',
                              -label => 1,
			      -bgcolor => 'red',
			      -sort_order  => 'high_score',
			      -min_score   => 80,
                              -max_score   => 100,#��identity��Ϊscore
			      -key         => 'domains',
			      -link        => '#$id',
			      -description => sub {
                                my $feature = shift;
                                return unless $feature->has_tag('description');
                                my ($description) = $feature->each_tag_value('description');
                                "$description";
                               },

                             );
	#ÿ��ѭ�����һ��feature��������ʾ
	foreach my $each_contig (@$align_info)
	{
		my @mm = split(/\t/, $each_contig);
		my $feature = Bio::Graphics::Feature->new(
							-name	      => $mm[0],
							-score        => $mm[5],#��identity��Ϊscore
							-start        => $mm[11],
							-end          => $mm[12],
							-primary_id   => $mm[0],
							);

		$track->add_feature($feature);
		
	}
#���	
	my ($url,$map,$mapname) = $panel->image_and_map(
						-root => $OUTPUT_DIR,
						-url=> "$$hit_id",
						);
						
	print $out qq(<img src="$url" border=0 usemap="#$mapname">),"\n";
	print $out $map;
}

sub output_aligments
{
	my $align_info = shift;
	# produce domains alignment for one protein
	my $break_length = 95;
	
	#д���ı���
	my $out_table = qq'
<style type="text/css">
td,th
{
  font-family : Arial, Sans-Serif, Helv, Helvetica, Verdana;
  font-size: 12px;
}
</style>
<table border=1 cellpadding=3 style=border-collapse:collapse bordercolor=#gray width=680 align=center>
  <tr bgcolor=#e2e8ec>
    <th>Order</th>
    <th>Query ID</th>
    <th>Query Start</th>
    <th>Query End</th>
    <th>Subjct Start</th>
    <th>Subjct End</th>
    <th>Identity</th>
    <th>E value</th>
    <th>Strand</th>
  </tr>
';

	my $k = 0;

	foreach my $each_contig (@$align_info)
	{
		$k++;

		my @cols = split(/\t/, $each_contig);#�õ�һ�����ݵ�������
		#print $each_contig."\n";		
		my @aligned_query = split(/,/, $cols[14]);#�õ�query������
		my @aligned_hit = split(/,/, $cols[15]);#�õ�hit������
		my @aligned_string = split(/,/, $cols[16]);#�õ�aligment������
		my $i=0;
		my $align;
		foreach my $each_seg (@aligned_query)
		{
			#print $each_seg."\n";			
			$align.="<span style=\"background-color:#99FFFF\";>$each_seg</span><br>";
			$align.="<span style=\"background-color:#99FFFF\";>$aligned_string[$i]</span><br>";
			$align.="<span style=\"background-color:#99FFFF\";>$aligned_hit[$i]</span><br> <br>";
			$i++;
		}

		$out_table .= qq' 
   <tr>
	<td>$k</td>
        <td><a name=$cols[0]>$cols[0]</a></td>
        <td>$cols[9]</td>
        <td>$cols[10]</td>
        <td>$cols[11]</td>
        <td>$cols[12]</td>
        <td>$cols[13]</td>
	<td>$cols[6]</td>
	<td>$cols[8]</td>
   </tr>
   <tr>
        <td colspan=9>
	<div style="padding:5px;">
        <b>Alignment:   </b><br>
	<span style="font-size:12px; color:#000; font-face: courier">
<pre>
$align
</pre>
	</span>
        </td>
	</div>
   </tr>';
	}

	$out_table.= qq'</table>';

	return $out_table;
}
