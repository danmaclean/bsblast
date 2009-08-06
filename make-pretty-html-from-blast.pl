#!/usr/bin/perl


use strict ;
use warnings ;
use Bio::SearchIO;
use Getopt::Long;
use Data::Dumper;

my $usage = "usage: $0 \n-b blastfile -d directory to search \n

-b provide a single blast file
-d provide a directory with blast files in, all files ending in .blast will be used
-l limit on number of proteins to look at (default = 30)
\n\n";


my $blast_file;
my $directory;
my $limit;
GetOptions(


 'b=s' => \$blast_file,	  
 'd=s'     => \$directory,
  'l-i'    => \$limit

) ;

die $usage unless $blast_file or $directory;
$limit = 5 unless defined $limit;

my @files;

if ($blast_file){

	@files = ($blast_file);

}
else{
	opendir(DIR, $directory) or die "can't opendir $directory: $!";
	while (defined(my $file = readdir(DIR))) {
		if ( $file =~ m/(.+)\.blast$/ ) {
			push @files, $file;
		}
	}
	die "No files found\n\n" unless $files[0];
}




my %data;

foreach my $file (@files){
    warn "Parsing $file: \n" ;
    
    my $parser = new Bio::SearchIO(-format => 'blast', -file => $file) ; 
    while (my $result = $parser->next_result) {
	my $query = $result->query_name;
      my $num_hits = $result->num_hits;

      while (my $hit = $result->next_hit) {

	my $desc = $hit->description();
	while(my $hsp = $hit->next_hsp) {

          my $escore = $hsp->evalue;
	#my $start = $hsp->start('hit') ;
	#my $end = $hsp->end('hit') ;
	#my $strand = $hsp->strand('hit') ;
	#my $query_string = $hsp->query_string ;
	#my $hit_string = $hsp->hit_string ;
          my $frac_identical = $hsp->frac_identical;
          my $query_length = $result->query_length;
          my $hit_length = $hit->length;
          my $hsp_length = $hsp->hsp_length;
          my $match_over_percentage = ($hsp_length / $query_length) * 100; ## percentage hsp over length of query ie denovo peptide
          my $match_over_percentage2 = ($hsp_length / $hit_length) * 100;
	#my $query = $result->query_name;
          my $hit_name = $hit->name;
	  my $accession = $hit->accession;
          my $description = $result->query_description;
	  my $hit_desc = $hit->description;
          my @hit_range =  $hsp->range('hit');
          my $hit_range = join('-',@hit_range);
          my @query_range = $hsp->range('query');
          my $query_range = join('-',@query_range);
          my $homology = $hsp->homology_string;
	  $homology =~ s/\s/0/g;
          my $query_string = $hsp->query_string;


	  $data{$hit_name}{'length'} = $hit_length;
	  $data{$hit_name}{'unique_peptides'}{$query} = 1;
	  $data{$hit_name}{'description'} = $hit_desc;	  
	  $data{$hit_name}{'accession'} = $accession;
	  my $id = scalar(keys %{$data{$hit_name}}) -2 + 1;
          $data{$hit_name}{$id}{'peptide'} = $query;
          $data{$hit_name}{$id}{'escore'} = $escore;
          $data{$hit_name}{$id}{'identities'} = int($frac_identical * 100);
          $data{$hit_name}{$id}{'match_over_percentage'} = int($match_over_percentage);
          $data{$hit_name}{$id}{'range'} = $hit_range;
          $data{$hit_name}{$id}{'homology'} = $homology;
          $data{$hit_name}{$id}{'query_string'} = $query_string;



      }
    }
  }
}
my $os = $^O;

my $dir = $$ . '_blast_result/';
$dir = $$. '_blast_result\\' if $os =~ m/mswin32/i;

warn "Parsing done ... assessing peptides..\n";
die "results directory exists .. aborting\n" if -d $dir; 
mkdir $dir || die "cannot make results directory $dir\n";
chdir $dir || die "cannot change to results directory $dir \n";

my $index = "../$$.index.html";
$index = "..\\$$.index.html" if $os =~ m/mswin32/i;


open INDEX, ">$index";

my @ordered_proteins = sort {keys %{$data{$b}{'unique_peptides'}}   <=> keys %{$data{$a}{'unique_peptides'}}   || keys %{$data{$b}} <=> keys %{$data{$a}} } keys %data ;
warn "peptides sorted\n\n\n";
my $s = start_html();
print INDEX $s;

warn "peptide similarity matrix computed\n\n";

my @prots = @ordered_proteins;
my $done = '0';

my %data2 = %data;
my %need_html_for;
#foreach my $p (@prots){
while (my $p = shift @prots){
	
	 ## limit the result output to the top 30... for now
	next unless exists $data{$p};
	warn $done, " of 30 best proteins\n\n";
	#warn $p, "\n";
	last if $done >= 30;	
	my %peptide_similarity = %{make_matrix_of_similarity(\%data, $p)};
	$need_html_for{$p} = 1;
	## get the 5 most similar proteins...
	my $html = $data2{$p}{'accession'} . '.html';
	#my %pr = %{$$peptide_similarity{$p}};
	
	#print Dumper %peptide_similarity;
	#die;
	my @most_similar = sort { $peptide_similarity{$b} <=> $peptide_similarity{$a} } keys %peptide_similarity;
	#print Dumper @most_similar;
	
	my @slice = ($most_similar[0],$most_similar[1], $most_similar[2], $most_similar[3], $most_similar[4]);
	#print Dumper "@slice";
	
	my %tmp;
	foreach my $s (@slice){
		next unless defined $s;
		$need_html_for{$s} = 1;
		$tmp{$s} = 1;
		delete $data{$s};
	}
	############################

	##remove the elements from the list
	#for (my $i = 0; $i <= scalar(@prots); $i++){
	#	splice @prots, $i, 1 if exists $tmp{$prots[$i]}; 
		

	#}

	#foreach my $s (@slice){
	#	next unless defined $s;
	#	@prots = grep { $_ ne $s } @prots;
	#	$need_html_for{$s} = 1;
	#}

	print INDEX "<div>\n";
	print INDEX "<h1>$p</h1>\n";
	print INDEX "<h3>$data{$p}{'description'}</h3>\n";

	my $u = scalar(keys %{$data{$p}{'unique_peptides'}});
	print INDEX "<h1>$u &nbsp; unique peptides</h1>\n";
	print INDEX "<a href=\"$dir$html\" target =\"_blank\">view results</a>\n";
	print INDEX "<h3>Top (5 max) proteins with most similar patterns of peptide hits</h3>\n\n";
	foreach my $s (@slice){
		my $link = $data2{$s}{'accession'}. '.html';
		print INDEX "<a href=\"$dir$link\" target=\"_blank\">$s</a>\n";
	

	}
	print INDEX "</div>";	
	delete $data{$p};
	
 	++$done;
}


foreach my $protein (keys %need_html_for){

	my %tracks;
	my %done_peps;
	my %pep_counts;
	$tracks{'1'} = ('-' x $data2{$protein}{'length'}); 
	foreach my $id (keys %{$data2{$protein}} ){
		next if $id =~ m/length/;
		next if $id =~ m/unique_peptides/;
		next if $id =~m/description/;
		next if $id =~m/accession/;
		#print Dumper %{$data2{$protein}};
		$pep_counts{$data2{$protein}{$id}{'peptide'} }++; 
		unless (defined $done_peps{$data2{$protein}{$id}{'peptide'}}){ #warn "IN HERE!!!\n";
			$done_peps{$data2{$protein}{$id}{'peptide'}} = 1;
			my @range = split(/-/, $data2{$protein}{$id}{'range'});
			my @num_tracks = sort {$b <=> $a} keys %tracks;
			my $max_track = shift @num_tracks; #warn "MAX TRACK $max_track\n";
			my @aas = split(//,$data2{$protein}{$id}{'homology'});
			#for (my $i = 2; $i<= $max_track + 1; $i++){
			#	my $free = 1;
			#	for (my $j = $range[0]; $j<=$range[1]; ++$j){ 	
			#		$free = 0 if defined $tracks{$i}{$j};
			#	}
			#	if ($free){
			#		for (my $j = $range[0]; $j<=$range[1]; ++$j){ 
			#			my $aa = $j-$range[0]; 	
			#			$tracks{$i}{$j} = $aas[$aa];
			#		}	
			#	}
				
				
			#}
			my $i = 2;
			my $free = 0;
			until ($free){
				$free = 1;
				for (my $j = $range[0]; $j<=$range[1]; ++$j){ 	
					$free = 0 if defined $tracks{$i}{$j};
				}				
				if ($free){
				for (my $j = $range[0]; $j<=$range[1]; ++$j){ 
						my $aa = $j-$range[0]; 	
						$tracks{$i}{$j} = $aas[$aa];
					}	
				}
				else{
					$i++; 
				}
			} 	
		}	 
	}
	#print Dumper %tracks;
	my $html = $data2{$protein}{'accession'} . '.html';
	open HTML, ">$html";
	print HTML $s;
#	print INDEX "<div>\n";
#	print INDEX "<h1>$protein</h1>\n";
#	print INDEX "<h3>$data2{$protein}{'description'}</h3>\n";

	my $u = scalar(keys %{$data2{$protein}{'unique_peptides'}});
#	print INDEX "<h1>$u &nbsp; unique peptides</h1>\n";
#	print INDEX "<a href=\"$html\" target =\"_blank\">view results</a>\n";
#	print INDEX "</div>";	


	print HTML "<div>\n";
	print HTML "<h1>$protein</h1>\n";
	print HTML "<h3>$data2{$protein}{'description'}</h3>\n";

	print HTML "<h2>$u &nbsp; unique peptides</h2>\n";
	#my @sorted = sort {$pep_counts{$b} <=> $pep_counts{$a} } keys %pep_counts;
	my %done;
	foreach my $id (keys %{$data2{$protein}} ){ #warn $id, "\n";
		next if $id =~ m/length/;
		next if $id =~ m/unique_peptides/;
		next if $id =~ m/description/;
		next if $id =~m/accession/;
		next if exists $done{$data2{$protein}{$id}{'peptide'}};
		$done{$data2{$protein}{$id}{'peptide'}} = 1;
		print HTML "<p>$data2{$protein}{$id}{'peptide'} &nbsp; Occurred $pep_counts{$data2{$protein}{$id}{'peptide'}} E = $data2{$protein}{$id}{'escore'} &nbsp; Identity = $data2{$protein}{$id}{'identities'} &nbsp; Percent = $data2{$protein}{$id}{'match_over_percentage'}</p>\n";

	}
	
	print HTML "<p>$tracks{'1'}</p>\n";

	for (my $i = 2; $i <= scalar(keys %tracks); $i++){
		#warn "in here\n";
		my $track_string = '';
		for (my $j = 1; $j< $data2{$protein}{'length'}; $j++ ){

			if (defined $tracks{$i}{$j}){
				$track_string .= $tracks{$i}{$j};
			}
			else {

				$track_string .= '-';

			}
		}
		print HTML "<p>$track_string</p>\n\n\n";
	}
	print HTML "</div>\n";
	my $e = end_html();
	print HTML $e;
	close HTML;
}


my $e = end_html();
print INDEX $e;
close INDEX;
###########################################

sub make_matrix_of_similarity{
	my $data = shift;
	my $p = shift;
	my %output;
	warn $p, "\n"; 
	my %list1 = %{$data{$p}{'unique_peptides'}};
	foreach my $protein (keys %{$data} ){
	#	warn "!$p!", "\t", "!$protein!", "\n";
		#foreach my $protein2 (keys %{$data}){
			if ($p ne $protein ){
			#warn $p, "\t", $protein, "\n";

			my %list2 = %{$data{$protein}{'unique_peptides'}};

			my $similarity = calc_sim(\%list1, \%list2);
			#warn "sim = $similarity\n";
			$output{$protein} = $similarity if $similarity > 0;
			}
		#}
	}
	#print Dumper %output;
	#die;
	return \%output;



}

sub calc_sim{
	my $list1 = shift;
	my $list2 = shift;

	if ( scalar (keys %$list1) > scalar(keys %{$list2}) ){
		my $i = 0;
		foreach my $k (keys %{$list1} ){

			if (exists $$list2{$k}){

				++$i;

			}

		}
		return  ( ($i / scalar (keys %{$list1}) ) * 100);

	}
	else {

		my $i = 0;
		foreach my $k (keys %{$list2} ){

			if (exists $$list1{$k}){

				++$i;

			}

		}
		return  ( ($i / scalar (keys %{$list2}) ) * 100);
	

	}


}


sub start_html{


return '


<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
  "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">

<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
  <head>
    <meta http-equiv="content-type" content="text/html;charset=UTF-8" />
        <title>de Novo peptide blast</title>


  </head>

  

  <body>', "\n\n";

}

sub end_html{


return "</body></html>\n\n";

}
