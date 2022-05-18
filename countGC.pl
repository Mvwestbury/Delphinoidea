#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use Bio::SeqIO; 
use Bio::AlignIO; 
use Bio::SimpleAlign;
use Bio::TreeIO;
use Bio::Index::Fasta; 


my $alnsDir = "out50k"; 

opendir (DIR, "$alnsDir") or die $!;
my @files = grep /.*\.phy$/, readdir(DIR);  
close(DIR);

foreach my $file (@files) {
    my $at_w = 0; 
    my $gc_w = 0; 
    my $ns_w = 0;
    my $gp_w = 0; 
    
    
    print "$file\t"; 
    
    my $alignInput = new Bio::AlignIO (-file   => "$alnsDir/$file", -format => "phylip");
    my $align = $alignInput->next_aln;
    
    foreach my $seq ( $align->each_seq() ) {  
        my $seqString = $seq->seq; 
        
        my $at_s = ($seqString =~ tr/AT/AT/); 
        my $gc_s = ($seqString =~ tr/GC/GC/); 
        my $ns_s = ($seqString =~ tr/N/N/); 
        my $gp_s = ($seqString =~ tr/-/-/); 
        
        $at_w = $at_w + $at_s;
        $gc_w = $gc_w + $gc_s;
        $ns_w = $ns_w + $ns_s;
        $gp_w = $gp_w + $gp_s;
                
    }
    
    my $total = $at_w + $gc_w; 
    my $gc_pc = $gc_w / $total; 
    
    
    print "$at_w\t";  
    print "$gc_w\t";  
    print "$ns_w\t"; 
    print "$gp_w\t";   
    print "$total\t";
    print "$gc_pc\n";   
    
    
}
print "\n"; 

