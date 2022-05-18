#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use Bio::SeqIO; 
use Bio::AlignIO; 
use Bio::SimpleAlign;
use Bio::TreeIO;
use Bio::Index::Fasta; 

#my $referenceIDs = $ARGV[0]; 
my $windw        = $ARGV[0];
my $outDIR       = $ARGV[1];

my @scids = qw(BELUGA BOTTLE FINLES HARBOU INDOBO LIPOTE NARWHA ORCAWL PILOTW WHITES); 

my $ids = qx(grep ">" dataRM/Lipote.fa); 
$ids =~ s/>LIPOTE_//g; 
$ids =~ s/ //g; 
my @ids = split("\n", $ids); 
my $counter = 0;


opendir (DIR, "dataRM") or die $!;
my @files = grep /.*\.fa$/, readdir(DIR);  
close(DIR);
# and index them
my $inx = Bio::Index::Fasta->new(-filename => "fastaINDX_RM",  -write_flag => 1);
foreach my $infile (@files) {
     print "indexing $infile\n";
     $inx->make_index("dataRM/$infile");
}


foreach my $scaffold (@ids) {
    my $lipote =  "LIPOTE" . "_" . $scaffold; 
    my $seq = $inx->get_Seq_by_id($lipote); 
    my $length = $seq->length(); 
    
    for (my $i=0; $i<($length-$windw); $i+=($windw+1000000)) {        
        $counter++; 
        $counter = sprintf("%05d", $counter);
        print "$counter\t$i\t"; 
                
        my %missingHash = ();
        my @ns;
        foreach my $scID (@scids) {
            # print "$scID "; 
            my $name = $scID . "_" . $scaffold; 
            #print "    name $name\t"; 
            my $seq = $inx->get_Seq_by_id($name); 
            my $seqString  = $seq->seq;
            my $subStr = substr($seqString,$i,$windw);
            my $ns = ($subStr =~ tr/N/N/); 
            push @ns, $ns; 
            #print "   ns: $ns\t"; 
            
            my $noUS = $scID;
            $noUS =~ s/_//;
            $ns = $ns/$windw *100;
            # print "$scID ns(pc): $ns\t"; 
            #print "$ns\t"; 
            my $rounded = sprintf("%.0f", $ns);  
            #print "   noUS $noUS \t"; 
            print "$rounded\t"; 
            $missingHash{$noUS} = $rounded; 
        }
                
        @ns = sort {$a<=>$b} @ns; 
        my $max = $ns[-1] / $windw;
        
        my $fileName; 
        if ($max < 0.5) {
            chdir "/raid6/stefanie/mick/windowTrees05"; 
            
            my $redcdAln = 0; 
            $fileName = $counter . "_" .  $scaffold . ".phy";
            $fileName =~ s/_\./\./; 
            print "$fileName\t";
            
            my $newAln = Bio::SimpleAlign->new();
            my $out = Bio::AlignIO->new(-file => ">$outDIR/$fileName", -format => 'fasta');
            
            my $outID; 
            foreach my $scID (@scids) {
                my $SCname = $scID . "_" . $scaffold; 
                #print "$SCname\n"; 
                my $seq = $inx->get_Seq_by_id($SCname); 
                my $seqString  = $seq->seq;
                my $subStr = substr($seqString,$i,$windw);                
                my $newID = $scID;
                $newID =~ s/_//; 
                my $s1 = new Bio::LocatableSeq(-seq => "$subStr", -id  => "$newID");
                $newAln->add_seq($s1);
            }
                        
            $out->write_aln($newAln); 
            
            chdir "$outDIR";
            my $outFile = $fileName;
            $outFile =~ s/\.fasta//; 
            
        } 
        else {
            print "alow50PComm\n"; 
        }
        %missingHash = ();    
    } 
    
    
    
    
}
system "rm fastaINDX"; 


