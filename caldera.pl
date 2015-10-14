#!/usr/bin/perl -w   
# 
# caldera.pl
# Identify interesting Y chromosome genes from thousands of Y chromosome sequences 
# Brantley Hall, August 2012  
#
#		     	      
# Problem - CQ-calculate is very good at finding Y sequences, but after we have thouands of candidate Y sequenes it is very difficult to find realistic Y genes 
# Solution - integrate a large amount of information that will allow the systematic discovery of Y genes
# Currently this is set up specific for the species and data available in the Tu Lab 
#
# Input 
# 1. CQ-calculate output
# 2. Multi-fasta file containing sequences listed in CQ-calculate
#
# Output Fields 
# 1. Name 
# 2. Length 
# 3. Female bowtie alignments 
# 4. Male bowtie alignments 
# 5. CQ 
# 6. Female blast alignments 
# 7. Male blast alignments 
# 8. rCQ 
# 9. 0-2hr blast alignments 
# 10. 2-4hr blast alignments 
# 11. 4-8hr blast alignments 
# 12. 8-12hr blast alignments 
# 13. Larva blast alignments 
# 14. Pupa blast alignments 
# 15. Female blast alignments 
# 16. Male blast blast alignments 
# 17. Blastn genome alignments 
# 18. Tblastx genome alignments 
# 19. Blastx alignments 
# 20. TE Words 
# 21. Top 5 blastx alignments 

use Getopt::Long;
use File::Basename;
use Bio::Seq;
use Bio::SeqIO;

GetOptions( "i=s"	=> \$input,
            "cq=s"      => \$CQ_output,    
	    "ste"	=> \$stephensi,
	    "gam"       => \$gambiae, 
	    "aeg"       => \$aegypti,
	    "test"      => \$test,
	    "o=s"       => \$output,
	    "s=s"       => \$sequence,			 	
	    "h"		=> \$help,
	    "help"	=> \$help 	
	   );	

# Initial and database setup 	   

&setup; 

open (IN, "$CQ_output") or die "$!";

while (<IN>) { 

	if ($sequence == $w) { 

		@line = split ' ', $_; 
		$name = $line[0];
		$length = $line[1]; 
		$female_bt = $line[2];
		$male_bt = $line[3];
		$CQ = $line[4];
		        
		&retrieve; 
		
		if ($big) { 
		
		        &generic_blastn("$small_female_database");
		        $small_female_lines = $lines;
		
		        if ($small_female_lines < 25) { 
		        
		                &sensitive_blastn("$female_database");
		                $female_lines = $lines;
		                
		                &sensitive_blastn("$male_database"); 
		                $male_lines = $lines; 
		                 
		                &rCQ_calculate;   
		        }
		        
		        else { 
		        
		                $female_lines = "-";
		                $male_lines = "-";
		                $rCQ = "-"
		        }                         
		}
		
		else { 
		
		        &generic_blastn("$female_database");
		        $female_lines = $lines; 
		
		        if ($female_lines < 40) { 
		               
		                &sensitive_blastn("$male_database"); 
		                $male_lines = $lines; 
		                 
		                &rCQ_calculate;
		        }
		
		        else { 
		
		                $male_lines = "-";
		                $rCQ = "-"
		        }
		}                
		        
		&blastn_trx("$emb02hr"); 
		$emb02hr_lines = $lines; 
		     
		&blastn_trx("$emb24hr"); 
		$emb24hr_lines = $lines; 
		  
		&blastn_trx("$emb48hr"); 
		$emb48hr_lines = $lines; 
		        
		&blastn_trx("$emb812hr"); 
		$emb812hr_lines = $lines; 
		       
		&blastn_trx("$larva"); 
		$larva_lines = $lines; 
		        
		&blastn_trx("$pupa"); 
		$pupa_lines = $lines; 
		       
		&blastn_trx("$male_trx"); 
		$male_trx_lines = $lines; 
		       
		&blastn_trx("$female_trx"); 
		$female_trx_lines = $lines; 
		     
		&how_repetitive_s_blastn("$genome");
		$how_repetitive_s_blastn = $lines;

		&how_repetitive_r_blastn("$genome");
		$how_repetitive_r_blastn = $lines;
		
		# If a sequence does not match anything in the A. aegypti genome then blast it against the A. aegypti trace file 
		if ($aegypti && $how_repetitive_s_blastn == 0 && $how_repetitive_r_blastn < 10) { 
	
		        &blast_trace("$trace");
		        $trace_lines = $lines; 
		}       
		
		else { 
		
		        $trace_lines = "-";
		}
		
		&repeatdb("$repeats");
		$repeats_alignments = $lines;
	
		&ygenedb("$ygenes");
		$ygenes_alignments = $lines; 
	
		&knownysequences("$knowny"); 
		$knowny_alignments = $lines; 

		if ($emb02hr_lines || $female_trx_lines > 20) { 
		
		        $TE_words = 0;
                        $blastx_alignments = 0; 
                        $top5 = '';
		        
		}
		
		else { 
		
		        &blastx("$nr");
                }                 
		
		# On the first cycle print lables 
		if (! $printed) { 
		
		        print OUT "Name\t Length\t Female-g-bt\t Male-g-bt\t CQ\t Female-g-blastn\t Male-g-blastn\t rCQ\t Emb0-2hr_trx\t Emb2-4hr_trx\t Emb4-8hr_trx\t Emb8-12hr_trx\t Larva_trx\t Pupa_trx\t Female_trx\t Male_trx\t How-Repetitive-s-blastn\t How-Repetitive-r-blastn\t Perc-identity\t Trace-lines\t Repeats-alignments\t Y-gene-alignments\t Which-Y-gene\t Known-Y-alignments\t blastx-alignmnets\t TE-words\t Top-5-blastx-alignments\n";   
		        $printed = 1;
		}
		
		print OUT "$name\t $length\t $female_bt\t $male_bt\t $CQ\t $female_lines\t $male_lines\t $rCQ\t $emb02hr_lines\t $emb24hr_lines\t $emb48hr_lines\t $emb812hr_lines\t $larva_lines\t $pupa_lines\t $female_trx_lines\t $male_trx_lines\t $how_repetitive_s_blastn\t $how_repetitive_r_blastn\t $variable\t $trace_lines\t $repeats_alignments\t $ygenes_alignments\t $which_y_gene\t $knowny_alignments\t $blastx_alignments\t $TE_words\t $top5\n";
		
		$w += 1; 
		$sequence += 1;
	}

	else { 
        
                $w += 1; 
        }
}

sub setup { 
	
	# List of mandatory options when user does not specify all mandatory options
	$help_message = "\nYou seem to be missing a mandatory option.\n\nUSAGE: caldera.pl [options]\nMandatory Options:\n\t-i\t [multi-fasta file containing reference sequences]\n\t\-cq\t [CQ output file]\n\n";

	if ($help) { 

		print "$help_message\n"; 
		exit;
	}
	 
	# Check to see if mandatory options are defined by user 
	if (! $input) { 

		print "$help_message\n";
		exit; 
	}

	# Check to see if the input filq, male, and female databases specified by the user exist
	if (! -e $input) { 
	
		print "\nThe input file: $input does not exist.\n\n"; 
		exit; 
	} 
	
	# Check to see if mandatory options are defined by user 
	if (! $CQ_output) { 

		print "$help_message\n";
		exit; 
	}

	# Check to see if the input filq, male, and female databases specified by the user exist
	if (! -e $CQ_output) { 
	
		print "\nThe input file: $CQ_output does not exist.\n\n"; 
		exit; 
	} 
	
	# Check to see if the input filq, male, and female databases specified by the user exist
	if (! $output) { 
	
		print "\nPlease provide an output file.\n\n"; 
		exit; 
	} 

        # Open the ouput file 
        open (OUT, ">$output") or die "Can't open the output file: : $!\n";
	
	# Setup variables 
	$blastx_alignments = 0; 
        $TE_words = 0;
        $top5 = '';
        # For start in middle 
	$w = 1;
	$sequence = 1;
	
	# Setup database sources 
	if ($stephensi) { 
		
		print "stephensi\n";
		$male_database = "/a/stephensi/DNAseq/ucDavis_illumina_84pbMP/stephensi_male_illumina_ucdavis.fa";
		$female_database = "/a/stephensi/DNAseq/ucDavis_illumina_84pbMP/stephensi_female_illumina_ucdavis.fa";
		$emb02hr = "/a/stephensi/RNAseq/10_embryo_0-1h_mRNA_VCB.fasta";
		$emb24hr = "/a/stephensi/RNAseq/11_embryo_2-4h_mRNA_VCB.fasta";
		$emb48hr = "/a/stephensi/RNAseq/01_embryo_4-8h_mRNA_Iowa.fasta";
		$emb812hr = "/a/stephensi/RNAseq/02_embryo_8-12h_mRNA_Iowa.fasta";
		$larva = "/a/stephensi/RNAseq/03_larvae_mRNA_Iowa.fasta";
		$pupa = "/a/stephensi/RNAseq/04_pupae_mRNA_Iowa.fasta";
		$male_trx = "/a/stephensi/RNAseq/06_male_1-5d_mRNA_Iowa.fasta";
		$female_trx = "/a/stephensi/RNAseq/05_female_1-5d_mRNA_Iowa.fasta";
		$genome = "/a/stephensi/assembled_genome/scaffolds-stephensi-genome-010114.fa";
		$repeats = "/a/gambiae/repeatScout_output/AgamP3-scout";
		$ygenes = "/home/brantley/db/caldera/data/ygenedb";
	        $knowny = "/a/gambiae/genome_vectorbase/gamY";
	}
	
	if ($gambiae) { 
	
	        print "gambiae\n";
		$male_database = "/a/gambiae/DNAseq/iowa_Dec2011_illumina/maleillu.fa"; 
		$female_database = "/a/gambiae/DNAseq/iowa_Dec2011_illumina/femaleillu.fa";
		$emb02hr = "/a/gambiae/RNAseq/0-2hrEmb.fa";
		$emb24hr = "/a/gambiae/RNAseq/2-4hrEmb.fa";
		$emb48hr = "/a/gambiae/RNAseq/4-8hrEmb.fa";
		$emb812hr = "/a/gambiae/RNAseq/8-12hrEmb.fa";
		$larva = "/a/gambiae/RNAseq/larva.fa";
		$pupa = "/a/gambiae/RNAseq/pupa.fa";
		$male_trx = "/a/gambiae/RNAseq/male.fa";
		$female_trx = "/a/gambiae/RNAseq/female.fa";
		$genome = "/a/gambiae/genome_vectorbase/chromosomes-PEST.AgamP3.fa";
		$repeats = "/a/gambiae/repeatScout_output/AgamP3-scout";
		$ygenes = "/home/brantley/db/caldera/data/ygenedb";
	        $knowny = "/a/gambiae/genome_vectorbase/gamY";  
	}
	
	if ($aegypti) { 
	
			# Test to speed things up because the Aedes aegypti is data is so big 
	        $big = 1; 
	        
	
	        print "aegypti\n"; 
	        $male_database = "/b/aegypti/DNAseq/vbi-male-female/fasta/male.fa";
		$small_female_database = "/b/aegypti/DNAseq/vbi-male-female/fasta/1r1/female-1r1.fa";
		$female_database = "/b/aegypti/DNAseq/vbi-male-female/fasta/female.fa";
		$emb02hr = "/b/aegypti/RNAseq/01_aegypti_mRNA_0-2H_emb_33bp.fasta";
		$emb24hr = "/b/aegypti/RNAseq/02_aegypti_mRNA_2-4H_emb_33bp.fasta";
		$emb48hr = "/b/aegypti/RNAseq/03_aegypti_mRNA_4-8H_emb_33bp.fasta";
		$emb812hr = "/b/aegypti/RNAseq/04_aegypti_mRNA_8-12H_emb_33bp.fasta";
		$larva = "/b/aegypti/RNAseq/05_aegypti_mRNA_larvae_39bp.fasta";
		$pupa = "/b/aegypti/RNAseq/06_aegypti_mRNA_pupae_41bp.fasta";
		$male_trx = "/b/aegypti/RNAseq/07_aegypti_mRNA_1-5D_male_38bp.fasta";
		$female_trx = "/b/aegypti/RNAseq/08_aegypti_mRNA_0-1D_ovary_39bp.fasta";
		$genome = "/b/aegypti/genome_vectorbase/contigs-Liverpool.AeagL1.fa";
		$repeats = "/b/aegypti/tefam-repeats/tefam1091";
		$ygenes = "/home/brantley/db/eldorado/candidate_gene_db/candidate_genes";
	        $knowny = "/a/gambiae/genome_vectorbase/gamY";
	        $trace = "/b/aegypti/trace_NCBI/aegypti_trace.fa";
		
	}	
	        
	if ($test) { 
	
	        print "Test Dataset\n\n";
	        $male_database = "/home/brantley/Dropbox/caldera/mtest"; 
		$female_database = "/home/brantley/Dropbox/caldera/ftest";
		$emb02hr = "/home/brantley/Dropbox/caldera/mtest"; 
		$emb24hr = "/home/brantley/Dropbox/caldera/ftest";  
		$emb48hr = "/home/brantley/Dropbox/caldera/mtest";  
		$emb812hr = "/home/brantley/Dropbox/caldera/ftest"; 
		$larva = "/home/brantley/Dropbox/caldera/mtest";  
		$pupa = "/home/brantley/Dropbox/caldera/ftest";  
		$male_trx = "/home/brantley/Dropbox/caldera/mtest";  
		$female_trx = "/home/brantley/Dropbox/caldera/ftest"; 
		$genome = "/home/brantley/Dropbox/caldera/mtest"; 
		$repeats = "/a/gambiae/repeatScout_output/AgamP3-scout"; 
		$ygenes = "/home/brantley/db/caldera/data/ygenedb";
	        $knowny = "/a/gambiae/genome_vectorbase/gamY";
	} 
	        
	$nr = "/a/nr_NCBI_May-05-2014/nr"; 
}

sub retrieve { 

        system ("blastdbcmd -entry $name -db $input -out seq");  
}        

sub generic_blastn { 
        
        # Subroutine to perform a generic blastn against male and female Illumina databases 
        $which_data = $_[0];
                          
       
        system ("blastn -query seq -out bseq -max_target_seqs 50 -evalue 1e-5 -num_threads 12 -outfmt 6 -db $which_data");

        $lines = 0;
        open(FILE, "bseq" ) or die "Can't open $!";
        while (sysread FILE, $buffer, 4096) {
        
                $lines += ($buffer =~ tr/\n//);
        }
        close FILE;
}

sub sensitive_blastn { 
        
        # Subroutine to perform a sensitive blastn (rCQ) to get understand the close paralogs of the candidate genes 
        $which_data = $_[0];
                  
        system ("blastn -query seq -out bseq -max_target_seqs 10000 -evalue 1e-5 -num_threads 12 -outfmt 6 -db $which_data");

        $lines = 0;
        open(FILE, "bseq" ) or die "Can't open $!";
        while (sysread FILE, $buffer, 4096) {
        
                $lines += ($buffer =~ tr/\n//);
        }
        close FILE;
}


sub rCQ_calculate { 

		# Subroutine to calculate relaxed CQ 
        if ($male_lines > 0) { 
        
                # Calculate the ratio of female_score divided by male_score 
	        $pre_quotient = $female_lines/$male_lines;
	        $rCQ = sprintf("%.3f", $pre_quotient);
	} 
	
	else { 
	
	        $rCQ = "N/A";
        }        							
}

sub blastn_trx { 
        
        # Subroutine to perform a blastn against the Illumina RNA-seq databases 	
         $which_data = $_[0];
                  
         # Male
        system ("blastn -query seq -out bseq -max_target_seqs 1000 -evalue 1e-5 -perc_identity 100 -num_threads 12 -outfmt 6 -db $which_data ");

        $lines = 0;
        open(FILE, "bseq" ) or die "Can't open $!";
        while (sysread FILE, $buffer, 4096) {
        
                $lines += ($buffer =~ tr/\n//);
        }
        close FILE;
}

sub how_repetitive_s_blastn { 

		# Subroutine to determine how many times the candidate gene is present in the genome of interest 
        $which_data = $_[0];
        
        system ("blastn -query seq -out bseq -max_target_seqs 10000 -evalue 1e-25 -num_threads 12 -outfmt 6 -perc_identity 90 -db $which_data");
        
        $lines = 0;
        open(FILE, "bseq" ) or die "Can't open $!";
        while (sysread FILE, $buffer, 4096) {
        
                $lines += ($buffer =~ tr/\n//);
        }
        close FILE;

	# This is to get the percent identity 
	system ("head -5 bseq > bseq2"); 
	open(FILE, "bseq2" ) or die "Can't open $!";
	
	
	@array = ();
	while (<FILE>) { 

		@banshee = split ' ', $_;  
		$perc_identity = "$banshee[1],  $banshee[2],  $banshee[3], $banshee[10] | ";
                push @array,$perc_identity;
	}
        
        $variable = join " ",@array;
	close FILE; 
}   

sub how_repetitive_r_blastn { 

		# Subroutine to determine how many times the candidate gene is present in the genome of interest
        $which_data = $_[0];
        
        system ("blastn -query seq -out bseq -max_target_seqs 10000 -evalue 1e-10 -num_threads 12 -outfmt 6 -db $which_data");
        
        $lines = 0;
        open(FILE, "bseq" ) or die "Can't open $!";
        while (sysread FILE, $buffer, 4096) {
        
                $lines += ($buffer =~ tr/\n//);
        }
        close FILE;
}       

sub blast_trace { 

		# Subroutine to blastn the input the genes of interest against Sanger trace databases from genome projects to see if the gene was not in the assembled genome but present in the traces 
        $which_data = $_[0];
        
        system ("blastn -query seq -out bseq -max_target_seqs 100 -evalue 1e-50 -num_threads 12 -outfmt 6 -db $which_data");
        
        $lines = 0;
        open(FILE, "bseq" ) or die "Can't open $!";
        while (sysread FILE, $buffer, 4096) {
        
                $lines += ($buffer =~ tr/\n//);
        }
        close FILE;
}

sub how_repetitive_tblastx { 

        $which_data = $_[0];
        
        #tblastx
        system ("tblastx -query seq -out bseq -max_target_seqs 1000 -evalue 1e-20 -num_threads 12 -outfmt 6 -db $which_data");
        
        $lines = 0;
        open(FILE, "bseq" ) or die "Can't open $!";
        while (sysread FILE, $buffer, 4096) {
        
                $lines += ($buffer =~ tr/\n//);
        }
        close FILE;
}

sub repeatdb { 

        # Compare the sequence to the repeat scout database 

	$which_data = $_[0];
        
        system ("blastn -query seq -out bseq -max_target_seqs 10000 -evalue 1e-25 -num_threads 12 -outfmt 6 -perc_identity 80 -db $which_data");
        
        $lines = 0;
        open(FILE, "bseq" ) or die "Can't open $!";
        while (sysread FILE, $buffer, 4096) {
        
                $lines += ($buffer =~ tr/\n//);
        }
        close FILE;
}

sub ygenedb { 

        # Compare the sequence to the database of known Y genes 

        $which_data = $_[0];
        
        # Need to add Y gene DB specific sequences 
        system ("blastn -query seq -out bseq -max_target_seqs 10000 -evalue 1e-25 -num_threads 12 -outfmt 6 -perc_identity 80 -db $which_data");
        
        $lines = 0;
        open(FILE, "bseq" ) or die "Can't open $!";
        while (sysread FILE, $buffer, 4096) {
        
                $lines += ($buffer =~ tr/\n//);
        }
        close FILE;
        
        if ($lines > 0) { 
 
		open (F, "bseq") or die "Can't open $!"; 
		while (<F>) { 
	
			# If a sequence matches a Y gene then figure out which Y gene it matches 
			@ar = split ' ', $_;
			$which_y_gene = $ar[1];
		}
	}

	else { 

		$which_y_gene = "-";
	}               
}

sub knownysequences { 

        # Compare the sequence to the database of known An. gambiae Y genes  

        $which_data = $_[0];
        
        # Need to add Y gene DB specific sequences 
        system ("blastn -query seq -out bseq -max_target_seqs 10000 -evalue 1e-25 -num_threads 12 -outfmt 6 -perc_identity 80 -db $which_data");
        
        $lines = 0;
        open(FILE, "bseq" ) or die "Can't open $!";
        while (sysread FILE, $buffer, 4096) {
        
                $lines += ($buffer =~ tr/\n//);
        }
        close FILE;    
}

sub blastx { 

        $which_data = $_[0];

        # Subroutine specific setup 
        $TE_words = 0;
        $blastx_alignments = 0; 
        undef $significant;
        $top5 = '';

        system ("blastx -query seq -out bseq -evalue 1e-5 -num_threads 12 -db $which_data"); 
        
        open (BX, "bseq");
        
        while (<BX>) { 

	        $d = () = ($_ =~ /retrovirus|pol\b|polyprotein|transpon|gag\b|reverse|transcriptase|transposase|mariner/gi);
	        $e = () = ($_ =~ /e-/g);
	        
	        $TE_words = $TE_words + $d;
	        $blastx_alignments = $blastx_alignments + $e;
	        
	        if (/Sequences producing significant alignments/) { 
	        
	                $significant = 1; 
	        }                                       
        }
        
        close BX;       
        
        if ($significant) { 
               
                system ("head -28 bseq | tail -5 > top5");
                open (FH, "top5");
                while (<FH>) { 
                
                        s/\n/\t/g;  
                        $top5 = "$top5" . "$_";
                } 
                close FH;        
        }
        
        else { 
        
                $top5 = "None";
        }
}
