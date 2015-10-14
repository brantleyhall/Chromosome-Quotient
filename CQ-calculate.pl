#!/usr/bin/perl -w   
# CQ-calculate.pl
# Locate Y Chromosome Sequences in an Assembled Genome Using Sex-Specific Illumina Databases 
# Brantley Hall, October 2011 
#
# System Requirements 
# 1. A computer with a Unix/Linux operating system (Tested on Fedora 16, CentOS, and Ubuntu) 
# 2. Perl (http://www.perl.org) 
# 3. Bowtie installed and in your path (http://bowtie-bio.sourceforge.net/index.shtml)
#
# To run CQ-calculate you need the following data
# 1. A multi-fasta file containing the reference sequences  
# 2. A male-specific illumina database in fastq or fasta format (Preferably fastq)  
# 3. A female-specific illumina database in fastq or fasta format (Preferably fastq)
#
# Summary 
# 1. Check for the existence of the bowtie index for the input file.    
# 2. Run bowtie against the assembly of interest with male and female databases and optionally transcriptome databases
# 3. Process the bowtie output files
# 4. Calculate the number of male and female hits for each reference sequence 
# 5. Determine the chromosome quotient for each reference sequence   
#
# Output 
# 1. Default output is in CQ-calculate_output_(-n) (Default test) 
# 2. The output suffix can be changed with the -n option 
# 3. The columns for the output file are [Sequence ID] [Female Hits] [Male Hits] [Chromosome Quotient]
#
# USAGE:
# perl CQ-calculate.pl [options]
# Mandatory Options 
# -i [multi-fasta file containing reference sequences]
# -f [female illumina database]
# -m [male illumina database]
# Advanced Options 
# -n [output suffix] (Default test)
# -t [minimum number of hits] (Default 20, must be greater than 0) 
# -l [length of illumina reads] (Default 84) 
# -p [number of processors for bowtie] (Default 1)	
# -fa [use FASTA data instead of FASTQ data]
# -k [keep bowtie output files]
# -norm [normalize the results by providing the mean chromosome quotient of autosomal sequences] 
# -log [print bowtie log]  
# -h [print help message]
# -help [print help message] 
# 
# Example: 
# perl CQ-calculate.pl -i example_genome -f female_illumina.fq -m male_illumina.fq -n example -p 8 -t 50
# File containing results with be CQ-calculate_output_example 

use Getopt::Long;
use File::Basename;
use Bio::Seq;
use Bio::SeqIO;

# Set the default values of advanced variables 
$name = "test"; # The deafult output suffix is test 
$alignment_threshold = 10; # The default threshold for male hits is 20, for more accurate results this should be increased 
$format = "q"; # The default format is fastq can be changed to fasta with -fa option  
$num_processors = 1; # The default number of processors is 1. If your computer has more than one processor, this should should be increased to improve performance.  

GetOptions( "i=s"	=> \$input,
	    "m=s"	=> \$male_database,
	    "f=s"	=> \$female_database,
	    "t=s"	=> \$trx_database,
	    "n=s"       => \$name,
	    "a=i"	=> \$alignment_threshold,
	    "norm=f"	=> \$normal, 		
	    "p=i"	=> \$num_processors,
	    "fa"	=> \$fasta,
	    "log"	=> \$log,
	    "k"		=> \$keep, 			 	
	    "h"		=> \$help,
	    "help"	=> \$help 	
	   );	


# Check for problems with the command specified by the user 
&setup; 

# Test to see if the bowtie index for the input file exists 
&ebwt_exist; 

# Run bowtie against the male database
&run_bowtie("m"); 

# Parse the male bowtie output 
&count_everything("m_vs_$filename"); 

# Run bowtie against the female database
&run_bowtie("f");

# Parse the female bowtie output 
&count_everything("f_vs_$filename");

# Parse the trx 
if ($trx_database) { 

        # Run bowtie against the male database
        &run_bowtie("t");         
        
	# Parse the transcriptome database 
	&count_everything("t_vs_$filename");
}

# Create fequency hashes 
&hash_setup; 

# Count the length of all the sequences in the input file 
&get_length; 

# Calculate the chromosome quotient for each sequence in the input file 
&CQ_calculate;
 
print "Finished\n\n"; 
 
exit 0; 

#################
#  Subroutines  #
#################

sub setup { 

	# Set up up the help documentation 
	$help_message = "\nUSAGE: perl CQ-calculate.pl [options]\n\nOptions:\n\t-i\t [multi-fasta file containing reference sequences]\n\t\-f\t [female illumina database]\n\t-m\t [male illumina database]\n\t-c\t [use upper CQ output option]\n\t-r\t [use range output option]\n\t\-n\t [output suffix]\n\t-t\t [minimum number of hits]\n\t-l\t [length of illumina reads]\n\t-p\t [number of processors]\n\t-fa\t [use fasta format illumina data]\n\t-k\t [keep bowtie output files]\n\t-norm\t [normalize the results by providing the mean chromosome quotient of autosomal sequences]\n\t-log\t [print bowtie log]\n\t-h\t [print help message]\n\t-help\t [print help message]\n"; 

	# List of mandatory options when user does not specify all mandatory options
	$mandatory_options = "\nYou seem to be missing a mandatory option.\n\nUSAGE: perl CQ-calculate.pl [options]\nMandatory Options:\n\t-i\t [multi-fasta file of sequences to be analyzed]\n\t\-f\t [female illumina database]\n\t-m\t [male illumina database]\n\nType CQ-calculate.pl -h for detailed usage.\n";

	if ($help) { 

		print "$help_message\n"; 
		exit 0;
	}
	 
	# Check to see if mandatory options are defined by user 
	if (! $input || ! $male_database || ! $female_database) { 

		print "$mandatory_options\n";
		exit 0; 
	}

	# Check to see if the input filq, male, and female databases specified by the user exist
	if (! -e $input) { 
	
		print "\nThe input file: $input does not exist.\n\n"; 
		exit 0; 
	} 

	if (! -e $male_database) { 
	
		print "\nThe male database: $male_database does not exist.\n\n"; 
		exit 0; 
	} 

	if (! -e $female_database) { 
	
		print "\nThe female database: $female_database does not exist.\n\n";
		exit 0;  
	}

	if ($alignment_threshold <= 0) { 

		print "\nThe threshold (-t) must be greater than zero.\n\n";
		exit 0; 
	} 

	# Retrieve the basename of the input file 	
	$filename = basename($input);
	
	
	print "\n";
}

sub ebwt_exist { 

	# If the bowtie index of the input file does not exist make the index with bowtie-build 
	
	$bowtie_index = $input . ".1.ebwt"; 
	if (! -e $bowtie_index) { 
	
		print "Running bowtie-build...\n";
		system ("bowtie-build $input $input > torm");
		unlink "torm"; 
	} 
} 
		
sub run_bowtie { 

	# Run bowtie with the options supplied by the user 

    # Retrieve the input argument 
    $which_data = $_[0];

	# Option to use fasta instead of fastq format 	
	if ($fasta) { 
	
		$format = "f"; 
	}

        if ($which_data eq "m") {
        
                $database = $male_database;
                $bowtie_output_file = "m_vs_$filename";
        }
        
        if ($which_data eq "f") {
        
                $database = $female_database; 
                $bowtie_output_file = "f_vs_$filename";
        }
        
        if ($which_data eq "t") {
        
                $database = $trx_database; 
                $bowtie_output_file = "t_vs_$filename";
        } 
        
        print "Aligning reference to $database...\n";       

	# Run bowtie against the male database 
	system ("bowtie -a -p $num_processors -v 0 $input --suppress 1,2,4,5,6,7,8,9 -$format $database $bowtie_output_file 2> bowtie_log"); 
	
	# If requested by the user print the bowtie output 
	if ($log) { 

		open (BT, "bowtie_log"); 
		
		while (<BT>) { 

			print $_; 
		}
		print "\n"; 		
		close BT; 
	}	  				
}

sub count_everything { 

	# Open the bowtie output and count the occurrence of every subject ID * Check name 
	# Outputs the format: 
	# Sequence ID	Number of Occurrences 
	# Example: 
	# contig1	45 
	
	print "Counting alignments...\n\n";
	 
	$count_input = $_[0];  
	$count_output = "cc_" . $count_input;  	
 
	open (COUNT_IN, "$count_input") or die "$!"; 
	open (COUNT_OUT, ">$count_output") or die "$!"; 

	# Reset the count hash 	
	%count_hash = ();	

	while(<COUNT_IN>) {
		
		map {$count_hash{$_}++} (split /\s+/)
	}
	print COUNT_OUT "$_ $count_hash{$_}\n" foreach sort keys %count_hash;

	close COUNT_IN; 
	close COUNT_OUT;
	
	# Remove bowtie output file to save space 
	if (! $keep) { 	
		
		unlink "$bowtie_output_file";
	}
}

sub hash_setup { 

	# Create three hashes that contain the sequence ID and the number of hits
	# Add keys and values to male hash  
	open (MALE_IN, "cc_m_vs_$filename") or die "$!"; 		
	while (<MALE_IN>) {

		@msplit = split ' ', $_; 
		$key = $msplit[0]; 
		$value = $msplit[1]; 
	
		$male_hash{$key} = $value; 
	}  
	close MALE_IN;

	# Add keys and values to female hash  
	open (FEMALE_IN, "cc_f_vs_$filename") or die "$!"; 		
	while (<FEMALE_IN>) {

		@fsplit = split ' ', $_; 
		$key = $fsplit[0]; 
		$value = $fsplit[1]; 
	
		$female_hash{$key} = $value; 
	}  
	close FEMALE_IN;
	
	if ($trx_database) { 
	
		# Add keys and values to trx hash  
		open (TRX_IN, "cc_t_vs_$filename") or die "$!"; 		
		while (<TRX_IN>) {
	
			@tsplit = split ' ', $_; 
			$key = $tsplit[0]; 
			$value = $tsplit[1]; 
		
			$trx_hash{$key} = $value; 
		}  
		close TRX_IN;
	}		
}

sub get_length { 

	print "Counting Length...\n";

	$seqio_obj = Bio::SeqIO->new(-file => "$input", -format => "fasta" ) or warn "$!";

	while ($seq_obj = $seqio_obj->next_seq) {
	
		$name_count = $seq_obj->id;	
		$length_count = $seq_obj->length;

		$count_hash{$name_count} = $length_count;
	}
}

sub CQ_calculate {

	# Calculate the ratio of female to male alignments 

	print "Calculating Chromosome Quotients...\n";
		
	# Open the output file 
	$passed_output = "CQ-calculate_output_" . "$name";
	open (OUTFILE, ">$passed_output") or die "$!";

	# Open the assembly of interest to find reference sequence names 	   
	open (reference_IN, "$input") or die "$!";

	# Open the assembly of interest and examine calculate the CQ of each seuqence 
	while (<reference_IN>) { 

		if (/\>/) { 

			# Find the reference sequence name 	
			# Remove > 		
			$reference_name_line = substr $_, 1;
			# Remove any descriptions		
			@line = split ' ', $reference_name_line; 
			$reference_name = $line[0];
		
			# If the reference sequence name exists in the male hash return the key 
			if (exists($male_hash{"$reference_name"})) { 
		
				$male_score = $male_hash{"$reference_name"};
			}
		
			# If the reference sequence name does not exist in the has then assign a value of 0
			else {
			
				$male_score = 0; 
			}  
		
			# If the reference sequence name exists in the female hash return the key 
			if (exists($female_hash{"$reference_name"})) { 
		
				$female_score = $female_hash{"$reference_name"};
			}
		
			# If the reference sequence name does not exist in the has then assign a value of 0
			else {
			
				$female_score = 0; 
			}  

			if ($trx_database) { 
				
				# If the reference sequence name exists in the female hash return the key 
				if (exists($trx_hash{"$reference_name"})) {
			
					$trx_score = $trx_hash{"$reference_name"};
				}
			
				# If the reference sequence name does not exist in the has then assign a value of 0
				else {
				
					$trx_score = 0; 
				}
			}

			# If the reference sequence name exists in the count hash return the key 
			if (exists($count_hash{"$reference_name"})) { 
		
				$length = $count_hash{"$reference_name"};
			}
		
			# If the reference sequence name does not exist in the has then assign a value of 0
			else {
			
				$length = 0; 
			} 
		
			# To prevent dividing by 0 CQ-calcalculate only processes reference seuqences with male_score greater than a defined threshold. (Default 20) 
			if ($male_score >= $alignment_threshold) { 
		
				# Calculate the ratio of female_score divided by male_score 
				$pre_quotient = $female_score/$male_score;
								
				# Round the ratio to three places past the decimal 
				if (! $normal) { 
			
					$quotient = sprintf("%.3f", $pre_quotient);
				}

				# Option to normalize data to known autosomal chromosome quotient 
				if ($normal) {  

					$normal_quotient = $pre_quotient/$normal; 
					$quotient = sprintf("%.3f", $normal_quotient);
				} 
		
				if (! $trx_database) { 
				
					# Print all the output 
					print OUTFILE "$reference_name\t $length\t $female_score\t $male_score\t $quotient\n";
				}
				
				if ($trx_database) { 
				
					# Print all the output 
					print OUTFILE "$reference_name\t $length\t $female_score\t $male_score\t $quotient\t $trx_score\n";
				}
			}
		} 
	}
	
	close reference_IN;
}



