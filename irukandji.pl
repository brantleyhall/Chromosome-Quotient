#!/usr/bin/perl -w

# Irukandji v1 
# Brantley Hall, November 2013 

# Irukandji is a new way to find Y chromosome sequences using unqiue male kmers
# It works by comparing male and female Illumina sequencing databases and finding kmers specific to the male database
# Then it finds sequence with greater than a threshold of male-specific kmers and reports them as Y chromosome sequences

# Compared to CQ-calculate, Irukandki can unambiguously identify identify repetitive Y sequences 

# Inputs 
# 1. multi-fasta file that will be searched for unique male kmers 
# 2. file of unique male kmers Format: kmer-sequence frequency #### Make sure that there are reverse complements of all these kmers 
# Output 
# 1. The number of male-specific kmers that align to each input sequence 

# No Split Option 
# 1. Read a file of male-specific kmers into a hash
# 2. Read in a multi-fasta file and split into individual seqs 
# 3. Read input multi-fasta file into seq_hash where the name is the key and the sequence is the value   
# 4. One sequence from the imput multi-fasta file at a time, split the sequence into kmers, and read into a hash removing duplicated  
# 5. Compare this hash to the msk_hash and count the number of matches to the msk_hash  

# Split option 
# The perl substr function is extremely slow when dealing with long strings 
# Therefore, I have made an option to split the input sequence into many 3,000 bp pieces so substring will be faster 

# 1. Read a file of male-specific kmers into a hash
# 2. Read in a multi-fasta file and split into individual seqs 
# 3. Read input multi-fasta file into seq_hash where the name is the key and the sequence is the value  

# Exceptions 
# Input sequences in the multi-fasta file can't have the same names 

# Set the default value for kmer length
$k = 20;  

use Getopt::Long;
use File::Basename;
use Bio::Perl;

GetOptions( "i=s"	=> \$input,
	    "o=s"	=> \$output,
	    "m=s"	=> \$msks,
       	"s"   	=> \$split,
	    "ps"	=> \$pre_split		
	   );

&setup; 

&process_msks;

&input_setup; 

if (! $split && ! $pre_split) { 
	
	&compare_kmers;
} 

if ($split) { 

	&compare_kmers_split;
} 

if ($pre_split) { 

	print "Pre-split\n";
	&compare_kmers_split; 
}

exit; 

#################
#  Subroutines  #
#################

sub setup { 

	# Open the output 
	open (OUT, ">$output") or die "Can't open $output: $!";	
}

sub process_msks { 

	# Open the file of male-specific kmers and read it into a hash 

	print "Processing msks\n";
	# Open the msk file 
	open (MIN, $msks) or die "male-specific kmer file not found\n $!"; 
	
	while (<MIN>) { 

		# Read the msk file into a hash 
		@msk_array = split ' ', $_;

		$key = $msk_array[0];
		$value = ''; 
		$msk_hash{$key} = $value;  
	} 	
	close MIN; 	
} 


sub input_setup { 

    # Retrieve the basename of the input file 	
	$filename = basename($input);

	# To fix slowness for long sequences there is a split option here
	# It used bp_split_seq.pl that I have modified slightly to output to a single file and not split sequences short than split 
	if (! $split) { 
	
        &process_input("$input"); 

        }	

	if ($split) { 

		print "Spitting input sequences\n";

        	# Remove only irukandji split 
        	system ("rm -rf split_$filename.fa"); 
	
        	# Split the input file into smaller chunks
        	system ("/home/brantley/db/sw/modified_bp_split_seq.pl -S 1 -f $input -c 3000 -out split_$filename -o 50> split_log"); # The -o 50 here is extremely important. 
    	
        	&process_input("split_$filename");            
	}	
}

sub process_input { 
	
	# Get the input from the subroutine invocation 
   	$pass_input = $_[0];      

	print "Processing input\n";
	
	# Reads the sequences to be cut into kmers into a hash removing the >header line 
	%seq_hash = (); # key = seq_name, value = seq;

	# redefine the record separator
	local $/ = ">";
	open (IN, $pass_input) or die "Can't open $input: $!";
	$in_line = <IN>; # toss the first record

	while ( $in_line = <IN> ) {

		chomp $in_line; # remove the ">" character in the end 
		my ($seq_name, $seq) = split( /\n/, $in_line, 2 );
		$seq =~ tr/\t\n\r//d;    # Remove whitespace
		$seq_hash{$seq_name} = $seq; # Coverts everything to uppercase 
	}
	close IN;
}

sub compare_kmers { 
	
	print "Generating kmers and comparing\n";

	# The memory footprint is decreased by processing the sequences one by one  
	foreach $seq_name (keys %seq_hash) { 
	    
		chomp $k; # where k is the length of the kmer to be used 

		%in_hash = (); 

		# Get the length of the sequence in the seq_hash and cut into kmers while shorter than the length 
		# Add these to a single hash for the sequence to prevent duplicates  
		while (length($seq_hash{$seq_name})  >= $k) {
		
			# Covert everything to uppercase 
			$uppercase_sequence = uc($seq_hash{$seq_name}); 
	
			# substring the first 20 bp to get a kmer 
			$mer = substr($uppercase_sequence, 0, $k);   
		
			# If the kmer already exists in the kmer hash 
			if (! exists($in_hash{"$mer"})) { 

				# Add the kmer and contig name to the hash 
				$in_hash{$mer} = ''; 
			} 	
		
			# Remove the first bp	        
			$seq_hash{$seq_name}= substr($seq_hash{$seq_name}, 1, length($seq_hash{$seq_name})-1);
		}

		$count = 0; # Reset count 
		# Check each unique kmer from this fasta sequence to see if it is in the msk_hash 
		foreach $umer (keys %in_hash) { 		

			if (exists($msk_hash{"$umer"})) {  
		
				$count += 1; 
			} 
		}

    		if ($count > 0) { # Replace this with a threshold soon...  
    
    			print OUT "$seq_name\t $count\n";
    		}
    
        }
	close OUT;
}

sub compare_kmers_split { 
	
	# This is a special subroutine for when split is active 
	# The input sequences are processed as normal 
	# Instead of $count determining how many msks match the input sequence, all the msks that match are put into a hash with the (split) sequence name 
	# This way, all the sequences that match a single input sequence (before the split) will be put into a single hash 
	# Then I split this hash into an array by _irukandjisplit and the read that into a hash so it is unique 
	# Then I count the number of entries in that hash 	

	print "Generating kmers and comparing (Split)\n";

    	%kmerlist_hash = ();

	# The memory footprint is decreased by processing the sequences one by one  
	foreach $seq_name (keys %seq_hash) { 
	    
		chomp $k; # where k is the length of the kmer to be used 

		%in_hash = (); 

		# Get the length of the sequence in the seq_hash and cut into kmers while shorter than the length 
		# Add these to a single hash for the sequence to prevent duplicates  
		while (length($seq_hash{$seq_name})  >= $k) {
		
			# Covert everything to uppercase 
			$uppercase_sequence = uc($seq_hash{$seq_name}); 
	
			# substring the first 20 bp to get a kmer 
			$mer = substr($uppercase_sequence, 0, $k);   
		
			# If the kmer already exists in the kmer hash 
			if (! exists($in_hash{"$mer"})) { 

				# Add the kmer and contig name to the hash 
				$in_hash{$mer} = ''; 
			} 	
		
			# Remove the first bp	        
			$seq_hash{$seq_name}= substr($seq_hash{$seq_name}, 1, length($seq_hash{$seq_name})-1);
		}

		# Split off the _coords-coords generated by bp_split_seq
		@split_array = split /_irukandjisplit/, $seq_name; 
		    
		# Have a hash of previous seq names, and if it is in that hash, add the count 
		$split_name = $split_array[0];

		$count = 0; # Reset count 
		# Check each unique kmer from this fasta sequence to see if it is in the msk_hash 
		foreach $umer (keys %in_hash) { 		

			if (exists($msk_hash{"$umer"})) { 				

				$count += 1;

				if (exists($kmerlist_hash{"$split_name"})) {  
				 
					$mer_list = $kmerlist_hash{$split_name};
					$mer_list = $mer_list . "_irukandjisplit" . $umer;
					$kmerlist_hash{$split_name} = $mer_list; 
				}

				 if (! exists($kmerlist_hash{"$split_name"})) {  
				 
					$kmerlist_hash{$split_name} = $umer;
				}
			} 
		}
	}

	# After reading through everything, print out the matches 
        # read through the hash of seq_names and print everything out
        foreach $entry (keys %kmerlist_hash) {

		%splitk_hash = ();  

        	@howmany_array = split /_irukandjisplit/, $kmerlist_hash{$entry}; 
		
		foreach $bigk (@howmany_array) { 

			$splitk_hash{$bigk} = ''; 
		}
		$length = keys (%splitk_hash); 		
		print OUT "$entry\t$length\n";
	} 
	close OUT;
}
