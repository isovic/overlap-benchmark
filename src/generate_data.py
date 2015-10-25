#! /usr/bin/python

import os;
import sys;
import glob;
import math;
import subprocess;
import numpy as np;

SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__));



def peek(fp, num_chars):
	data = fp.read(num_chars);
	if len(data) == 0:
		return '';
	fp.seek(num_chars * -1, 1);
	return data;

def get_single_read(fp):
	lines = '';
	
	line = fp.readline();
	header = line;
	lines += line;
	next_char = peek(fp, 1);
	
	num_lines = 1;
	while len(next_char) > 0 and next_char != lines[0] or (next_char == '@' and num_lines < 4):
		line = fp.readline();
		lines += line;
		next_char = peek(fp, 1);
		num_lines += 1;
		
	return [header.rstrip(), lines.rstrip()];

def CreateFolders(folder_path):
	if not os.path.exists(folder_path):
		sys.stderr.write(('Creating output folder on path: "%s".\n' % (folder_path)));
		os.makedirs(folder_path);

def EstimateCoverageForNumReads(genome_path, genome_filename, mean_read_length, num_reads):
	fp_in = None;
	complete_genome_path = genome_path + '/' + genome_filename + '.fa'

	try:
		fp_in = open(complete_genome_path, 'r');
	except IOError:
		sys.stderr.write(('ERROR: Could not open file "%s" for reading!' % complete_genome_path) + '\n');
		exit(1);
	
	total_genome_length = 0;
	
	while True:
		[header, read] = get_single_read(fp_in);
		
		if (len(header) == 0):
			break;
		
		seq = read[1];
		total_genome_length += len(seq);
		
	fp_in.close();
	
	coverage = int(math.ceil((float(num_reads) / (float(total_genome_length) / float(mean_read_length))))) * 3;
	
	#sys.stderr.write((num_reads) + '\n');
	#sys.stderr.write((total_genome_length) + '\n');
	#sys.stderr.write((mean_read_length) + '\n');
	#print (float(total_genome_length) / float(mean_read_length))
	#print float(num_reads)
	#print coverage
	
	return coverage;

def subsample_alignments_from_fastq(reads_path_prefix, subsampled_set):
	reads_path = reads_path_prefix + '.fq';
	reads_path_fasta = reads_path_prefix + '.fa';
	if (os.path.exists(reads_path) == True):
		complete_reads_path = reads_path_prefix + '-complete_dataset.fq';
		#complete_reads_path = os.path.dirname(reads_path) + '/' + os.path.splitext(os.path.basename(reads_path))[0] + '-complete_dataset' + os.path.splitext(reads_path)[1];
		
		sys.stderr.write(('Renaming file "%s" to "%s"...' % (reads_path, complete_reads_path)) + '\n');
		
		os.rename(reads_path, complete_reads_path);
		
		fp_in = open(complete_reads_path, 'r');
		fp_out = open(reads_path, 'w');

		subsampled_headers = [];
		current_subsample = 0;
		num_read_pairs = 0;
		i = 0;
		while (True):
			[header, read] = get_single_read(fp_in);
			
			if (len(read) == 0):
				break;

			if (i == subsampled_set[current_subsample]):
				fp_out.write(read + '\n');
				subsampled_headers.append(header[1:]);
				current_subsample += 1;
				if (current_subsample >= len(subsampled_set)):
					break;
			
			i += 1;
		
		fp_in.close();
		fp_out.close();

		shell_command = 'rm %s' % (complete_reads_path);
		sys.stderr.write('Removing intermediate file: "%s"\n' % complete_reads_path);
		sys.stderr.write(('Executing command: "%s"\n\n' % shell_command));
		subprocess.call(shell_command, shell=True);

		# sys.stderr.write('Converting the FASTQ file to FASTA...');
		# fastqparser.convert_to_fasta(reads_path, reads_path_fasta);
		# sys.stderr.write('done!\n\n');

		return subsampled_headers;
	else:
		sys.stderr.write('ERROR: Reads file "%s" does not exist!' % (reads_path) + '\n');
		exit(1);

def subsample_alignments_from_sam(reads_path_prefix, subsampled_set):
	sam_path = reads_path_prefix + '.sam';
	if (os.path.exists(sam_path) == True):
		complete_sam_path = reads_path_prefix + '-complete_dataset.sam';
		#complete_reads_path = os.path.dirname(reads_path) + '/' + os.path.splitext(os.path.basename(reads_path))[0] + '-complete_dataset' + os.path.splitext(reads_path)[1];
		
		sys.stderr.write(('Renaming file "%s" to "%s"...' % (sam_path, complete_sam_path)) + '\n');
		
		os.rename(sam_path, complete_sam_path);
		
		fp_in = open(complete_sam_path, 'r');
		fp_out = open(sam_path, 'w');
		
		current_subsample = 0;
		current_num_alignments = 0;
		
		sys.stderr.write(('num_alignments_to_extract = %d' % len(subsampled_set)) + '\n');

		subsampled_qnames = [];

		for line in fp_in:
			if (len(line.strip()) == 0 or line.startswith('@') == True):
				fp_out.write(line);
			else:
				if (current_num_alignments == subsampled_set[current_subsample]):
					fp_out.write(line);
					subsampled_qnames.append(line.split('\t')[0]);
					current_subsample += 1;
					if (current_subsample >= len(subsampled_set)):
						break;

				current_num_alignments += 1;
		
		sys.stderr.write(('current_subsample = %d, subsampled_set[current_subsample] = %d' % (current_subsample, subsampled_set[-1]) ) + '\n');

		fp_in.close();
		fp_out.close();

		shell_command = 'rm %s' % (complete_sam_path);
		sys.stderr.write('Removing intermediate file: "%s"\n' % complete_sam_path);
		sys.stderr.write(('Executing command: "%s"\n\n' % shell_command));
		subprocess.call(shell_command, shell=True);

		return subsampled_qnames;
	else:
		sys.stderr.write(('ERROR: Reads file "%s" does not exist!' % (reads_path)) + '\n');
		exit(1);

def subsample_generated_reads(out_file_prefix, num_reads_to_generate):
	num_generated_reads = 0;
	reads_path = out_file_prefix + '.fq';
	if (os.path.exists(reads_path) == True):
		fp = open(reads_path, 'r');
		while (True):
			[header, read] = get_single_read(fp);
			if (len(read) == 0):
				break;
			num_generated_reads += 1;

		sys.stderr.write('num_generated_reads = %d\n' % num_generated_reads);
		sys.stderr.write('num_reads_to_generate = %d\n' % num_reads_to_generate);
		subsampled_set = sorted(np.random.choice(num_generated_reads, num_reads_to_generate, replace=False));

		subsampled_headers = subsample_alignments_from_fastq(out_file_prefix, subsampled_set);
		subsampled_qnames = subsample_alignments_from_sam(out_file_prefix, subsampled_set);

		# Check the validity of subsampled sets (i.e. if the same reads are included).
		# subsampled_headers = sorted(subsampled_headers);
		# subsampled_qnames = sorted(subsampled_qnames);
		sys.stderr.write('Performing sanity checks on subsampled output...\n');
		if (len(subsampled_headers) != len(subsampled_qnames)):
			sys.stderr.write('ERROR: Subsampling of generated reads and generated SAM file did not produce the same sequences in output! Number of reads in these files differs!\n');
			sys.stderr.write('len(subsampled_headers) = %d\n' % len(subsampled_headers));
			sys.stderr.write('len(subsampled_qnames) = %d\n' % len(subsampled_qnames));
			exit(1);

		current_header = 0;
		current_qname = 0;
		while (current_header < len(subsampled_headers) and current_qname < len(subsampled_qnames)):
			if (subsampled_headers[current_header] != subsampled_qnames[current_qname]):
				sys.stderr.write('ERROR: Subsampling of generated reads and generated SAM file did not produce the same sequences in output! Check the ordering of reads in these files!\n');
				sys.stderr.write('subsampled_headers[current_header] = "%s"\n' % subsampled_headers[current_header]);
				sys.stderr.write('subsampled_qnames[current_qname] = "%s"\n' % subsampled_qnames[current_qname]);
				exit(1);
			current_header += 1;
			current_qname += 1;

		sys.stderr.write('All sanity checks passed!\n');
	else:
		sys.stderr.write('ERROR: Reads not generated! Cannot subsample!\n');
		exit(1);



def GetPBSimRefName(ref_file):
	try:
		fp = open(ref_file, 'r');
	except IOError:
		sys.stderr.write(('ERROR: Could not open file "%s" for reading!' % ref_file) + '\n');
		exit(1);

	header = fp.readline();
	header = header.rstrip();
	fp.close();
	
	if not header.startswith('>'):
		sys.stderr.write(("ERROR: PBsim's ref file does not start with a FASTA header!") + '\n');
		exit(1);

	ref_name = header[1:];
	trimmed_ref_name = ref_name.split()[0];
	
	return [ref_name, trimmed_ref_name];

# --difference-ratio   ratio of differences. substitution:insertion:deletion.
def GeneratePacBio(reference_path, output_path, fold_coverage=20, length_mean=3000, length_sd=2300.0, length_min=100, length_max=25000, accuracy_mean=0.78, accuracy_sd=0.02, accuracy_min=0.75, difference_ratio='10:60:30', num_reads_to_generate=-1):
#	machine_name = 'PacBio';
	# CreateFolders(READS_SIMULATED_ROOT_ABS + '/' + machine_name, genome_filename);
	# CreateFolders(EVALUATION_PATH_ROOT_ABS + '/' + READS_SIMULATED_ROOT + '/' + machine_name, genome_filename);
	
	# complete_genome_path = genome_path + '/' + genome_filename + '.fa'
	complete_genome_path = reference_path;
#	out_file_prefix = output_path + '/' + machine_name + '/' + genome_filename + '/' + 'reads';
	out_file_prefix = '%s/reads' % (output_path);
	CreateFolders(output_path);
	
	simulator_path = '%s/pbsim-1.0.3-Linux-amd64' % (SCRIPT_PATH);
	simulator_bin = '%s/pbsim-1.0.3-Linux-amd64/Linux-amd64/bin/pbsim' % (SCRIPT_PATH);
	maf_convert_path = '%s/last-534/scripts/maf-convert' % (SCRIPT_PATH);

	if (not os.path.exists(simulator_bin)):
		sys.stderr.write('ERROR: Could not find the PBsim binaries! Run "%s setup" first. Exiting.\n' % (sys.argv[0]));
		exit(1);
	
	final_sam_file = out_file_prefix + '.sam';
	fp = open(final_sam_file, 'w');
	fp.close();
	
	final_fastq_file = out_file_prefix + '.fq';
	fp = open(final_fastq_file, 'w');
	fp.close();
# SIMULATOR_PATH=/home/ivan/work/eclipse-workspace/golden-bundle/tools/pbsim-1.0.3-Linux-amd64
# REFERENCE_PATH=/home/ivan/work/eclipse-workspace/golden-bundle/reference-genomes/escherichia_coli.fa
# $SIMULATOR_PATH/Linux-amd64/bin/pbsim --data-type CLR --depth 20 --model_qc $SIMULATOR_PATH/data/model_qc_clr --seed 1234567890 --prefix $OUT_FILE_PREFIX $REFERENCE_PATH
# /last-460/scripts/maf-convert.py sam $MAF_FILE $SAM_FILE

	#random_seed = '1234567890';
	random_seed = '32874638';

	# Data type:
	#	Continuous Long Read (CLR) : long and high error rate.
	#	Circular consensus Read (CCS) : short and low error rate.
	# --difference-ratio   ratio of differences. substitution:insertion:deletion.
	shell_command = simulator_bin + r' --data-type CLR --depth ' + str(fold_coverage) + \
					' --model_qc ' + simulator_path + \
					'/data/model_qc_clr ' + \
					' --length-mean ' + str(length_mean) + \
					' --length-sd ' + str(length_sd) + \
					' --length-min ' + str(length_min) + \
					' --length-max ' + str(length_max) + \
					' --accuracy-mean ' + str(accuracy_mean) + \
					' --accuracy-sd ' + str(accuracy_sd) + \
					' --accuracy-min ' + str(accuracy_min) + \
					' --difference-ratio ' + difference_ratio + \
					' --prefix ' + out_file_prefix + \
					r' ' + complete_genome_path;
				
					#' --seed ' + random_seed + \
	
	sys.stderr.write(('Simulating PacBio reads using PBsim') + '\n');
	sys.stderr.write(('Executing command: "%s"' % shell_command) + '\n');
	
	#exit(1);
	subprocess.call(shell_command, shell=True);

	sys.stderr.write((' ') + '\n');
	
	sys.stderr.write(('Converting generated *.maf files to corresponding SAM files.') + '\n');

	maf_files = glob.glob(out_file_prefix + '*.maf');
	maf_files = sorted(maf_files);
	sys.stderr.write('\n'.join(maf_files) + '\n');
	
	for maf_file in maf_files:		
		# Convert the maf file to SAM format.
		sam_file = maf_file[0:-3] + 'sam';
		shell_command = maf_convert_path + ' sam ' + maf_file + ' > ' + sam_file;
		sys.stderr.write(('Converting MAF to SAM ("%s" -> "%s")' % (maf_file, sam_file)) + '\n');
		sys.stderr.write(('Executing command: "%s"' % shell_command) + '\n');
		subprocess.call(shell_command, shell=True);
		
		# Use Bash's sed to replace the 'ref' keyword in the generated maf files with the actual name of the reference sequence, and concatenate all the separate SAM files to one file.
		[reference_name, trimmed_reference_name] = GetPBSimRefName(maf_file[0:-3] + 'ref');
		# Here we escape the special characters so that SED command runs properly if any of these characters should to appear in a FASTA header.
		escape_chars = r'\/()[].*^$';
		reference_name = ''.join([('\\' + char) if char in escape_chars else char for char in reference_name]);
		shell_command = r'cat ' + sam_file + r" | sed 's/^\(.*\)ref/\1" + reference_name + r"/' >> " + final_sam_file
		sys.stderr.write(('Replacing PBsim\'s "ref" keyword with actual FASTA header') + '\n');
		sys.stderr.write(('Executing command: "%s"' % shell_command) + '\n');
		subprocess.call(shell_command, shell=True);
		
		fastq_file = maf_file[0:-3] + 'fastq';
		shell_command = r'cat ' + fastq_file + ' >> ' + final_fastq_file;
		sys.stderr.write(('Concatenating FASTQ file to the total reads file ("%s" -> "%s")' % (fastq_file, final_fastq_file)) + '\n');
		sys.stderr.write(('Executing command: "%s"' % shell_command) + '\n');
		subprocess.call(shell_command, shell=True);
		sys.stderr.write((' ') + '\n');
		
		sys.stderr.write((' ') + '\n');

		shell_command = 'rm %s' % (maf_file);
		sys.stderr.write('Removing intermediate file: "%s"\n' % maf_file);
		sys.stderr.write(('Executing command: "%s"\n\n' % shell_command));
		subprocess.call(shell_command, shell=True);
		shell_command = 'rm %s' % (sam_file);
		sys.stderr.write('Removing intermediate file: "%s"\n' % sam_file);
		sys.stderr.write(('Executing command: "%s"\n\n' % shell_command));
		subprocess.call(shell_command, shell=True);
		shell_command = 'rm %s' % (fastq_file);
		sys.stderr.write('Removing intermediate file: "%s"\n' % fastq_file);
		sys.stderr.write(('Executing command: "%s"\n\n' % shell_command));
		subprocess.call(shell_command, shell=True);

		ref_file = maf_file[0:-3] + 'ref';
		shell_command = 'rm %s' % (ref_file);
		sys.stderr.write('Removing intermediate file: "%s"\n' % ref_file);
		sys.stderr.write(('Executing command: "%s"\n\n' % shell_command));
		subprocess.call(shell_command, shell=True);

#	if (GENERATE_FIXED_AMOUNT_OF_READS == True):
	# if (num_reads_to_generate > 0):
	# 	ExtractNReadsFromFile(out_file_prefix, num_reads_to_generate);
	# 	ExtractNAlignmentsFromSAM(out_file_prefix, num_reads_to_generate);
	if (num_reads_to_generate > 0):
		subsample_generated_reads(out_file_prefix, num_reads_to_generate);

def GenerateGridTest(reference_path, out_path, coverage=30, error_rates=[0.0, 0.05, 0.10, 0.15, 0.20]):
	##### OXFORD NANOPORE DATA #####
	# --difference-ratio   ratio of differences. substitution:insertion:deletion.
	# GenerateOxfordNanoporeFromObservedStatistics('caenorhabditis_elegans', num_reads_to_generate=num_reads_to_generate);

	reference_path = os.path.abspath(reference_path);
	out_path = os.path.abspath(out_path);

	reference_basename = os.path.splitext(os.path.basename(reference_path))[0];

	# sys.stderr.write(('num_reads_to_generate = %d' % num_reads_to_generate) + '\n');

	for error_rate in error_rates:
		# coverage = EstimateCoverageForNumReads(REFERENCE_GENOMES_ROOT_ABS, genome_filename, mean_read_length, num_reads_to_generate) + 1;
		machine_name = '%s-e%.2f' % (reference_basename, error_rate);
		output_folder = '%s/%s' % (out_path, machine_name);

		# These are simulations for difference_ratio obtained from LAST's alignments of E. Coli R7.3 reads (from Nick Loman).
		# 1d reads:
		# [CIGAR statistics - individual indels]
		#                       	mean	std	median	min	max
		# Error rate stats:     	0.41	0.05	0.40	0.10	0.60
		# Insertion rate stats: 	0.05	0.02	0.05	0.00	0.23
		# Deletion rate stats:  	0.16	0.05	0.15	0.00	0.49
		# Mismatch rate stats:  	0.20	0.03	0.20	0.03	0.32
		# Match rate stats:     	0.75	0.04	0.75	0.59	0.97
		# Read length stats:    	3629.76	3294.04	2438.00	57.00	31299.00
		# Difference ratio: 51:11:38 (mismatch:insertion:deletion)
		GeneratePacBio(reference_path, output_folder, fold_coverage=coverage, length_mean=5600, length_sd=3500, length_min=50, length_max=100000,
																	accuracy_mean=(1.0 - error_rate), accuracy_sd=0.09, accuracy_min=(1.0 - 0.60), difference_ratio='55:17:28',
																	num_reads_to_generate=0);


def download_and_install():
	sys.stderr.write('Downloading and unpacking PBsim.\n');
	command = 'cd %s; wget http://pbsim.googlecode.com/files/pbsim-1.0.3-Linux-amd64.tar.gz; tar -xzvf pbsim-1.0.3-Linux-amd64.tar.gz' % (SCRIPT_PATH);
	subprocess.call(command, shell='True');
	sys.stderr.write('\n');

	sys.stderr.write('Downloading and unpacking LAST aligner. Its scripts are needed to convert from MAF to SAM.\n');
	command = 'cd %s; wget http://last.cbrc.jp/last-534.zip; unzip last-534.zip' % (SCRIPT_PATH);
	subprocess.call(command, shell='True');
	sys.stderr.write('\n');

def verbose_usage_and_exit():
	sys.stderr.write('Simulates reads from a given genome.\n\n');
	sys.stderr.write('Usage:\n');
	sys.stderr.write('\t%s mode [<reference_file> <output_path> coverage]\n' % sys.argv[0]);
	sys.stderr.write('\n');
	sys.stderr.write('\t- mode - either "run" or "setup". Is "install" other parameters can be ommitted.\n');

	exit(0);

if __name__ == "__main__":
	if (len(sys.argv) < 2 or len(sys.argv) > 7):
		verbose_usage_and_exit();

	if (sys.argv[1] == 'setup'):
		download_and_install();
		exit(0);

	elif (sys.argv[1] == 'run'):
		if (len(sys.argv) != 5):
			verbose_usage_and_exit();

		reference_file = sys.argv[2];
		output_path = sys.argv[3];
		coverage = int(sys.argv[4]);

		# GenerateGridTest('%s/../genomes/escherichia_coli.fa' % (SCRIPT_PATH), '%s/../reads-simulated/' % (SCRIPT_PATH), coverage=30, error_rates=[0.0]);
		GenerateGridTest(reference_file, output_path, coverage=coverage);

	else:
		verbose_usage_and_exit();
