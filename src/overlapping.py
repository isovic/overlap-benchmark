#! /usr/bin/python

import os;
import sys;
import glob;
import math;
import subprocess;
import numpy as np;

SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__));
UTILITY_PATH = os.path.join(os.path.dirname(os.path.realpath(__file__)), '../utility');
SAM_TO_M4_BIN = '%s/blasr/pbihdfutils/bin/samtom4' % (UTILITY_PATH);
TOOLS_PATH = os.path.join(os.path.dirname(os.path.realpath(__file__)), '../tools');

DRY_RUN = False;

#################################################
#################################################
#################################################
def execute_command(dry_run, command):
	sys.stderr.write('[Executing] "%s"\n' % command);
	if (dry_run == False):
		subprocess.call(command, shell=True);
	sys.stderr.write('\n');

def execute_command_get_stdout(command):
	sys.stderr.write('[Executing] "%s"\n' % (command));
	p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
	[out, err] = p.communicate()
	sys.stderr.write('\n');
	return [out, err];

def measure_command_wrapper(out_filename):
	return '/usr/bin/time --format "Command line: %%C\\nReal time: %%e s\\nCPU time: -1.0 s\\nUser time: %%U s\\nSystem time: %%S s\\nMaximum RSS: %%M kB\\nExit status: %%x" --quiet -o %s ' % out_filename;
#################################################
#################################################
#################################################

def setup_graphmap():
	sys.stderr.write('Cloning GraphMap.\n');
	execute_command(DRY_RUN, 'cd %s; git clone https://github.com/isovic/graphmap.git && cd graphmap && make -j' % (TOOLS_PATH));
	sys.stderr.write('Done.\n\n');

def setup_mhap():
	sys.stderr.write('Cloning MHAP.\n');
	execute_command(DRY_RUN, 'cd %s; git clone https://github.com/marbl/MHAP.git && cd MHAP && ant' % (TOOLS_PATH));
	sys.stderr.write('Done.\n\n');

def setup_daligner():
	pass;

def setup_minimap():
	sys.stderr.write('Cloning Minimap.\n');
	execute_command(DRY_RUN, 'cd %s; git clone https://github.com/lh3/minimap.git && cd minimap && make -j' % (TOOLS_PATH));
	execute_command(DRY_RUN, 'cd %s; git clone https://github.com/lh3/miniasm.git && cd miniasm && make -j' % (TOOLS_PATH));
	sys.stderr.write('Done.\n\n');
#################################################
#################################################
#################################################

def run_graphmap(reads_file, out_overlaps_file):
	memtime_file = os.path.splitext(out_overlaps_file)[0] + '.memtime';
	bin_file = '%s/graphmap/bin/Linux-x64/graphmap' % (TOOLS_PATH);
	execute_command(DRY_RUN, '%s %s -w owler -r %s -d %s -o %s -L mhap' % (measure_command_wrapper(memtime_file), bin_file, reads_file, reads_file, out_overlaps_file));

def run_mhap_default(reads_file, out_overlaps_file):
	memtime_file = os.path.splitext(out_overlaps_file)[0] + '.memtime';
	bin_file = '%s/MHAP/target/mhap-1.6.jar' % (TOOLS_PATH);
	execute_command(DRY_RUN, '%s java -Xmx32g -server -jar %s -s %s > %s' % (measure_command_wrapper(memtime_file), bin_file, reads_file, out_overlaps_file));

def run_mhap_nanopore_fast(reads_file, out_overlaps_file):
	memtime_file = os.path.splitext(out_overlaps_file)[0] + '.memtime';
	bin_file = '%s/MHAP/target/mhap-1.6.jar' % (TOOLS_PATH);
	execute_command(DRY_RUN, '%s java -Xmx32g -server -jar %s -s %s --nanopore-fast > %s' % (measure_command_wrapper(memtime_file), bin_file, reads_file, out_overlaps_file));

def run_mhap_pacbio_fast(reads_file, out_overlaps_file):
	memtime_file = os.path.splitext(out_overlaps_file)[0] + '.memtime';
	bin_file = '%s/MHAP/target/mhap-1.6.jar' % (TOOLS_PATH);
	execute_command(DRY_RUN, '%s java -Xmx32g -server -jar %s -s %s --pacbio-fast > %s' % (measure_command_wrapper(memtime_file), bin_file, reads_file, out_overlaps_file));

def run_mhap_pacbio_sensitive(reads_file, out_overlaps_file):
	memtime_file = os.path.splitext(out_overlaps_file)[0] + '.memtime';
	bin_file = '%s/MHAP/target/mhap-1.6.jar' % (TOOLS_PATH);
	execute_command(DRY_RUN, '%s java -Xmx32g -server -jar %s -s %s --pacbio-sensitive > %s' % (measure_command_wrapper(memtime_file), bin_file, reads_file, out_overlaps_file));

def run_minimap_default(reads_file, out_overlaps_file):
	memtime_file = os.path.splitext(out_overlaps_file)[0] + '.memtime';
	bin_file = '%s/minimap/minimap' % (TOOLS_PATH);
	execute_command(DRY_RUN, '%s %s %s %s > %s.paf' % (measure_command_wrapper(memtime_file), bin_file, reads_file, reads_file, out_overlaps_file));
	execute_command(DRY_RUN, '%s/miniasm/misc/paf2mhap.pl %s %s.paf > %s' % (TOOLS_PATH, reads_file, out_overlaps_file, out_overlaps_file));

def run_minimap_github_params(reads_file, out_overlaps_file):
	memtime_file = os.path.splitext(out_overlaps_file)[0] + '.memtime';
	bin_file = '%s/minimap/minimap' % (TOOLS_PATH);
	execute_command(DRY_RUN, '%s %s -Sw5 -L100 -m0 %s %s > %s.paf' % (measure_command_wrapper(memtime_file), bin_file, reads_file, reads_file, out_overlaps_file));
	execute_command(DRY_RUN, '%s/miniasm/misc/paf2mhap.pl %s %s.paf > %s' % (TOOLS_PATH, reads_file, out_overlaps_file, out_overlaps_file));

#################################################
#################################################
#################################################

def evaluate_overlaps(overlaps_file, truth_overlaps):
	bin_file = '%s/MHAP/target/mhap-1.6.jar edu.umd.marbl.mhap.main.EstimateROC' % (TOOLS_PATH);
	execute_command(DRY_RUN, 'java -cp %s %s %s %s 2>&1 | tee %s.eval.txt' % (bin_file, truth_overlaps, overlaps_file, reads_file, overlaps_file));

def parse_results(summary_file, memtime_file):
	[sensitivity, specificity, ppv] = parse_summary(summary_file);
	[cmdline, realtime, cputime, usertime, systemtime, maxrss, time_unit, mem_unit] = parse_memtime(memtime_file);
	return [sensitivity, specificity, ppv, cputime, maxrss];

def parse_summary(summary_file):
	fp_summary = None;
	try:
		fp_summary = open(summary_file, 'r');
	except:
		sys.stderr.write('ERROR: Could not open file "%s" for reading!\n' % (summary_file));
		return;

	sensitivity = '-';
	specificity = '-';
	ppv = '-';

	for line in fp_summary:
		split_line = line.strip().split(':');
		if (split_line[0] == 'Estimated sensitivity'):
			sensitivity = float(split_line[-1].strip());
		elif (split_line[0] == 'Estimated specificity'):
			specificity = float(split_line[-1].strip());
		elif (split_line[0] == 'Estimated PPV'):
			ppv = float(split_line[-1].strip());

	return [sensitivity, specificity, ppv];

def parse_memtime(memtime_path):
	cmdline = '';
	realtime = 0;
	cputime = 0;
	usertime = 0;
	systemtime = 0;
	maxrss = 0;
	rsscache = 0;
	time_unit = '';
	mem_unit = '';

	try:
		fp = open(memtime_path, 'r');
		lines = [line.strip() for line in fp.readlines() if (len(line.strip()) > 0)];
		fp.close();
	except Exception, e:
		# sys.stderr.write('ERROR: Could not find memory and time statistics in file "%s".\n' % (memtime_path));
		return [cmdline, realtime, cputime, usertime, systemtime, maxrss, time_unit, mem_unit];

	for line in lines:
		if (len(line.strip().split(':')) != 2):
			continue;

		if (line.startswith('Command line:')):
			cmdline = line.split(':')[1].strip();
		elif (line.startswith('Real time:')):
			split_line = line.split(':')[1].strip().split(' ');
			realtime = float(split_line[0].strip());
			try:
				time_unit = split_line[1].strip();
			except:
				sys.stderr.write('ERROR: No time unit specified in line: "%s"!\n' % (line.strip()));
		elif (line.startswith('CPU time:')):
			split_line = line.split(':')[1].strip().split(' ');
			cputime = float(split_line[0].strip());
			try:
				time_unit = split_line[1].strip();
			except:
				sys.stderr.write('ERROR: No time unit specified in line: "%s"!\n' % (line.strip()));
		elif (line.startswith('User time:')):
			split_line = line.split(':')[1].strip().split(' ');
			usertime = float(split_line[0].strip());
			try:
				time_unit = split_line[1].strip();
			except:
				sys.stderr.write('ERROR: No time unit specified in line: "%s"!\n' % (line.strip()));
		elif (line.startswith('System time:')):
			split_line = line.split(':')[1].strip().split(' ');
			systemtime = float(split_line[0].strip());
			try:
				time_unit = split_line[1].strip();
			except:
				sys.stderr.write('ERROR: No time unit specified in line: "%s"!\n' % (line.strip()));
		elif (line.startswith('Maximum RSS:')):
			split_line = line.split(':')[1].strip().split(' ');
			maxrss = float(split_line[0].strip());
			try:
				mem_unit = split_line[1].strip();
			except:
				sys.stderr.write('ERROR: No time unit specified in line: "%s"!\n' % (line.strip()));

	if (cputime == -1.0):
		cputime = usertime + systemtime;

	return [cmdline, realtime, cputime, usertime, systemtime, maxrss, time_unit, mem_unit];



#################################################
#################################################
#################################################

def run_overlap(reads_file, truths_file, output_path):
	if (not os.path.exists(output_path)):
		os.makedirs(output_path);

	run_graphmap(reads_file, '%s/overlaps-graphmap.mhap' % (output_path));
	run_mhap_default(reads_file, '%s/overlaps-mhap-default.mhap' % (output_path));
	run_mhap_nanopore_fast(reads_file, '%s/overlaps-mhap-nanopore_fast.mhap' % (output_path));
	run_mhap_pacbio_fast(reads_file, '%s/overlaps-mhap-pacbio_fast.mhap' % (output_path));
	run_mhap_pacbio_sensitive(reads_file, '%s/overlaps-mhap-pacbio_sensitive.mhap' % (output_path));
	run_minimap_default(reads_file, '%s/overlaps-minimap-default.mhap' % (output_path));
	run_minimap_github_params(reads_file, '%s/overlaps-minimap-github_params.mhap' % (output_path));

	evaluate_overlaps('%s/overlaps-graphmap.mhap' % (output_path), truths_file);
	evaluate_overlaps('%s/overlaps-mhap-default.mhap' % (output_path), truths_file);
	evaluate_overlaps('%s/overlaps-mhap-nanopore_fast.mhap' % (output_path), truths_file);
	evaluate_overlaps('%s/overlaps-mhap-pacbio_fast.mhap' % (output_path), truths_file);
	evaluate_overlaps('%s/overlaps-mhap-pacbio_sensitive.mhap' % (output_path), truths_file);
	evaluate_overlaps('%s/overlaps-minimap-default.mhap' % (output_path), truths_file);
	evaluate_overlaps('%s/overlaps-minimap-github_params.mhap' % (output_path), truths_file);

	results = '';
	results += '%s\n' % (output_path);
	results += 'overlapper\tsensitivity\tspecificity\tppv\tcputime\tmaxrss\n';
	results += 'graphmap\t' + '\t'.join([str(val) for val in parse_results('%s/overlaps-graphmap.mhap.eval.txt' % (output_path), '%s/overlaps-graphmap.memtime' % (output_path))]) + '\n';
	results += 'mhap-nanopore_fast\t' + '\t'.join([str(val) for val in parse_results('%s/overlaps-mhap-nanopore_fast.mhap.eval.txt' % (output_path), '%s/overlaps-mhap-nanopore_fast.memtime' % (output_path))]) + '\n';
	results += 'mhap-pacbio_fast\t' + '\t'.join([str(val) for val in parse_results('%s/overlaps-mhap-pacbio_fast.mhap.eval.txt' % (output_path), '%s/overlaps-mhap-pacbio_fast.memtime' % (output_path))]) + '\n';
	results += 'mhap-pacbio_sensitive\t' + '\t'.join([str(val) for val in parse_results('%s/overlaps-mhap-pacbio_sensitive.mhap.eval.txt' % (output_path), '%s/overlaps-mhap-pacbio_sensitive.memtime' % (output_path))]) + '\n';
	results += 'mhap-default\t' + '\t'.join([str(val) for val in parse_results('%s/overlaps-mhap-default.mhap.eval.txt' % (output_path), '%s/overlaps-mhap-default.memtime' % (output_path))]) + '\n';
	results += 'minimap-default\t' + '\t'.join([str(val) for val in parse_results('%s/overlaps-minimap-default.mhap.eval.txt' % (output_path), '%s/overlaps-minimap-default.memtime' % (output_path))]) + '\n';
	results += 'minimap-github_params\t' + '\t'.join([str(val) for val in parse_results('%s/overlaps-minimap-github_params.mhap.eval.txt' % (output_path), '%s/overlaps-minimap-github_params.memtime' % (output_path))]) + '\n';
	results += '\n';
	fp = open('%s/summary.csv' % (output_path), 'w');
	fp.write(results);
	fp.close();



#################################################
#################################################
#################################################

def download_and_install():
	if (not os.path.exists(TOOLS_PATH)):
		os.makedirs(TOOLS_PATH);

	setup_graphmap();
	setup_mhap();
	setup_daligner();
	setup_minimap();

def verbose_usage_and_exit():
	sys.stderr.write('Runs overlapping on simulated (or real) datasets.\n\n');
	sys.stderr.write('Usage:\n');
	sys.stderr.write('\t%s mode [<reads.fasta> <truth_overlaps.m4> <output_path>]\n' % sys.argv[0]);
	sys.stderr.write('\n');
	sys.stderr.write('\t- mode - either "run" or "setup". If "setup" other parameters can be ommitted.\n');
	sys.stderr.write('\n');
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

		reads_file = sys.argv[2];
		truths_file = sys.argv[3];
		output_path = sys.argv[4];

		run_overlap(reads_file, truths_file, output_path);

	else:
		verbose_usage_and_exit();
