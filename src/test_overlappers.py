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
	execute_command(DRY_RUN, '%s %s -r %s -d %s -o %s -L mhap' % (measure_command_wrapper(memtime_file), bin_file, reads_file, reads_file, out_overlaps_file));

def run_mhap(reads_file, out_overlaps_file):
	memtime_file = os.path.splitext(out_overlaps_file)[0] + '.memtime';
	bin_file = '%s/MHAP/target/mhap-1.6.jar' % (TOOLS_PATH);
	execute_command(DRY_RUN, '%s java -Xmx32g -server -jar %s -s %s > %s' % (measure_command_wrapper(memtime_file), bin_file, reads_file, out_overlaps_file));

def run_minimap(reads_file, out_overlaps_file):
	memtime_file = os.path.splitext(out_overlaps_file)[0] + '.memtime';
	bin_file = '%s/minimap/minimap' % (TOOLS_PATH);
	execute_command(DRY_RUN, '%s %s -Sw5 -L100 -m0 %s %s > %s.paf' % (measure_command_wrapper(memtime_file), bin_file, reads_file, out_overlaps_file));
	execute_command(DRY_RUN, '%s/miniasm/misc/paf2mhap.pl %s %s.paf > %s' % (reads_file, out_overlaps_file, out_overlaps_file));

def evaluate_overlaps(overlaps_file, truth_overlaps):
	bin_file = '%s/MHAP/target/mhap-1.6.jar edu.umd.marbl.mhap.main.EstimateROC' % (TOOLS_PATH);
	execute_command(DRY_RUN, 'java -cp %s %s %s' % (bin_file, overlaps_file, truth_overlaps));

def run_overlap(reads_file, truths_file, output_path):
	run_graphmap(reads_file, '%s/overlaps-graphmap.mhap' % (output_path));
	run_mhap(reads_file, '%s/overlaps-mhap.mhap' % (output_path));
	run_minimap(reads_file, '%s/overlaps-minimap.mhap' % (output_path));

	evaluate_overlaps('%s/overlaps-graphmap.mhap' % (output_path), truths_file);
	evaluate_overlaps('%s/overlaps-mhap.mhap' % (output_path), truths_file);
	evaluate_overlaps('%s/overlaps-minimap.mhap' % (output_path), truths_file);

#################################################
#################################################
#################################################

def download_and_install():
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
