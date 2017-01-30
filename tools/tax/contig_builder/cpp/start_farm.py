#!/usr/bin/python
import subprocess
import os
import stat
import glob
import sys
import shell

def main():
	if __name__ != "__main__":
		return

	if len(sys.argv) < 2:
		print "need <accession> [<skip reads file>]"
		return

#	align_it = False
	acc = sys.argv[1]
#	spots = get_total_spots(acc)
#	if spots > 70*1000*1000:
#		print "run is too large -", spots, "spots"
#		return

	skip_reads_file = ""
	if len(sys.argv) >= 3:
		skip_reads_file = sys.argv[2] 

	script_dir = os.path.dirname(os.path.realpath(__file__))
	cmdline = "qsub " + qsub_options(acc, script_dir, skip_reads_file) + " " + script_dir + "/start_farm.sh"
	print cmdline
	print shell.execute_local(cmdline)

def qsub_options(acc, script_dir, skip_reads_file):
	spots = get_total_spots(acc)
	print "run spots", spots

	million = 1000000

	mem = 20  # todo: tune memory params

	if spots > 2 * million:
		mem = 40

	if spots > 15 * million:
		mem = 100

	if spots > 60 * million:
		mem = 110

	job_name = acc

	minute = 60
	hour = 60*minute
	max_task_life_time_seconds = 2*hour
#	max_task_life_time_seconds = 16*minute
#	options = "-terse -sync n -cwd -P unified -j n -N " + job_name + " -b n -l mem_free=" + str(mem) + "G,h_vmem=" + str(mem) + "G,h_rt=" + str(max_task_life_time_seconds) +" -v accession=" + acc +",script_dir=" + script_dir
	#todo: tune reserve_mem and mem_free
	options = "-terse -sync n -cwd -P unified -j n -N " + job_name + " -b n -l mem_free=" + str(mem) + "G,h_vmem=" + str(mem) + "G,reserve_mem=" + str(mem) + "G,h_rt=" + str(max_task_life_time_seconds) +" -v accession=" + acc +",script_dir=" + script_dir + ",skip_reads=" + skip_reads_file
#	options += " -m n" #no e-mail # todo: uncomment
	return options

def get_total_spots(acc):
	out = shell.execute_local("/net/snowman/vol/projects/trace_software/vdb/linux/release/x86_64/bin/sra-stat -x --quick " + acc)
	out = out.split('spot_count="')[1].split('"')[0]
	return int(out)

main()
