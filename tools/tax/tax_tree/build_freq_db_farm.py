#!/usr/bin/python
import sys
import subprocess
import perf
import os
import stat
import glob

def execute_async_result(p, check_error):
	out, err = p.communicate()
#	print out
#	if sh_filename:
#		os.remove(sh_filename)
#		qsub_cleanup(sh_filename)
	
	if check_error and len(err) > 0:
		print err
		raise Exception(err)

	return out

def execute_local(cmdline, check_error = True):
	print "exec_local cmdline:", cmdline
	p = subprocess.Popen([cmdline], stdout = subprocess.PIPE, stderr = subprocess.PIPE, stdin = subprocess.PIPE, shell = True) # todo: try w/o shell
	return execute_async_result(p, check_error)

def load_files_list(filename):
	files = []
	f = open(filename)
	for line in f:
		name = line.rstrip().split('\t')[-1]
		files.append(name)

	return files

def qsub_options(task_name):
	minute = 60
	hour = 60*minute
	max_task_life_time_seconds = 24*hour
	options = "-sync n -cwd -P unified -j n -N " + task_name + " -m n -b n -l mem_free=4G,h_vmem=8G,h_rt=" + str(max_task_life_time_seconds)
#	print options
#	options += " -m n" #no e-mail # todo: uncomment
	return options

def grant_exec(filename):
	os.chmod(filename, os.stat(filename).st_mode | stat.S_IEXEC)

def get_task_name(fasta):
	tax_id = int(fasta.split(".")[-2].split('/')[-1])
	return "build_freq_db_" + str(tax_id)

def quoted(filename):
	return "'" + filename.replace("'", "'\\''") + "'"
#	return "'" + filename + "'"

def start_processing(fasta):
	print fasta
	shell_file = fasta + ".sh"
	f = open(shell_file, 'w')
	path = os.path.dirname(os.path.realpath(__file__))
	f.write(path + "/search/build_freq_db " + quoted(fasta) + " " + quoted(fasta + ".freq9.amino") + " 9 18")
	f.close()
	grant_exec(shell_file)
	execute_local("qsub " + qsub_options(get_task_name(fasta)) + " " + quoted(shell_file))

def main():
	if __name__ != "__main__":
		return

	if len(sys.argv) < 2:
		print "need <files.list>"
		return

	files_list = load_files_list(sys.argv[1])

	for fasta in files_list:
		start_processing(fasta)
	
main()

