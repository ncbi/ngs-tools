#!/usr/bin/python
import subprocess
import perf
import os
import stat
import glob

#@perf.timed
def execute(cmdline, check_error, used_out_file, farm):
	p, sh_filename = execute_async_start(cmdline, used_out_file, farm)
	return execute_async_result(p, check_error, sh_filename)

def execute_async_start(cmdline, used_out_file, farm):
	used_out_file = os.path.abspath(used_out_file)
	work_dir = "./" #folder_of(used_out_file)
	sh_full_filename = get_sh_full_filename(used_out_file)
#	print cmdline

	open(sh_full_filename, 'w').write(cmdline + '\n')
	grant_exec(sh_full_filename)
	cmdline = sh_full_filename
	if farm:
		cmdline = "qsub " + qsub_options(work_dir, name_without_folder(sh_full_filename)) + " " + sh_full_filename

	p = subprocess.Popen([cmdline], stdout = subprocess.PIPE, stderr = subprocess.PIPE, stdin = subprocess.PIPE, cwd = work_dir, shell = True) # todo: try w/o shell
	return p, sh_full_filename

def execute_async_still_running(p):
	return p.poll() == None

def execute_async_result(p, check_error, sh_filename):
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
#	print "exec_local cmdline:", cmdline
	p = subprocess.Popen([cmdline], stdout = subprocess.PIPE, stderr = subprocess.PIPE, stdin = subprocess.PIPE, shell = True) # todo: try w/o shell
	return execute_async_result(p, check_error, None)

def grant_exec(filename):
	os.chmod(filename, os.stat(filename).st_mode | stat.S_IEXEC)

def get_sh_full_filename(used_out_file):
	return os.path.abspath(used_out_file) + ".sh"

def folder_of(filename):
	return os.path.split(filename)[0] + '/'

def name_without_folder(filename):
	return os.path.split(filename)[1]

def qsub_cleanup(filename):
	files = glob.glob(filename + "*")
	for f in files:
		if os.stat(f).st_size == 0:
			os.remove(f)

def qsub_options(work_dir, job_name):
	minute = 60
	hour = 60*minute
	max_task_life_time_seconds = 2*hour
#	print "work dir", work_dir
	options = "-terse -sync y -wd " + work_dir + " -P unified -j n -N " + job_name + " -b n -l mem_free=40G,h_vmem=100G,h_rt=" + str(max_task_life_time_seconds) # todo: tune memory params
#	options += " -m n" #no e-mail # todo: uncomment
	return options

