#!/usr/bin/python
import subprocess
import perf
import os
import stat
import glob

def execute_async_result(p, check_error = True):
	out, err = p.communicate()
	
	if check_error and len(err) > 0:
#		print err
		raise Exception(err)

	return out

def execute_async_start(cmdline):
#	print cmdline
	return subprocess.Popen([cmdline], stdout = subprocess.PIPE, stderr = subprocess.PIPE, stdin = subprocess.PIPE, shell = True) # todo: try w/o shell

def execute_local(cmdline, check_error = True):
	return execute_async_result(execute_async_start(cmdline), check_error)
 
#	return execute_async_result(p, check_error, None)

