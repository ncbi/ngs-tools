#!/opt/python-2.7/bin/python
import sys
import argparse
import subprocess

def get_analysis_results(hours): #server, db_name, user, password):
#	separator = '@'
	p = subprocess.Popen(
		["sqsh-ms -S SRA_PRIMARY -D SRA_Main -U sra_sa -P sra_sa_pw -b -h "], #-s " + separator], 
		stdout = subprocess.PIPE, 
		stdin = subprocess.PIPE, 
		shell = True)

	input = "select spots, bases, DATEDIFF(minute, load_time, GETDATE()) from SRA_Track..LoadResults(nolock) lr join SRA_Main..SRAFiles(nolock) sf  on lr.link_uid = sf.file_id and lr.link_type = 'SRAFile' join SRA_Main..Run(nolock) r on r.acc = sf.acc where job_type = 'tax_analysis' and ok = 1 and load_time > DATEADD(HOUR, " + str(-hours) + ", GETDATE())"
	input += "\ngo\nexit\n"
	out, err = p.communicate(input)
	if err:
		raise Exception("fetching sql data error: " + err)

	out = out.split('\n')
	out = [line.split() for line in out if line != ""]
	out = [(int(spots), int(bases), int(min_ago)) for (spots, bases, min_ago) in out if spots != "NULL" and bases != "NULL"]
	return out

def group_results(r):
	grouped = dict()
	GROUP_BUCKET = 60
	for spots, bases, min_ago in r:
		group = min_ago/GROUP_BUCKET # todo: think
		if min_ago > 0 and min_ago % GROUP_BUCKET == 0:
			group -= 1

		if not group in grouped:
			grouped[group] = (0, 0, 0)

		runs, sum_spots, sum_bases = grouped[group]
		grouped[group] = (runs + 1, sum_spots + spots, sum_bases + bases)

#	grouped = [grouped[group] for group in grouped]
	return grouped

def print_perf(r):
	print "hrs ago\truns\tMspots\tGbases"
	for h in xrange(group_count(r)):
		runs, sum_spots, sum_bases = r[h]
		print str(h) + '\t' + str(runs) + '\t' + str(sum_spots/1000000) + '\t' + str(sum_bases/1000000000)

def y_match(v, y, top, height):
	return (v * height / top) == y # todo: check for round errors

def draw_graph(name, v):
	v = v[::-1]
	top = max(v)
	HEIGHT = 20
	print name
#	print top
	SCREEN_WIDTH = 220
	fit_times = 1
	if len(v) > 0:
		fit_times = SCREEN_WIDTH/len(v)

	skip_len = max([1, min([3, fit_times])]) - 1

	for i in xrange(HEIGHT + 1):
		y = HEIGHT - i
		if i % 4 == 0:
			print top*y/HEIGHT, '\t|',
		else:
			print '\t|',

		for x in xrange(len(v)):
			if y_match(v[x], y, top, HEIGHT):
				sys.stdout.write('*')
			else:
				sys.stdout.write(' ')

			for t in xrange(skip_len):
				sys.stdout.write(' ')

		print ""

	print '\t-',
	for x in xrange(len(v) * (skip_len + 1)):
		sys.stdout.write('-')

	print ""

def group_count(r):
	return max( [x for x in r] ) + 1

def to_array(r, part):
	v = []
	for h in xrange(group_count(r)):
		if not h in r:
			v.append(0)
		else:
			v.append(r[h][part])

	return v

def draw_perf(r):
	draw_graph("runs", to_array(r, 0))
	spots = to_array(r, 1)
	spots = [s/1000000 for s in spots]
	draw_graph("megaspots", spots)
	bases = to_array(r, 2)
	bases = [b/1000000000 for b in bases]
	draw_graph("gigabases", bases)

def main():
	parser = argparse.ArgumentParser(description='tax analysis report generator')
	parser.add_argument('-b', '--hours-back', metavar='N', type=int, help='hours back', default = 12)
	parser.add_argument('--draw', action='store_true', help='draw graph')

	args = parser.parse_args()
	analysis_results = get_analysis_results(args.hours_back)
#	print analysis_results
#	return
	grouped_results = group_results(analysis_results)
#	print grouped_results
#	return

	if args.draw:
		draw_perf(grouped_results)
	else:
		print_perf(grouped_results)

if __name__ == "__main__":
	main()
