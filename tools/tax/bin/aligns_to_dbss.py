#!/bin/env python2.7

import os
import sys
import argparse
import subprocess
import tempfile
import threading
import logging
import re
import functools
import collections
import datetime
import signal
import time

logger = logging.getLogger('aligns_to_dbss')
TIMESTAMP_FORMAT = '%Y-%m-%d %H:%M:%S'

class GlueError(Exception):
    def __init__(self, message, rc=3):
        Exception.__init__(self, message)
        self.rc = rc

def exit_on_exception(fn):
    @functools.wraps(fn)
    def wrapped(*args, **kwargs):
        try:
            return fn(*args, **kwargs)
        except GlueError as e:
            logger.error('%s', e)
            os._exit(e.rc)
        except Exception as e:
            logger.exception('Unhandled exception')
            os._exit(5)
    return wrapped

SIGNAL_NAMES = dict((number, name) for name, number in signal.__dict__.items()
                    if name.startswith('SIG') and not name.startswith('SIG_'))

def format_returncode(returncode):
    if returncode < 0:
        return '%s (rc=%s)' % (SIGNAL_NAMES.get(-returncode, 'unknown signal'), returncode)
    else:
        return 'rc=%s' % returncode

def extract_tax_list(input_file, output_file):
    spot_count = 0
    last_spot_id = None
    tax_set = collections.Counter()
    
    for line in input_file:
        parts = line.split('\t')
        spot_id = parts[0]
        assert spot_id >= last_spot_id, 'aligns_to output is not sorted'
        if spot_id != last_spot_id:
            last_spot_id = spot_id
            spot_count += 1
        tax_ids = parts[1:]
        tax_set.update(tax_ids)

    only_one_hit = 0
    for tax_id in sorted(tax_set):
        output_file.write('%s\n' % tax_id)
        if tax_set[tax_id] == 1:
            only_one_hit += 1
    output_file.flush()
    logger.info('Extracted %s tax ids from %s identified spots (%s tax ids with only 1 hit)', len(tax_set), spot_count, only_one_hit)
    
@exit_on_exception
def stderr_processor(process_stderr, first_pass):
    last_progress = None
    while True: # for line in file does not work because of broken buffernig
        line = process_stderr.readline()
        if not line:
            break
        # strip caret returns, they break matching and pass1 / pass2 prefixes
        while line.startswith('\r'):
            line = line[1:]
        m = re.match(r'^(\d+)% processed$', line)
        if m: # report total progress
            percent = int(m.group(1))
            if first_pass:
                percent /= 2
            else:
                percent = 50 + percent / 2
            if last_progress != percent:
                last_progress = percent
                logger.info('%s%% processed', percent)
        else: # pass prefix insertion
            timestamp = datetime.datetime.now().strftime(TIMESTAMP_FORMAT)
            sys.stderr.write(timestamp)
            if first_pass:
                sys.stderr.write(' pass1: ')
            else:
                sys.stderr.write(' pass2: ')
            sys.stderr.write(line)
    sys.stderr.flush()

@exit_on_exception
def process(path, dbs, dbss, unaligned_only, first_step_output, tax_list):
    self_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
    aligns_to = os.path.join(self_dir, 'aligns_to')
    for f in [dbs, dbss]:
        if not os.path.exists(f):
            raise GlueError('%s does not exist' % f)

    logger.info('Starting first pass')
    cmdline = [aligns_to, '-dbs', dbs, path]
    if unaligned_only:
        cmdline.append('-unaligned_only')
    logger.debug('Cmdline: %s', ' '.join(cmdline))
    process = subprocess.Popen(cmdline, stdout=first_step_output, stderr=subprocess.PIPE)
    thread_stderr = threading.Thread(target=stderr_processor, name='stderr_processor', args=(process.stderr, True))
    thread_stderr.start()
    thread_stderr.join()
    process.wait()
    first_step_output.flush()
    logger.info('First pass done')

    subprocess.check_call(['sort', first_step_output.name, '-o', first_step_output.name])
    logger.info('First pass output sorted')

    first_step_output.seek(0)
    extract_tax_list(first_step_output, tax_list)
    tax_list.flush()
    logger.info('Tax list extracted')
    
    if process.returncode != 0:
        msg = 'first pass finished with %s' % format_returncode(process.returncode)
        if process.returncode < 0: # forward signals to self
            logger.info('%s, forwarding to self', msg)
            os.kill(os.getpid(), -process.returncode)
            time.sleep(0.1)
        raise GlueError(msg, 1)
    if tax_list.tell() == 0:
        logger.info('Tax list is empty, no point in doing second step, quitting')
        return

    logger.info('Starting second pass')
    cmdline = [aligns_to, '-dbss', dbss, '-tax_list', tax_list.name, path]
    if unaligned_only:
        cmdline.append('-unaligned_only')
    logger.debug('Cmdline: %s', ' '.join(cmdline))
    process = subprocess.Popen(cmdline, stderr=subprocess.PIPE)
    thread_stderr = threading.Thread(target=stderr_processor, name='stderr_processor', args=(process.stderr, False))
    thread_stderr.start()
    thread_stderr.join()
    process.wait()
    logger.info('Second pass done')
    if process.returncode != 0:
        msg = 'second pass finished with %s' % format_returncode(process.returncode)
        if process.returncode < 0: # forward signals to self
            logger.info('%s, forwarding to self', msg)
            os.kill(os.getpid(), -process.returncode)
            time.sleep(0.1)
        raise GlueError(msg, 2)

if __name__ == '__main__':
    # logging.basicConfig(level=logging.DEBUG)
    # extract_tax_list(sys.stdin, sys.stdout)
    # sys.exit(0)
    parser = argparse.ArgumentParser(description='two step alignment helper')
    parser.add_argument('-v', action='store_true', help='verbose output')
    parser.add_argument('-unaligned_only', action='store_true', help='verbose output')
    parser.add_argument('-dbs', required=True, help='path to dbs database')
    parser.add_argument('-dbss', required=True, help='path to dbss database')
    parser.add_argument('-first_step_output', help='first step output location, if not presented, temporary file is used')
    parser.add_argument('-tax_list', help='tax list output location, if not presented, temporary file is used')
    parser.add_argument('path', metavar='PATH', help='path to file to process')
    args = parser.parse_args()
    
    level = logging.DEBUG if args.v else logging.INFO
    logging.basicConfig(level=level, format='%(asctime)s glue : %(message)s', datefmt=TIMESTAMP_FORMAT)
    
    if args.tax_list:
        tax_list = open(args.tax_list, 'w+b')
    else:
        tax_list = tempfile.NamedTemporaryFile()
    if args.first_step_output:
        first_step_output = open(args.first_step_output, 'w+b')
    else:
        first_step_output = tempfile.NamedTemporaryFile()
        
    process(args.path, args.dbs, args.dbss, args.unaligned_only, first_step_output, tax_list)
