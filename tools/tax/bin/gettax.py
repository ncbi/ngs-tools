#!/opt/python-2.7env/bin/python

import os
import tempfile
import tarfile
import shutil
import subprocess
import sqlite3
import logging
import collections
import errno
import argparse
import time
from datetime import datetime

logger = logging.getLogger('gettax')

PROGRESS_GRANULARITY = 100000

Taxon = collections.namedtuple('Taxon', ['tax_id', 'parent_tax_id', 'rank', 'scientific_name'])

class Connection(object):
    def __init__(self, sqlite_conn):
        self._sqlite_conn = sqlite_conn
        self.close = self._sqlite_conn.close

    def __enter__(self):
        self._sqlite_conn.__enter__()
        return self

    def __exit__(self, *args):
        return self._sqlite_conn.__exit__(*args)

    def gettax(self, tax_id):
        return _gettax(tax_id, self._sqlite_conn)

def _parse_line(line):
    """parses original dump line"""

    assert line.endswith('\t|\n'), repr(line)
    line = line[:-3]
    parts = line.split('\t|\t')
    return parts

def _parse_nodes(conn, path):
    """parses nodes dump"""

    count = 0
    with open(path) as f:
        for line in f:
            parts = _parse_line(line)
            assert len(parts) >= 3, repr(line)
            tax_id, parent_tax_id, rank = parts[:3]
            tax_id = int(tax_id)
            parent_tax_id = int(parent_tax_id)
            if rank == 'no rank' or not rank:
                rank = None
            conn.execute('insert into taxons values (?, ?, ?, null)', [tax_id, parent_tax_id, rank])
            count += 1
            if count % PROGRESS_GRANULARITY == 0:
                logger.debug('parsed %s nodes', count)
    return count

def _parse_names(conn, path):
    """parses names dump"""

    count = 0
    cur = conn.cursor()
    with open(path) as f:
        for line in f:
            parts = _parse_line(line)
            assert len(parts) >= 4, repr(line)
            tax_id, name, _, name_class = parts[:4]
            tax_id = int(tax_id)
            assert name, repr(line)
            assert name_class, repr(line)
            if name_class == 'scientific name':
                cur.execute('update taxons set scientific_name = ? where tax_id = ?', [name, tax_id])
                assert cur.rowcount == 1
                count += 1
                if count % PROGRESS_GRANULARITY == 0:
                    logger.debug('parsed %s names', count)
    cur.close()
    return count

def _parse_merged(conn, path):
    count = 0
    #cur = conn.cursor()
    with open(path) as f:
        for line in f:
            parts = _parse_line(line)
            assert len(parts) == 2, repr(line)
            old_id, new_id = parts
            conn.execute('''
insert into taxons
select ?, parent_tax_id, rank, scientific_name
from taxons
where tax_id = ?
            ''', [old_id, new_id])
            count += 1
    return count

def create_db(src, dst):
    """converts dump into sqlite database"""

    logger.debug('creating sqlite db: %s -> %s', src, dst)
    assert os.path.exists(src), 'dump does not exist: %s' % src

    md5 = src + '.md5'
    if os.path.exists(md5):
        logger.info('md5 check')
        src_dir = os.path.dirname(src)
        subprocess.check_output(['md5sum', '-c', md5], cwd=src_dir)

    before = datetime.now()
    conn = sqlite3.connect(dst, isolation_level='EXCLUSIVE')
    temp_dir = tempfile.mkdtemp()
    try:
        conn.execute('''
create table if not exists taxons (
    tax_id int primary key,
    parent_tax_id int,
    rank text,
    scientific_name text
    )
        ''')
        conn.execute('''
create table if not exists rebuilds (
        timestamp datetime,
        duration real
    )
        ''')
        conn.execute('delete from taxons')

        logger.info('extracting')
        tar = tarfile.open(src)
        tar.extractall(temp_dir)

        logger.info('parsing nodes')
        nodes = os.path.join(temp_dir, 'nodes.dmp')
        node_count = _parse_nodes(conn, nodes)
        logger.info('parsed %s nodes', node_count)

        logger.info('parsing names')
        names = os.path.join(temp_dir, 'names.dmp')
        name_count = _parse_names(conn, names)
        logger.info('parsed %s names', name_count)

        logger.info('sanity check')
        cur = conn.cursor()
        cur.execute('select count(*) from taxons where scientific_name is null')
        rows = cur.fetchall()
        assert len(rows) == 1 and len(rows[0]) == 1
        assert rows[0][0] == 0, '%s taxons without name' % rows[0][0]
        assert node_count == name_count, 'node/name count mismatch: %s != %s' % (node_count, name_count)
        cur.close()

        logger.info('parsing merged')
        merged = os.path.join(temp_dir, 'merged.dmp')
        merged_count = _parse_merged(conn, merged)
        logger.info('parsed %s merged ids', merged_count)

        logger.info('fixing root')
        conn.execute('update taxons set parent_tax_id = 0 where tax_id = parent_tax_id')

        logger.info('writing metainfo')
        duration = datetime.now() - before
        conn.execute('insert into rebuilds values (?, ?)', [before.isoformat(' '), duration.total_seconds()])

        logger.info('committing')
        conn.commit()

        logger.info('done')
    finally:
        shutil.rmtree(temp_dir)
        conn.close()

def _mtime(path):
    """returns path mtime or none if path does not exist"""

    try:
        return os.path.getmtime(path)
    except OSError as e:
        if e.errno != errno.ENOENT:
            raise
        return None

def _gettax(tax_id, conn):
    """returns corresponding Taxon instance by tax_id"""
    tax_id = int(tax_id)

    cur = conn.cursor()
    cur.execute('select * from taxons where tax_id = ?', [tax_id])
    rows = cur.fetchall()
    cur.close()

    if rows:
        assert len(rows) == 1
        return Taxon(*rows[0])

def connect(taxdump, sqlite_cache=None, rebuild_timeout=None, connection_timeout=None):
    if sqlite_cache is None:
        assert taxdump, 'netiher taxdump nor cache was provided'
        sqlite_cache_file = tempfile.NamedTemporaryFile(prefix='taxdump', suffix='.sqlite')
        sqlite_cache = sqlite_cache_file.name
        cache_mtime = None
    else:
        cache_mtime = _mtime(sqlite_cache)

    if taxdump:
        rebuild = False
        dump_mtime = _mtime(taxdump)
        if cache_mtime is None:
            logger.debug('rebuilding cache because it does not exist')
            rebuild = True
        elif cache_mtime < dump_mtime:
            time_since_last_rebuild = time.time() - cache_mtime

            if rebuild_timeout is None or time_since_last_rebuild > rebuild_timeout:
                logger.debug('rebuilding cache because it is out of date')
                rebuild = True
        if rebuild:
            try:
                create_db(taxdump, sqlite_cache)
            except sqlite3.OperationalError as e:
                if 'database is locked' not in str(e):
                    raise
                logger.info('database is already being rebuilt, skipping')

    assert os.path.exists(sqlite_cache), 'cache does not exist: %s' % sqlite_cache
    logger.debug('Using cache: %s', sqlite_cache)
    kwargs = {}
    if connection_timeout is not None:
        kwargs = dict(timeout=connection_timeout)
    conn = sqlite3.connect(sqlite_cache, **kwargs)
    return Connection(conn)

def format_taxon(taxon):
    """creates output similar to original gettax output"""

    if not taxon:
        return 'tax id not found'

    attrs = taxon._fields
    taxon = taxon._asdict()
    taxon['rank'] = taxon['rank'] or 'no rank'
    lines = []
    longest_attr = max(attrs, key=len)
    for attr in attrs:
        value = taxon[attr]
        attr_view = attr.replace('_', ' ').rjust(len(longest_attr))
        line = '%s: %s' % (attr_view, value)
        lines.append(line)
    return '\n'.join(lines)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--tax-dump')
    parser.add_argument('-c', '--sqlite-cache')
    parser.add_argument('-f', '--force-rebuild', action='store_true')
    parser.add_argument('-r', '--rebuild-timeout', type=float, help='delay between rebuilds in seconds')
    parser.add_argument('--connection-timeout', type=float, help='connection timeout, to wait until rebuilt is completed')
    parser.add_argument('-d', '--debug', action='store_true', help='enable debug logging')
    parser.add_argument('-q', '--quiet', action='store_true', help='disable info logging')
    parser.add_argument('-a', '--attribute', help='output only the given attribute')
    parser.add_argument('tax_id', metavar='TAX_ID', nargs='*', type=int, help='list of tax ids')

    args = parser.parse_args()
    if args.attribute:
        args.attribute = args.attribute.strip().replace(' ', '_')
    if not args.tax_dump and not args.sqlite_cache:
        parser.error('need tax_dump or sqlite_cache or both to work')
    if args.debug and args.quiet:
        parser.error('debug and quiet options are incompatible')
    if args.force_rebuild and not (args.tax_dump and args.sqlite_cache):
        parser.error('force rebuild needs both taxdump and cache paths')
    if not args.force_rebuild and not args.tax_id:
        parser.error('nothing to do')

    level = logging.INFO
    if args.debug:
        level = logging.DEBUG
    elif args.quiet:
        level = logging.WARNING
    logging.basicConfig(level=level)

    if args.force_rebuild:
        create_db(args.tax_dump, args.sqlite_cache)

    if args.tax_id:
        with connect(args.tax_dump, args.sqlite_cache, args.rebuild_timeout, args.connection_timeout) as conn:
            for tax_id in args.tax_id:
                taxon = conn.gettax(tax_id)
                if args.attribute:
                    assert args.attribute in taxon._fields, 'unknown attribute: %s' % args.attribute
                    print taxon._asdict()[args.attribute]
                else:
                    print
                    print format_taxon(taxon)
    else:
        logger.debug('no tax_ids to process')
