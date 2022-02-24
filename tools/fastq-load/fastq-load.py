#!/opt/python-3.9/bin/python
#===========================================================================
#
#                            PUBLIC DOMAIN NOTICE
#               National Center for Biotechnology Information
#
#  This software/database is a "United States Government Work" under the
#  terms of the United States Copyright Act.  It was written as part of
#  the author's official duties as a United States Government employee and
#  thus cannot be copyrighted.  This software/database is freely available
#  to the public for use. The National Library of Medicine and the U.S.
#  Government have not placed any restriction on its use or reproduction.
#
#  Although all reasonable efforts have been taken to ensure the accuracy
#  and reliability of the software and data, the NLM and the U.S.
#  Government do not and cannot warrant the performance or results that
#  may be obtained by using this software or data. The NLM and the U.S.
#  Government disclaim all warranties, express or implied, including
#  warranties of performance, merchantability or fitness for any particular
#  purpose.
#
#  Please cite the author in any work or product based on this material.
#
#===========================================================================
#
#

"""
fastq-load.py --output=<archive path> <other options> <fastq files> | general-loader

Options:

    output:         Output archive path

    offset:         For interpretation of ascii quality (offset=33 or offset=64) or
                    indicating the presence of numerical quality (offset=0).
                    (accepts PHRED_33 and PHRED_64, too)

    quality:        Same as offset (different from latf due to requirement for '=')

    readLens:       For splitting fixed length sequence/quality from a single fastq
                    Number of values must = # of reads (comma-separated). Must be
                    consistent with number of read types specified

    readTypes:      For specifying read types using B (Biological) or T (Technical)
                    Use sequence like TBBT - no commas. Defaults to BB. Must be
                    consistent with number of values specified for read lengths.
                    If you want the read sequence to be used as the spot group or
                    part of the spot group use G (Group). Multiple reads incorporated
                    into the spot group will be concatenated with an '_' separator.

    spotGroup:      Indicates spot group to associate with all spots that are loaded
                    Overrides barcodes on deflines.

    orphanReads:    File or files contain orphan reads or fragments mixed with non-orphan
                    reads. If all files are either completely fragments or completely
                    paired (e.g. split-3 fastq-dump output) then this option is
                    unnecessary. However, for pre-split 454 this option would probably be
                    necessary.

    logOdds:        Input fastq has log odds quality (set prior to 'offset' or 'quality'
                    option)

    ignAndDiscardNames:   For when names repeat, go on position only. Specify eight line fastq
                    if that is the case. Does not work with orphanReads. Does not work
                    for split seq/qual files.

    useAndDiscardNames:   Too many names to store but still useful for determination of
                    read pairs. So names are used and then discarded.

    ignoreNames:    Determination of pairs via names will not work but first read name
                    retained. Pairing based on position only.

    read1PairFiles: Filenames containing read1 seqs ordered to correspond with read2 pair
                    files. Required with --orphanReads option. Files paths must still be
                    provided on the command line in addition to this option. Must be
                    specified in conjunction with read2PairFiles option. Comma-separated.
                    Also useful if ignoring/discarding names. Put '\\' before spaces in
					file names.

    read2PairFiles: Filenames containing read2 seqs. Required for --orphanReads option.
                    Include a filename from the read1 files if read2 is in the same file
                    using corresponding positions. Files paths must still be provided on
                    the command line in addition to this option. Comma-separated. Must
                    be specified in conjunction with read1PairFiles option.
                    Also useful if ignoring/discarding names. Put '\\' before spaces in
					file names.

    read3PairFiles: Filenames containing read3 seqs. (NO ORPHANS)

    read4PairFiles: Filenames containing read4 seqs. (NO ORPHANS)

    read5PairFiles: Filenames containing read5 seqs. (NO ORPHANS)

    read6PairFiles: Filenames containing read6 seqs. (NO ORPHANS)

    read1QualFiles: Filenames containing read1 quals ordered to correspond with read1 pair
                    files. Files paths must still be provided on the command line in addition
                    to this option. Must be specified in conjunction with read2QualFiles option.
                    Comma-separated.. Put '\\' before spaces in file names.
					(Provide only if ignoring or discarding names)

    read2QualFiles: Filenames containing read2 quals ordered to correspond with read2 pair
                    files. Files paths must still be provided on the command line in addition
                    to this option. Must be specified in conjunction with read1QualFiles option.
                    Comma-separated. Put '\\' before spaces in file names.
					(Provide only if ignoring or discarding names)

    filesToIgnore:  Filenames to ignore when loading. Put '\\' before spaces in file names.
                    (comma-separated list of filenames)

    filesToTruncate:  Filenames to truncate along with the truncate line numbers.
	                  Put '\\' before spaces in file names.
                      (comma-separated list containing entries of 'Filename:LineNumber'.)

    fileSpotGroups: Comma-separated list of barcodes based on alphabetic order of files. Empty
                    values are represented by back-to-back commas. Only applies to 1st file in
                    a paired set and not to qual-only files

    useFilenameForSG: Use file name without suffixes for spot group. Only applies to 1st file in
                    a paired set and not to qual-only files

    platform:       454, Pacbio, Illumina, ABI, etc.

    readLabels:     Rev, Fwd, Barcode, etc., comma-separated (no whitespaces)
                    For reads used for spot groups, specify an empty string (e.g. ',,')

    mixedDeflines:  Indicates mixed defline types exist in one of the fastq files.
                    Results in slower processing of deflines.

    genericDeflines:Forces retention of all non-whitespace characters after defline character
                    as spot name (no spot group or read number extraction performed)

    mixedTypes:     Allow for mixed types of fastq to be processed in one run
                    (e.g. normal and multiline)

    ignLeadChars:   Set # of leading defline characters to ignore for pairing

    ignTailChars:   Set # of tailing defline characters to ignore for pairing

    ignCharsAfterPlus: Ignore characters that occur after qual defline plus character

    discardBarcodes For cases where too many barcodes exist (>30000)

    discardBCchars: Set # of leading or trailing (negative value) barcode characters to discard

    appendBCtoName: Append barcode to spot name

    concatPairFiles Indicate that provided files contain concatenated pair files.
                    Fastq format only, one pair set per file, no multiline seq/qual,
                    should use with uncompressed files like the SRA production pipeline,
                    and platform should not be nanopore. Put '\\' before spaces in file names.

    removeBadChars: Remove control/unprintable characters found in fastq lines to allow
                    processing the fastq files (e.g. nulls, extra linefeeds from windows)

    removeSeqSpaces Remove space characters found in sequence

    removeLastChar  Remove last character in each line

    spotNumAtNameStart Put spot number at start of spot name followed by an underscore
                       Works for column-based submissions only to retain non-unique
                       column values as part of the spot name (e.g. seq occurrence counts)

    allowEarlyFileEnd For 3 or more files, complete load at early end of one of the files
                      (because that is all they have)

    seqOnlyFile:    File contains sequence only (filename minus suffixes is used for spot name)

    genbankFile:    File is in a GenBank format

    noPairFileCheck Prevents check for file pairings (useful for large number of files that
                    don't need this check.)

    convertDeflineSpaces Converts each whitespace group in a defline to an underscore in
                    the defline (in order to incorporate additional information into the
                    spot name)

    addFileToName   Prepends the fastq file name (minus suffixes) to the spot name

    ignoreSGforPairing   Ignores spot group when checking for file pairing

    useSharq        Run the binary sharq if the files are compatible

    doNotUseSharq   Do not run the binary sharq even if the files seem compatible

    schema:         Set vdb schema to use during load

    z|xml-log:      XML version of stderr output

    h|help:         Displays this message

    V|version:      Displays version of fastq-load.py

You may have to provide read1 and read2 lists if pairing of files is necessary and
orphan reads exist. If both reads exist in the same file, put the same filename
in both lists.

Log odds will be recognized if the quality offset is determined and not provided.
If offset is provided and data is log odds, then set 'logodds' option prior to
setting 'offset' or 'quality' options.

"""
# pylint: disable=line-too-long,invalid-name,missing-function-docstring
# pylint: disable=too-many-lines,consider-using-with,consider-using-in
# pylint: disable=too-many-boolean-expressions,too-many-statements,unused-variable
# pylint: disable=too-many-nested-blocks,too-many-return-statements,anomalous-backslash-in-string
# pylint: disable=too-many-return-statements,too-many-instance-attributes
# pylint: disable=bare-except,too-many-branches,too-many-public-methods
# pylint: disable=too-many-locals,consider-using-dict-items,global-statement
# pylint: disable=too-many-arguments,consider-using-generator,superfluous-parens

import sys
import array
from enum import Enum
import logging
import os
import re
import copy
import gzip
import bz2
import datetime
import time
import subprocess
import getpass
import cProfile
import pstats
import sh

try:
    from backports import lzma
    from pyth.plugins.rtf15.reader import Rtf15Reader
    from pyth.plugins.plaintext.writer import PlaintextWriter
    from collections import OrderedDict
    from operator import itemgetter
    from fastqLoadCLib import cParseQual
    import GeneralWriter
except ModuleNotFoundError:
    if __name__ == '__main__':
        raise

############################################################
# Environment globals
############################################################

filePaths = {}    # map from filename to file path
fileHandles = {}  # map from filename to file handle
fileObjects = {}  # map from filename to file object for determining percent read
fileSizes = {}    # map from filename to file size for determining percent read
fileSkip = {}     # hash of filenames that have been paired or found to be quality files
fileTypes = {}    # map from filename to file type
                  # normal, singleLine, eightLine, multiLine, multiLineEightLine,
                  # seqQual, multiLineSeqQual, fasta, multiLineFasta
fileReadPairs = {}   # map from filename to file that contains read pairs
fileReadPairs2 = {}  # map from filename to file that contains read pairs (if three reads)
fileReadPairs3 = {}  # map from filename to file that contains read pairs (if four reads)
fileReadPairs4 = {}  # map from filename to file that contains read pairs (if five reads)
fileReadPairs5 = {}  # map from filename to file that contains read pairs (if six reads)
filePore2D = {}      # map from nanopore template file to nanopore 2D file
fileQualPairs = {}   # map from filename to file containing qualities
fileLabels = {}      # map from filename to label for file
fileConcatStart = {} # retain location where second read information starts in file
fileWithDupReads = {}# retain files that contain duplicate reads.
fileSpotGroups = {}  # map from filenames to barcode (1st file in pairs only, no qual files)
filesToTruncate = {} # hash of files to truncate

sw = None
statusWriter = None
version = "1.1.2"

############################################################
# Output usage statement along with message if specified
############################################################

def usage ( message, state ):
    usageMessage =  """
Usage: fastq-load.py\n
          [ --output=<archive>        (optional) ]
          [ --offset=<number>         (optional) ]
          [ --readLens=<number(s)>    (optional) ]
          [ --readTypes=<B|T(s)>      (optional) ]
          [ --spotGroup=<string>      (optional) ]
          [ --orphanReads             (optional) ]
          [ --logOdds                 (optional) ]
          [ --ignoreNames             (optional) ]
          [ --ignAndDiscardNames      (optional) ]
          [ --useAndDiscardNames      (optional) ]
          [ --fabricateNames          (optional) ]
          [ --read1PairFiles=<string> (optional) ]
          [ --read2PairFiles=<string> (optional) ]
          [ --read3PairFiles=<string> (optional) ]
          [ --read4PairFiles=<string> (optional) ]
          [ --read5PairFiles=<string> (optional) ]
          [ --read6PairFiles=<string> (optional) ]
          [ --read1QualFiles=<string> (optional) ]
          [ --read2QualFiles=<string> (optional) ]
          [ --filesToignore=<string>  (optional) ]
          [ --filesToTruncate=<string>(optional) ]
          [ --fileSpotGroups=<string> (optional) ]
          [ --useFilenameForSG        (optional) ]
          [ --platform=<string>       (optional) ]
          [ --readLabels=<labels>     (optional) ]
          [ --maxErrorCount=<count>   (optional) ]
          [ --maxErrorOutput=<count>  (optional) ]
          [ --maxSearchCount=<count>  (optional) ]
          [ --maxDeflineLen=<count>   (optional) ]
          [ --maxSeqLineCount=<count> (optional) ]
          [ --mixedDeflines           (optional) ]
          [ --genericDeflines         (optional) ]
          [ --mixedTypes              (optional) ]
          [ --isMultiLine             (optional) ]
          [ --ignLeadChars=<count>    (optional) ]
          [ --ignTailChars=<count>    (optional) ]
          [ --ignCharsAfterPlus       (optional) ]
          [ --discardBarcodes         (optional) ]
          [ --discardBCchars          (optional) ]
          [ --appendBCtoName          (optional) ]
          [ --bcSpaceOffset           (optional) ]
          [ --concatPairFiles         (optional) ]
          [ --removeBadChars          (optional) ]
          [ --removeSeqSpaces         (optional) ]
          [ --badFlatNumQual          (optional) ]
          [ --convertEmptyDeflines    (optional) ]
          [ --convertDeflineSpaces    (optional) ]
          [ --ignoreSGforPairing      (optional) ]
          [ --spotNumAtNameStart      (optional) ]
          [ --allowEarlyFileEnd       (optional) ]
          [ --seqOnlyFile             (optional) ]
          [ --genbankFile             (optional) ]
          [ --noPairFileCheck         (optional) ]
          [ --addFileToName           (optional) ]
          [ --fastqWithGtrThan        (optional) ]
          [ --fastqWithoutAtSign      (optional) ]
          [ --removeLastChar          (optional) ]
          [ --useSharq                (optional) ]
          [ --doNotUseSharq           (optional) ]
          [ --abiLastPrimerBase=<T|G> (optional) ]
          [ --schema=<string>         (optional) ]
          [ --nameColumns=<string>    (optional) ]
          [ --seqColumns=<string>     (optional) ]
          [ --qualColumns=<string>    (optional) ]
          [ --groupColumns=<string>   (optional) ]
          [ -z|--xml-log=<string>     (optional) ]
          [ -h|--help                 (optional) ]
          [ -V|--version              (optional) ]
          [ fastq-file(s)                        ]

"""
    if message:
        sys.stderr.write( "\n{}\n".format(message) )
    sys.stderr.write(usageMessage)
    sys.stderr.write(__doc__)
    sys.exit(state)

############################################################
# Defline StatusWriter class
############################################################

class StatusWriter:
    """ Outputs status to stderr and optionally to an xml log file """
    def __init__(self, vers):
        self.vers = vers
        self.xmlLogHandle = None
        self.pid = os.getpid()

    ############################################################
    # Open xml log file if provided
    ############################################################

    def setXmlLog ( self, xmlLogFile ):
        xmlLogFile = xmlLogFile.strip()
        try:
            self.xmlLogHandle = open(xmlLogFile, 'w')
            self.xmlLogHandle.write("<Log>\n")
        except OSError:
            sys.exit( "\nFailed to open {} for writing\n\n".format(xmlLogFile) )

    ############################################################
    # Close Xml Log file
    ############################################################
    def closeXmlLog ( self ):
        self.xmlLogHandle.write("</Log>\n")
        self.xmlLogHandle.close()

    ############################################################
    # Output status message
    ############################################################
    def outputInfo ( self, message, percent=None ):
        dateTime = self.getTime()
        if self.xmlLogHandle:
            if percent:
                self.xmlLogHandle.write('<info app="fastq-load.py" message="{}" pid="{}" timestamp="{}" version="{}" percent="{}"/>\n'
                                        .format(self.escape(message),self.pid,dateTime,self.vers,percent))
            else:
                self.xmlLogHandle.write('<info app="fastq-load.py" message="{}" pid="{}" timestamp="{}" version="{}"/>\n'
                                        .format(self.escape(message),self.pid,dateTime,self.vers))
            self.xmlLogHandle.flush()
        sys.stderr.write("{} fastq-load.py.{} info: {}\n".format(dateTime,self.vers,message) )
        sys.stderr.flush()

    ############################################################
    # Output status message
    ############################################################
    def outputWarning ( self, message ):
        dateTime = self.getTime()
        if self.xmlLogHandle:
            self.xmlLogHandle.write('<warning app="fastq-load.py" message="{}" pid="{}" timestamp="{}" version="{}"/>\n'
                                    .format(self.escape(message),self.pid,dateTime,self.vers))
            self.xmlLogHandle.flush()
        sys.stderr.write("{} fastq-load.py.{} warn: {}\n".format(dateTime,self.vers,message) )
        sys.stderr.flush()

    ############################################################
    # Output status message and exit
    ############################################################
    def outputErrorAndExit (self, message):
        dateTime = self.getTime()
        if self.xmlLogHandle:
            self.xmlLogHandle.write('<error app="fastq-load.py" message="{}" pid="{}" timestamp="{}" version="{}"/>\n'
                                    .format(self.escape(message),self.pid,dateTime,self.vers))
            self.xmlLogHandle.flush()
            self.closeXmlLog()
        sys.exit( "\n{} fastq-load.py.{} Error: {}".format(dateTime,self.vers,message) )

    ############################################################
    # Escape message (wrote my own instead of importing sax escape)
    ############################################################
    @staticmethod
    def escape(message):
        message = message.replace('&', '&amp;')
        message = message.replace('"', '&quot;')
        message = message.replace("'", '&apos;')
        message = message.replace('<', '&lt;')
        message = message.replace('>', '&gt;')
        return message

    ############################################################
    # Get formatted date time
    ############################################################
    @staticmethod
    def getTime():
        now = datetime.datetime.utcnow()
        now = now.replace(microsecond=0)
        return now.isoformat()

############################################################
# Define Platform and convert platform to SRA enum
############################################################

class Platform(Enum):
    """ Enumeration class for SRA platforms """
    SRA_PLATFORM_UNDEFINED         = 0
    SRA_PLATFORM_454               = 1
    SRA_PLATFORM_ILLUMINA          = 2
    SRA_PLATFORM_ABSOLID           = 3
    SRA_PLATFORM_COMPLETE_GENOMICS = 4
    SRA_PLATFORM_HELICOS           = 5
    SRA_PLATFORM_PACBIO_SMRT       = 6
    SRA_PLATFORM_ION_TORRENT       = 7
    SRA_PLATFORM_CAPILLARY         = 8
    SRA_PLATFORM_OXFORD_NANOPORE   = 9

    @classmethod
    def convertPlatformString ( cls, platformString ):
        platformString = platformString.upper()
        if platformString in ("454", "LS454"):
            return cls.SRA_PLATFORM_454
        if platformString == "ILLUMINA":
            return cls.SRA_PLATFORM_ILLUMINA
        if platformString in ("ABI", "SOLID", "ABSOLID", "ABISOLID"):
            return cls.SRA_PLATFORM_ABSOLID
        if platformString in ("PACBIO", "PACBIO_SMRT"):
            return cls.SRA_PLATFORM_PACBIO_SMRT
        if platformString in ("CAPILLARY", "SANGER"):
            return cls.SRA_PLATFORM_CAPILLARY
        if platformString == "NANOPORE":
            return cls.SRA_PLATFORM_OXFORD_NANOPORE
        if platformString == "HELICOS":
            return cls.SRA_PLATFORM_HELICOS
        if platformString == "ION_TORRENT":
            return cls.SRA_PLATFORM_ION_TORRENT
        if platformString in ("BGI", "UNDEFINED", "MIXED"):
            return cls.SRA_PLATFORM_UNDEFINED
        return None

############################################################
# Defline class (not an enum)
############################################################

class Defline:
    """ Retains information parsed from fastq defline """

    ILLUMINA_NEW                = 1
    ILLUMINA_OLD                = 2
    PACBIO                      = 3
    ABSOLID                     = 4
    QIIME_ILLUMINA_NEW          = 5
    QIIME_ILLUMINA_NEW_BC       = 6
    ILLUMINA_NEW_DOUBLE         = 7
    LS454                       = 8
    UNDEFINED                   = 9
    ILLUMINA_OLD_BC_RN          = 10
    ILLUMINA_NEW_WITH_SUFFIX    = 11
    QIIME_ILLUMINA_OLD          = 12
    QIIME_ILLUMINA_OLD_BC       = 13
    QIIME_454                   = 14
    QIIME_454_BC                = 15
    QIIME_ILLUMINA_NEW_DBL      = 16
    QIIME_ILLUMINA_NEW_DBL_BC   = 17
    QIIME_GENERIC               = 18
    READID_BARCODE              = 19
    ILLUMINA_NEW_OLD            = 20
    NANOPORE                    = 21
    HELICOS                     = 22
    ION_TORRENT                 = 23
    SANGER_NEWBLER              = 24
    ALTNUM                      = 25
    ABSOLID2                    = 26
    ILLUMINA_NEW_DATA_GRP       = 27
    BGI                         = 28

    def __init__(self, deflineString):
        self.deflineType = None
        self.deflineStringOrig = ''
        self.deflineString = ''
        self.saveDeflineType = True
        self.name = ''
        self.platform = "NotSet"
        self.isValid = False
        self.readNum = ""
        self.filterRead = 0
        self.spotGroup = ''
        self.prefix = ''
        self.lane = ''
        self.tile = ''
        self.x = ''
        self.y = ''
        self.numDiscards = 0
        self.foundRE = None
        self.dateAndHash454 = ''
        self.region454 = ''
        self.xy454 = ''
        self.qiimeName = ''
        self.filename = ''
        self.filenameTrunc = ''
        self.poreRead = ''
        self.poreFile = ''
        self.channel = ''
        self.readNo = 0 # Not to be confused with readNum which is like 1/2; readNo is a nanopore characteristic
        self.flowcell = ''
        self.field = ''
        self.camera = ''
        self.position = ''
        self.runId = ''
        self.row = ''
        self.column = ''
        self.dir = ''
        self.panel = ''
        self.tagType = ''
        self.suffix = ''
        self.abiTitle = ''
        self.statusWriter = None
        self.ignLeadCharsNum = None
        self.ignoredLeadChars = None
        self.ignoredTailChars = None
        self.ignTailCharsNum = 0
        self.ignCharsAfterPlus = False
        self.altNum = None # Alternately arrived at spot or read numbers
        self.retainAltNum = False # If large alternate numbers (>4) encountered, then altNum is retained as part of the read name, otherwise it is treated as a read number
        self.bcSpaceOffset = False
        self.forceGeneric = False
        self.addFileToName = False
        self.convertDeflineSpaces = False
        self.appendPoreReadToName = False
        self.fastqWithGtrThan = False
        self.fastqWithoutAtSign = False

        self.bgiSelected = None
        self.bgiOldFormat = re.compile(r"[@>+](\S{1,3}\d{9}\S{0,3})(L\d{1})(C\d{3})(R\d{3})([_]?\d{1,8})(#[!-~]*?|)(/[1234]\S*|)(\s+|$)")
        self.bgiNewFormat = re.compile(r"[@>+](\S{1,3}\d{9}\S{0,3})(L\d{1})(C\d{3})(R\d{3})([_]?\d{1,8})(\S*)(\s+|[_|-])([12345]|):([NY]):(\d+|O):?([!-~]*?)(\s+|$)")

        self.illuminaNewSelected = None
        self.illuminaNew = re.compile(r"[@>+]([!-~]+?)([:_])(\d+)([:_])(\d+)([:_])(-?\d+\.?\d*)([:_])(-?\d+\.\d+|\d+)(\s+|[_|-])([12345]|):([NY]):(\d+|O):?([!-~]*?)(\s+|$)")
        self.illuminaNewNoPrefix = re.compile(r"[@>+]([!-~]*?)(:?)(\d+)([:_])(\d+)([:_])(\d+)([:_])(\d+)(\s+|_)([12345]|):([NY]):(\d+|O):?([!-~]*?)(\s+|$)")
        self.illuminaNewWithSuffix = re.compile(r"[@>+]([!-~]+)([:_])(\d+)([:_])(\d+)([:_])(-?\d+\.?\d*)([:_])(-?\d+\.\d+|\d+)([!-~]+?\s+|[!-~]+?[:_|-])([12345]|):([NY]):(\d+|O):?([!-~]*?)(\s+|$)")
        self.illuminaNewWithUnderscores = re.compile(r"[@>+]([!-~]+?)(_)(\d+)(_)(\d+)(_)(\d+)(_)(\d+)(\s+|_)([12345]|)_([NY])_(\d+|O)_?([!-~]*?)(\s+|$)")
        self.illuminaNewSuffix = re.compile(r"(#[!-~]*?|)(/[12345]|\\[12345])?([!-~]*?)(#[!-~]*?|)(/[12345]|\\[12345])?([:_|]?)(\s+|$)")

        self.illuminaOldSelected = None
        self.illuminaOldColon = re.compile(r"[@>+]?([!-~]+?)(:)(\d+)(:)(\d+)(:)(-?\d+\.?\d*)([-:])(-?\d+\.\d+|-?\d+)_?[012]?(#[!-~]*?|)\s?(/[12345]|\\[12345])?(\s+|$)")
        self.illuminaOldUnderscore = re.compile(r"[@>+]?([!-~]+?)(_)(\d+)(_)(\d+)(_)(-?\d+\.?\d*)(_)(-?\d+\.\d+|-?\d+)(#[!-~]*?|)\s?(/[12345]|\\[12345])?(\s+|$)")
        self.illuminaOldNoPrefix = re.compile(r"[@>+]?([!-~]*?)(:?)(\d+)(:)(\d+)(:)(-?\d+\.?\d*)(:)(-?\d+\.\d+|-?\d+)(#[!-~]*?|)\s?(/[12345]|\\[12345])?(\s+|$)")
        self.illuminaOldWithSuffix = re.compile(r"[@>+]?([!-~]+?)(:)(\d+)(:)(\d+)(:)(-?\d+\.?\d*)(:)(-?\d+\.\d+|-?\d+)(#[!-~]+)(/[12345])[!-~]+(\s+|$)")
        self.illuminaOldWithSuffix2 = re.compile(r"[@>+]?([!-~]+?)(:)(\d+)(:)(\d+)(:)(-?\d+\.?\d*)(:)(-?\d+\.?\d*[!-~]+?)(#[!-~]*?|)\s?(/[12345]|\\[12345])?(\s+|$)")
        self.illuminaOldSuffix = re.compile(r"(-?\d+\.\d+|-?\d+)([!-~]*)") # Must have '*' and not '+'. Otherwise, name for pairing is truncated by one character.
#        self.illuminaOldSuffix2 = re.compile(r"(-?\d+\.\d+|-?\d+):([NY]):(\d+|O):?([!-~]*)")

        self.illuminaOldBcRnSelected = None
        self.illuminaOldBcRnOnly = re.compile(r"[@>+]([!-~]+?)(#[!-~]+?)(/[12345]|\\[12345])(\s+|$)")
        self.illuminaOldBcOnly = re.compile(r"[@>+]([!-~]+?)(#[!-~]+)(\s+|$)(.?)")
        self.illuminaOldRnOnly = re.compile(r"[@>+]([!-~]+?)(/[12345]|\\[12345])(\s+|$)(.?)")

        self.illuminaNewDataGroup = re.compile(r"[@>+]([!-~]+?)(\s+|[_|])([12345]|):([NY]):(\d+|O):?([!-~]*?)(\s+|$)")

        self.altNumSelected = None
        self.altNumOnly = re.compile(r"[@>+]([!-~]+?)([_\-/|])(\d+)(\s+|$).?")
        self.altNumOnly2 = re.compile(r"[@>+]([!-~]+?)([_\-/|])([FR])(\s+|$).?")

        self.qiimeSelected = None
        self.qiimeBc = re.compile(r"[@>+]([!-~]*).*?\s+orig_bc=[!-~]+\s+new_bc=([!-~]+)\s+bc_diffs=[01]")
        self.qiimeBc2 = re.compile(r"[@>+]([!-~]*).*?\s+([!-~]+)\s+orig_bc=[!-~]+\s+new_bc=([!-~]+)\s+bc_diffs=[01]")
        self.qiimeIlluminaNew = re.compile(r"[@>+]([!-~]*)\s+([!-~]*?)([:_])(\d+)([:_])(\d+)([:_])(\d+)([:_])(\d+)(\s+|_|:)([12345]):([NY]):(\d+|O):?([!-~]*?)(\s+|$)")
        self.qiimeIlluminaNew2 = re.compile(r"[@>+]([!-~]*)\s+([!-~]+)\s+([!-~]*?)([:_])(\d+)([:_])(\d+)([:_])(\d+)([:_])(\d+)(\s+|_|:)([12345]):([NY]):(\d+|O):?([!-~]*?)(\s+|$)")
        self.qiimeIlluminaNewPeriods = re.compile(r"[@>+]([!-~]*)\s+([!-~]*?)(\.)(\d+)(\.)(\d+)(\.)(\d+)(\.)(\d+)(\s+|_|:)([12345])\.([NY])\.(\d+|O)\.?([!-~]*?)(\s+|$)")
        self.qiimeIlluminaNewUnderscores = re.compile(r"[@>+]([!-~]*)\s+([!-~]*?)(_)(\d+)(_)(\d+)(_)(\d+)(_)(\d+)(\s+|_|:)([12345])_([NY])_(\d+|O)_?([!-~]*?)(\s+|$)")

#        self.qiimeIlluminaOld = re.compile(r"[@>+]([!-~]*\s?[!-~]*)\s+([!-~]+?)(:)(\d+)(:)(\d+)(:)(-?\d+)(:)(-?\d+)(#[!-~]*?|)\s?(/[12345]|\\[12345])?(\s+|$)")
        self.qiimeIlluminaOld = re.compile(r"[@>+]([!-~]*)\s+([!-~]+?)(:)(\d+)(:)(\d+)(:)(-?\d+)(:)(-?\d+)(#[!-~]*?|)\s?(/[12345]|\\[12345])?(\s+|$)")

        self.ls454 = re.compile(r"[@>+]([!-~]+_|)([A-Z0-9]{7})(\d{2})([A-Z0-9]{5})(/[12345])?(\s+|$)")
        self.qiime454 = re.compile(r"[@>+]([!-~]*)\s+([A-Z0-9]{7})(\d{2})([A-Z0-9]{5})(/[12345])?(\s+|$)")

        self.pacbioSelected = None
        self.pacbio = re.compile(r"[@>+](m\d{5,6}_\d{6}_[!-~]+?_c\d{33}_s\d+_[pX]\d/\d+/?\d*_?\d*|m\d{6}_\d{6}_[!-~]+?_c\d{33}_s\d+_[pX]\d[|/]\d+[|/]ccs[!-~]*?)(\s+|$)")
        self.pacbio2 = re.compile(r"[@>+]([!-~]*?m\d{5,6}\S{0,3}_\d{6}_\d{6}[/_]\d+[!-~]*?)(\s+|$)")
        self.pacbio3 = re.compile(r"[@>+]([!-~]*?m\d{5,6}\S{0,3}_\d{6}_\d{6}[/_]\d+/ccs[!-~]*?)(\s+|$)")
        self.pacbio4 = re.compile(r"[@>+]([!-~]*?m\d{5,6}\S{0,3}_\d{6}_\d{6}[/_]\d+/\d+_\d+[!-~]*?)(\s+|$)")

        self.nanoporeSelected = None
        self.nanopore = re.compile(r"[@>+]+?(channel_)(\d+)(_read_)?(\d+)?([!-~]*?)(_twodirections|_2d|-2D|_template|-1D|_complement|-complement|\.1C|\.1T|\.2D)?(:[!-~ ]+?_ch\d+_file\d+_strand.fast5)?(\s+|$)")
        self.nanopore2 = re.compile(r"[@>+]([!-~]*?ch)(\d+)(_file)(\d+)([!-~]*?)(_twodirections|_2d|-2D|_template|-1D|_complement|-complement|\.1C|\.1T|\.2D)(:[!-~ ]+?_ch\d+_file\d+_strand.fast5)?(\s+|$)")
        self.nanopore3 =   re.compile(r"[@>+]([!-~]*?)[: ]?([!-~]+?Basecall)(_[12]D[_0]*?|_Alignment[_0]*?|_Barcoding[_0]*?|)(_twodirections|_2d|-2D|_template|-1D|_complement|-complement|\.1C|\.1T|\.2D|)[: ]([!-~]*?)[: ]?([!-~ ]+?_ch)_?(\d+)(_read|_file)_?(\d+)(_strand\d*.fast5|_strand\d*.*|)(\s+|$)")
        self.nanopore3_1 = re.compile(r"[@>+]([!-~]+?)[: ]?([!-~]+?Basecall)(_[12]D[_0]*?|_Alignment[_0]*?|_Barcoding[_0]*?|)(_twodirections|_2d|-2D|_template|-1D|_complement|-complement|\.1C|\.1T|\.2D|)[: ]([!-~]*?)[: ]?([!-~ ]+?_read_)(\d+)(_ch_)(\d+)(_strand\d*.fast5|_strand\d*.*)(\s+|$)")
        self.nanopore4 = re.compile(r"[@>+]([!-~]*?\S{8}-\S{4}-\S{4}-\S{4}-\S{12}\S*[_]?\d?)[\s+[!-~ ]*?|]$")
        self.nanopore5 = re.compile(r"[@>+]([!-~]*?[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{12}_Basecall)(_[12]D[_0]*?|_Alignment[_0]*?|_Barcoding[_0]*?)(_twodirections|_2d|-2D|_template|-1D|_complement|-complement|\.1C|\.1T|\.2D)\S*?($)")

        self.helicos = re.compile(r"[@>+](VHE-\d+)-(\d+)-(\d+)-(\d)-(\d+)(\s+|$)")

        self.ionTorrentSelected = None
        self.ionTorrent = re.compile(r"[@>+]([A-Z0-9]{5})(:)(\d{1,5})(:)(\d{1,5})(/[12345]|\\[12345]|[LR])?(\s+|$)")
        self.ionTorrent2 = re.compile(r"[@>+]([A-Z0-9]{5})(:)(\d{1,5})(:)(\d{1,5})(\s+|[_|])([12345]|):([NY]):(\d+|O):?([!-~]*?)(\s+|$)")

        self.abSolid = re.compile(r"[@>+]([!-~]*?)(\d+)_(\d+)_(\d+)_(F3|R3|F5-BC|BC|F5-P2|F5-RNA|F5-DNA|FC[1-9])([!-~]*?)(\s+|$)")
        self.abSolid2 = re.compile(r"[@>+](\d{1,4})_(\d{1,4})_(\d{1,4})(\s+|$)")

        self.sangerNewbler = re.compile(r"[@>+]([!-~]+?)\s+template=([!-~]+)\s+dir=([!-~]+)(\s+|$)")

        self.readIdBarcode = re.compile(r"[@>+]([!-~]+)\s+read_id=([!-~]*?::|)([!-~]+)\s+barcode=([!-~]+).*(\s+|$)")

        self.generic = re.compile(r"[@>+]([!-~]+)(\s+|$)")

        self.sharqCompatible = False

        if deflineString:
            self.parseDeflineString(deflineString)

    ############################################################
    # Determine if deflines from two files represent paired reads
    # or paired seq/qual files. If testing for paired seq/qual,
    # then checkReadNum is True. Note that name could be '0'.
    ############################################################

    @classmethod
    def isPairedDeflines( cls, defline1, defline2, sameReadNum, ignoreSGforPairing=False, discardBCchars=0, useFileOrder=0 ):

        # This case is used to determine if a seq files goes with a
        # qual file. In this case the read numbers must be the same

        spotGroup1 = defline1.spotGroup
        spotGroup2 = defline2.spotGroup
        if discardBCchars != 0:
            if discardBCchars < 0:
                spotGroup1 = spotGroup1[0:len(spotGroup1) + discardBCchars]
                spotGroup2 = spotGroup2[0:len(spotGroup2) + discardBCchars]
            else:
                spotGroup1 = spotGroup1[discardBCchars:]
                spotGroup2 = spotGroup2[discardBCchars:]

        if ( sameReadNum and
             defline1.name and
             defline2.name and
             ( defline1.name == defline2.name or
               ( defline1.abiTitle and defline1.name == defline1.abiTitle + "_" + defline2.name ) or
               ( defline2.abiTitle and defline2.name == defline2.abiTitle + "_" + defline1.name ) ) and
             defline1.readNum == defline2.readNum ):

            if ( defline1.tagType and
                 defline2.tagType and
                 not defline1.tagType == defline2.tagType ):
                return 0
            else:
                return 1

        # This case is used to determine if two files should be
        # paired together (i.e. different read numbers)

        elif ( not sameReadNum and
               defline1.name and
               defline2.name and
               defline1.name == defline2.name and
               ( ignoreSGforPairing or spotGroup1 == spotGroup2 ) ):

            # Determine which defline is associated with first read

            if useFileOrder:
                if defline1.filename < defline2.filename:
                    return 1
                elif defline1.filename > defline2.filename:
                    return 2

            elif ( defline1.readNum and
                   defline2.readNum ):

                # Return 1 if first read associated with defline1

                if defline1.readNum <= defline2.readNum:
                    return 1

                # Return 2 if first read associated with defline2

                else:
                    return 2

            # Check for nanopore template (1D) reads

            elif ( defline1.poreRead and
                   defline1.poreRead == "template" ):
                return 1

            elif ( defline2.poreRead and
                   defline2.poreRead == "template" ):
                return 2

            # Not pairing complement and 2D nanopore reads

            elif ( defline1.poreRead or
                   defline2.poreRead ):
                return 0

            # Check for absolid reads

            elif ( defline1.tagType and
                   defline2.tagType ):

                if ( defline1.tagType == "F3" and
                     defline2.tagType != "F3" ):
                    return 1

                elif ( defline2.tagType == "F3" and
                       defline1.tagType != "F3" ):
                    return 2

                else:
                    return 1

            # Compare defline strings to determine order if not read numbers

            elif defline1.deflineString < defline2.deflineString:
                return 1

            elif defline1.deflineString > defline2.deflineString:
                return 2

            # Compare filenames

            elif defline1.filename < defline2.filename:
                return 1

            elif defline1.filename > defline2.filename:
                return 2

            # Defaulting to first defline associated with first read
            # if indeterminate via readNums and defline order

            else:
                return 1
        else:
            return 0

    ############################################################
    # Reset object variables
    ############################################################
    def reset(self):
        self.deflineStringOrig = ''
        self.deflineString = ''
        self.name = ''
        self.isValid = True
        self.readNum = ''
        self.filterRead = 0
        self.spotGroup = ''
        self.prefix = ''
        self.lane = ''
        self.tile = ''
        self.x = ''
        self.y = ''
        self.dateAndHash454 = ''
        self.region454 = ''
        self.xy454 = ''
        self.qiimeName = ''
        self.poreRead = ''
        self.poreFile = ''
        self.channel = 0
        self.readNo = 0
        self.runId = ''
        self.row = ''
        self.column = ''
        self.camera = ''
        self.field = ''
        self.position = ''
        self.flowcell = ''
        self.panel = ''
        self.tagType = ''
        self.suffix = ''
        self.dir = ''
        self.ignoredLeadChars = None
        self.foundRE = None

    ############################################################
    # Set status from spot writer
    ############################################################

    def setStatus(self,status):
        self.statusWriter = status.statusWriter
        self.retainAltNum = status.retainAltNum
        self.bcSpaceOffset = status.bcSpaceOffset
        self.addFileToName = status.addFileToName
        self.filenameTrunc = status.filenameTrunc
        self.fastqWithGtrThan = status.fastqWithGtrThan
        self.fastqWithoutAtSign = status.fastqWithoutAtSign
        self.convertDeflineSpaces = status.convertDeflineSpaces
        if status.ignoreNames:
            self.deflineType = self.UNDEFINED
        if status.mixedDeflines:
            self.saveDeflineType = False
        if status.genericDeflines:
            self.forceGeneric = True
        self.ignLeadCharsNum = status.ignLeadCharsNum
        self.ignTailCharsNum = status.ignTailCharsNum
        self.ignCharsAfterPlus = status.ignCharsAfterPlus
        self.appendPoreReadToName = status.appendPoreReadToName

    ############################################################
    # Parse deflines based on cases we have seen in the past
    ############################################################

    readNumberAtEnd = re.compile(r"/[1234]$")
    readNumberInSuffix = re.compile(r"^/[12]")
    floatInDiscards = re.compile(r"\.")
    readNumInSpotGroup = re.compile(r"/[12]")
    getPoreReadNo = re.compile(r"read[=_]?(\d+)")
    getPoreChannel = re.compile(r"ch[=_]?(\d+)")
    getPoreBarcode = re.compile(r"barcode=(\S+)")
    getPoreReadNo2 = re.compile(r"read_?(\d+)")
    getPorePass = re.compile(r'pass[/\\]')
    getPoreFail = re.compile(r"fail[/\\]")
    getPoreBarcode2 = re.compile(r'(NB\d{2}|BC\d{2}|barcode\d{2})([/\\])')
    barcodeStringPresent = re.compile(r"^barcode")
    pore2Dpresent = re.compile(re.escape(".2D."))
    poreTemplatePresent = re.compile(re.escape(".template."))
    poreComplementPresent = re.compile(re.escape(".complement."))
    alphanumericPresent = re.compile(r"[a-zA-Z0-9]")

    def parseDeflineString(self,deflineString):

        self.reset()
        self.deflineStringOrig = deflineString
        self.deflineString = self.deflineStringOrig.strip()

        # Find/remove occurrence of space after :, #, etc. for illumina.

        if self.deflineString:
            if self.bcSpaceOffset:
                self.deflineString = re.sub(r'([:#]) (\S+)$', r'\1\2', self.deflineString, 1)
            if self.convertDeflineSpaces:
                self.deflineString = re.sub(r'\s+','_',self.deflineString)
            if self.fastqWithGtrThan and self.deflineString[0] == '>':
                self.deflineString = '@' + self.deflineString[1:]
            if self.fastqWithoutAtSign and not self.deflineString[0] in "@+" :
                self.deflineString = '@' + self.deflineString

        # End of file produces empty defline string

        if not self.deflineString:
            self.isValid = False

        elif ( self.deflineString == '+' or
               ( self.ignCharsAfterPlus and
                 self.deflineString[0] == '+') ):
            self.isValid = True

        ############################################################
        # Force generic defline
        #
        # @m150529_110144_42160_c100830212550000001823183611251532_s1_p0|40117|ccs|CDK4|chr12:58141509-58146229:-1 (SRR6360536)
        # @m150529_110144_42160_c100830212550000001823183611251532_s1_p0|506|ccs|TUBB|chr6:30688156-30693194:1 (SRR6360536)
        ############################################################

        elif self.forceGeneric:

            # Capture name after '@' sign

            m = self.generic.match( self.deflineString )

            # Confirm regular expression succeeded

            if m is None :
                self.isValid = False
                return self.isValid

            # Get match value

            self.name = m.group(1)
            self.platform = "UNDEFINED"

        ############################################################
        # helicos
        #
        # @VHE-242383071011-15-1-0-2 (YYY034449)
        ############################################################

        elif ( self.deflineType == self.HELICOS or
               ( self.deflineType is None and
                 self.helicos.match(self.deflineString) ) ):

            m = self.helicos.match(self.deflineString)

            # Confirm regular expression succeeded

            if m is None :
                self.isValid = False
                return self.isValid

            # Get match values

            (self.flowcell, self.channel, self.field, self.camera, self.position, endSep) = m.groups()

            # Set defline name

            self.name = f'{self.flowcell}-{self.channel}-{self.field}-{self.camera}-{self.position}'

            # Retain defline type if desired

            if ( not self.deflineType and
                 self.saveDeflineType ):
                self.deflineType = self.HELICOS
                self.platform = "HELICOS"

        ############################################################
        # ab solid
        #
        # >3_189_730_F3 (XXX001662)
        # >3_189_730_R3 (XXX001662)
        # >461_28_1048_F3 (XXX001354)
        # >349_1793_467_F3 (XXX015374)
        # >1_98_123_BC (XXX1001434)
        # >427_21_101_F3 (XXX112693)
        # >427_17_22_F5-P2 (XXX112693)
        # >2_35_407_F3 (XXX3159530)
        # >2_35_407_F5-BC (XXX3159530)
        # @1_43_495_F3 (XXX645822)
        # @1_43_495_F5-BC (XXX645822)
        # >47_15_207_F5-RNA (YYY005001)
        # >47_15_207_F5-DNA (YYY005002)
        # >1_01_1_23_192_F3 (YYY3222665)
        # >1_23_192_F3_1_01 (ZZZ3222665)
        # >427_27_224_F3 1:0003213231 (ZZZ005000)
        # @1_11_357_FC1_01_Library21 (ERR2135869)
        # @1_3_235_FC3_01_Library38 (ERR2135856)
        # @2_78_182_BC4 (ERR2135858)
        # @1_106_125 (ERR1924251) - backed out implementation
        ############################################################

        elif ( self.deflineType == self.ABSOLID or
               ( self.deflineType is None and
                 self.abSolid.match(self.deflineString) ) ):

            m = self.abSolid.match(self.deflineString)

            # Confirm regular expression succeeded

            if m is None :
                self.isValid = False
                return self.isValid

            # Get match values (prefix and/or suffix may contribute to making name unique)

            (self.prefix, self.panel, self.x, self.y, self.tagType, self.suffix, endSep) = m.groups()

            # Set defline name

            if ( self.tagType[0:2] == "FC" or
                 not self.suffix == '' ) :
                self.suffix = f'_{self.tagType}{self.suffix}'
                self.tagType = ''

            if self.abiTitle:
                self.name = f'{self.abiTitle}_{self.prefix}{self.panel}_{self.x}_{self.y}'
            else:
                self.name = f'{self.prefix}{self.panel}_{self.x}_{self.y}'

            # Retain defline type if desired

            if ( not self.deflineType and
                 self.saveDeflineType ):
                self.deflineType = self.ABSOLID
                self.platform = "ABSOLID"

        ############################################################
        # BGI
        #
        ############################################################

        elif ( self.deflineType == self.BGI or
               ( self.deflineType is None and
                 ( self.bgiOldFormat.match ( self.deflineString ) or
                   self.bgiNewFormat.match ( self.deflineString ) ) ) ):

            # Capture information after '@' sign

            m = None
            if ( self.deflineType and
                 self.bgiSelected ) :
                m = self.bgiSelected.match ( self.deflineString )
                self.foundRE = self.bgiSelected

            if ( self.deflineType is None or
                 m is None ):
                if self.bgiNewFormat.match ( self.deflineString ) :
                    m = self.bgiNewFormat.match( self.deflineString )
                    self.foundRE = self.bgiNewFormat
                elif self.bgiOldFormat.match ( self.deflineString ) :
                    m = self.bgiOldFormat.match( self.deflineString )
                    self.foundRE = self.bgiOldFormat

                # Set here because last check may not be executed

                if m and self.saveDeflineType:
                    self.bgiSelected = self.foundRE
                    self.platform = "BGI"

            # Confirm regular expression succeeded

            if m is None :
                self.isValid = False
                return self.isValid

            # Get match values

            if self.foundRE == self.bgiOldFormat:
                (self.flowcell, self.lane, self.column, self.row, self.readNo, self.spotGroup, self.readNum, endSep) = m.groups()

                # Interpret readNum, if found

                if ( self.readNum and
                     self.readNum[0] == "/" ):
                        self.readNum = self.readNum[1:]

                # Remove # sign

                if self.spotGroup:
                    self.spotGroup = self.spotGroup[1:]

            else:
                (self.flowcell, self.lane, self.column, self.row, self.readNo, self.suffix, sep1,
                 self.readNum, self.filterRead, reserved, self.spotGroup, endSep) = m.groups()

                # Set filter value better

                if self.filterRead == 'Y':
                    self.filterRead = 1
                else:
                    self.filterRead = 0

            # Value of 0 for spot group is equivalent to no spot group

            if self.spotGroup == "0":
                self.spotGroup = ''

            # Set defline name

            self.name = f'{self.flowcell}{self.lane}{self.column}{self.row}{self.readNo}'

            if (not self.deflineType and
                self.saveDeflineType):
                self.deflineType = self.BGI
                self.platform = "BGI"
                self.sharqCompatible = True

        ############################################################
        # New Illumina
        #
        # @M00730:68:000000000-A2307:1:1101:14701:1383 1:N:0:1 (XXX574591)
        # @HWI-962:74:C0K69ACXX:8:2104:14888:94110 2:N:0:CCGATAT (XXX610048)
        # @HWI-ST808:130:H0B8YADXX:1:1101:1914:2223 1:N:0:NNNNNN-GGTCCA-AAAA (YYY1106612)
        # @HWI-M01380:63:000000000-A8KG4:1:1101:17932:1459 1:N:0:Alpha29 CTAGTACG|0|GTAAGGAG|0 (SRR1767413)
        # @HWI-ST959:56:D0AW4ACXX:8:1101:1233:2026 2:N:0: (XXX770604)
        # @DJB77P1:546:H8V5MADXX:2:1101:11528:3334 1:N:0:_I_GACGAC (WWW000015)
        # @HET-141-007:154:C391TACXX:6:1216:12924:76893 1:N:0 (WWW000006)
        # @DG7PMJN1:293:D12THACXX:2:1101:1161:1968_1:N:0:GATCAG (XXX998303)
        # @M01321:49:000000000-A6HWP:1:1101:17736:2216_1:N:0:1/M01321:49:000000000-A6HWP:1:1101:17736:2216_2:N:0:1 (WWW000016)
        # @MISEQ:36:000000000-A5BCL:1:1101:24982:8584;smpl=12;brcd=ACTTTCCCTCGA 1:N:0:ACTTTCCCTCGA (WWW000026)
        # @HWI-ST1234:33:D1019ACXX:2:1101:1415:2223/1 1:N:0:ATCACG (XXX1692309)
        # @aa,HWI-7001455:146:H97PVADXX:2:1101:1498:2093 1:Y:0:ACAAACGGAGTTCCGA (WWW000027)
        # @NS500234:97:HC75GBGXX:1:11101:6479:1067 1:N:0:ATTCAG+NTTCGC (WWW000028)
        # @M01388:38:000000000-A49F2:1:1101:14022:1748 1:N:0:0 (WWW000029)
        # @M00388:100:000000000-A98FW:1:1101:17578:2134 1:N:0:1|isu|119|c95|303 (WWW000030)
        # @HISEQ:191:H9BYTADXX:1:1101:1215:1719 1:N:0:TTAGGC##NGTCCG (WWW000031)
        # @HWI:1:X:1:1101:1298:2061 1:N:0: AGCGATAG (barcode is discarded) (WWW000039)
        # @8:1101:1486:2141 1:N:0:/1 (WWW000033)
        # @HS2000-1017_69:7:2203:18414:13643|2:N:O:GATCAG (WWW000045)
        # @HISEQ:258:C6E8AANXX:6:1101:1823:1979:CGAGCACA:1:N:0:CGAGCACA:NG:GT (SRR3156573)
        # @HWI-ST226:170:AB075UABXX:3:1101:1436:2127 1:N:0:GCCAAT (XXX104658)
        # @HISEQ06:187:C0WKBACXX:6:1101:1198:2254 1:N:0: (XXX788729)
        ############################################################

        elif ( self.deflineType in (self.ILLUMINA_NEW,
                                    self.ILLUMINA_NEW_DOUBLE,
                                    self.ILLUMINA_NEW_OLD) or
               ( self.deflineType is None and
                 ( self.illuminaNew.match(self.deflineString) or
                   self.illuminaNewWithSuffix.match(self.deflineString) or
                   self.illuminaNewNoPrefix.match(self.deflineString) or
                   self.illuminaNewWithUnderscores.match(self.deflineString) ) ) ):

            m = None
            if ( self.deflineType and
                 self.illuminaNewSelected ):
                m = self.illuminaNewSelected.match(self.deflineString)

            # Allowing for some variation in deflines for the same platform type (i.e. if not m)

            if ( self.deflineType is None or
                 m is None ):
                if self.illuminaNew.match(self.deflineString):
                    m = self.illuminaNew.match(self.deflineString)
                    self.foundRE = self.illuminaNew
                elif self.illuminaNewNoPrefix.match(self.deflineString):
                    m = self.illuminaNewNoPrefix.match(self.deflineString)
                    self.foundRE = self.illuminaNewNoPrefix
                elif self.illuminaNewWithSuffix.match(self.deflineString):
                    m = self.illuminaNewWithSuffix.match(self.deflineString)
                    self.foundRE = self.illuminaNewWithSuffix
                elif self.illuminaNewWithUnderscores.match(self.deflineString):
                    m = self.illuminaNewWithUnderscores.match(self.deflineString)
                    self.foundRE = self.illuminaNewWithUnderscores

                # Set here because last check may not be executed

                if m and self.saveDeflineType:
                    self.illuminaNewSelected = self.foundRE
                    self.platform = "ILLUMINA"

            # Confirm regular expression succeeded

            if m is None :
                self.isValid = False
                return self.isValid

            (self.prefix, sep1, self.lane, sep2, self.tile, sep3, self.x, sep4, self.y, sep5,
             self.readNum, self.filterRead, reserved, self.spotGroup, endSep ) = m.groups()

            # Capture info in sep5 as suffix

            if ( sep5 is not None and
                 len(sep5) > 2 ):
                s = self.illuminaNewSuffix.match ( sep5 )
                if s is not None:
                    self.suffix = s.group(3)

            # Value of 0 for spot group is equivalent to no spot group

            if self.spotGroup == "0":
                self.spotGroup = ''

            # Change hyphens to colons because Illumina name processing fails on this

            if sep1 == '-':
                sep1 = ':'
            if sep2 == '-':
                sep2 = ':'
            if sep3 == '-':
                sep3 = ':'
            if sep4 == '-':
                sep4 = ':'

            # Set defline name

            if self.prefix:
                self.name = f'{self.prefix}{sep1}{self.lane}{sep2}{self.tile}{sep3}{self.x}{sep4}{self.y}'
            else:
                self.name = f'{self.lane}{sep2}{self.tile}{sep3}{self.x}{sep4}{self.y}'

            # Check for doubled-up defline for both reads
            # Potential for mixed double-up and fragment deflines, too.
            # So setting to self.ILLUMINA_NEW_DOUBLE on first occurrence
            # Assuming single character separator between the two names.

            if ( self.deflineType == self.ILLUMINA_NEW_DOUBLE or
                 ( not self.deflineType and
                   len(self.spotGroup) > len(self.name) and
                   re.search( re.escape(self.name), self.spotGroup) ) ):
                start = str.find(self.spotGroup,self.name)
                if start != -1:
                    self.spotGroup = self.spotGroup[0:start-1]
                    if self.saveDeflineType:
                        self.deflineType = self.ILLUMINA_NEW_DOUBLE

            # Check for and remove read numbers after spot group

            elif ( self.deflineType == self.ILLUMINA_NEW_OLD or
                   ( not self.deflineType and
                     self.readNumberAtEnd.search(self.spotGroup) ) ):
                self.spotGroup = self.spotGroup[0:len(self.spotGroup)-2]
                if self.saveDeflineType:
                    self.deflineType = self.ILLUMINA_NEW_OLD

            # Set filter value better

            if self.filterRead == 'Y':
                self.filterRead = 1
            else:
                self.filterRead = 0

            # Save defline type if not previously set (must be here)

            if ( not self.deflineType and
                 self.saveDeflineType ):
                self.deflineType = self.ILLUMINA_NEW
                self.sharqCompatible = True

        ############################################################
        # Old Illumina
        #
        # @HWIEAS210R_0014:3:28:1070:11942#ATCCCG/1 (XXX799414)
        # @FCABM5A:1:1101:15563:1926#CTAGCGCT_CATTGCTT/1 (WWW000014)
        # @HWI-ST1374:1:1101:1161:2060#0/1 (XXX001101)
        # @ID57_120908_30E4FAAXX:3:1:1772:953/1 (XXX020199)
        # @HWI-EAS30_2_FC200HWAAXX_4_1_646_399^M (XXX013586)
        # @R1:1:221:197:649/1 (YYY003002)
        # @R16:8:1:0:1617#0/1^M (YYY003003)
        # @7:1:792:533 (YYY013547)
        # @7:1:164:-0 (YYY013547)
        # @HWUSI-BETA8_3:7:1:-1:14 (YYY015261)
        # @rumen9533:0:0:0:5 (YYY020795)
        # @IL10_334:1:1:4:606 (YYY036999)
        # @HWUSI-EAS517-74:3:1:1023:8178/1 ~ RGR:Uk6; (WWW000002)
        # @FCC19K2ACXX:3:1101:1485:2170#/1 (WWW000017)
        # @MB27-02-1:1101:1723:2171 /1 (WWW000003)
        # @AMS2007273_SHEN-MISEQ01:47:1:1:12958:1771:0:1#0 (WWW000008)
        # @AMS2007273_SHEN-MISEQ01:47:1:1:17538:1769:0:0#0 (WWW000008)
        # @120315_SN711_0193_BC0KW7ACXX:1:1101:1419:2074:1#0/1 (WWW000001)
        # @ATGCT_7_1101_1418_2098_1 (WWW000007)
        # HWI-ST155_0544:1:1:6804:2058#0/1 (SINGLE LINE FASTQ SO NO '@') (YYY003101)
        # @HWI-ST225:626:C2Y82ACXX:3:1208:2931:82861_1 (WWW000040)
        # >M01056:83:000000000-A7GBN:1:1108:17094:2684--W1 (WWW000019)
        # @D3LH75P1:1:1101:1054:2148:0 1:1 (WWW000032)
        # @HWI-IT879:92:5:1101:1170:2026#0/1:0 (WWW000037)
        # @HWI-ST225:626:C2Y82ACXX:3:1208:2222:82880_1 (WWW000040)
        # @FCA5PJ4:1:1101:14707:1407#GTAGTCGC_AGCTCGGT/1 (SRR2006030)
        # @SOLEXA-GA02_1:1:1:0:106 (ERR011021)
        # @HWUSI-EAS499:1:3:9:1822#0/1 (XXX000093)
        # @BILLIEHOLIDAY_1_FC20F3DAAXX:8:2:342:540 (XXX013565)
        # @BILLIEHOLIDAY_1_FC200TYAAXX_3_1_751_675 (XXX013571)
        # @HWI-EAS30_2_FC20416AAXX_7_1_116_317 (XXX013600)
        # >KN-930:1:1:653:356 (XXX014056)
        # USI-EAS50_1:6:1:392:881 (XXX014283)
        # @FC12044_91407_8_1_46_673 (XXX015015)
        # @HWI-EAS299_2_30MNAAAXX:5:1:936:1505/1 (XXX015076)
        # @HWI-EAS-249:7:1:1:443/1 (XXX037956)
        # @741:6:1:1204:10747/1 (XXX094419)
        # @HWI-EAS385_0086_FC:1:1:1239:943#0/1 (XXX488373)
        # @HWUSI-EAS613-R_0001:8:1:1020:14660#0/1 (XXX556206)
        # @ILLUMINA-D01686_0001:7:1:1028:14175#0/1 (XXX567550)
        # @1920:1:1:1504:1082/1 (XXX627950)
        # @HWI-EAS397_0013:1:1:1083:11725#0/1 (XXX651965)
        # >HWI-EAS6_4_FC2010T:1:1:80:366 (YYY001656)
        # @NUTELLA_42A08AAXX:4:001:0003:0089/1 (YYY003100)
        # @R16:8:1:0:875#0/1 (YYY014126)
        # HWI-EAS102_1_30LWPAAXX:5:1:1456:776 (YYY016872)
        # @HWI-EAS390_30VGNAAXX1:1:1:377:1113/1 (YYY020188)
        # @ID57_120908_30E4FAAXX:3:1:1772:953/1 (YYY020203)
        # HWI-EAS440_102:8:1:168:1332 (YYY029167)
        # @SNPSTER4_246_30GCDAAXX_PE:1:1:3:896/1 (YYY029194)
        # @FC42AUBAAXX:6:1:4:1280#TGACCA/1 (YYY030833)
        # @SOLEXA9:1:1:1:2005#0/1 (YYY037749)
        # @SNPSTER3_264_30JGGAAXX_PE:2:1:218:311/1 (YYY058403)
        # @HWUSI-EAS535_0001:7:1:747:14018#0/1 (YYY065453)
        # @FC42ATTAAXX:5:1:0:20481 (YYY066636)
        # @HWUSI-EAS1571_0012:8:1:1017:20197#0/1 (YYY089777)
        # @SOLEXA1_0052_FC:8:1:1508:1078#TTAGGC/1 (YYY171628)
        # @SN971:2:1101:15.80:103.70#0/1 (WWW000034)
        ############################################################

        elif ( self.deflineType == self.ILLUMINA_OLD or
               ( self.deflineType is None and
                 ( self.illuminaOldColon.match ( self.deflineString ) or
                   self.illuminaOldUnderscore.match ( self.deflineString ) or
                   self.illuminaOldNoPrefix.match ( self.deflineString ) or
                   self.illuminaOldWithSuffix.match ( self.deflineString ) or
                   self.illuminaOldWithSuffix2.match ( self.deflineString ) ) ) ):

            # For first time around check for need to identify appropriate separator
            # Retain defline type if desired and not set and count extra numbers in name
            # Retain appropriate pattern. Note that illumina old with junk comes first
            # because illumina old colon pattern will always match that case, too.

            m = None
            if ( self.deflineType and
                 self.illuminaOldSelected ):
                m = self.illuminaOldSelected.match ( self.deflineString )

            # Allowing for some variation (i.e. check for m is None)

            if ( self.deflineType is None or
                 m is None ):
                if self.illuminaOldWithSuffix.match ( self.deflineString ):
                    m = self.illuminaOldWithSuffix.match ( self.deflineString )
                    self.foundRE = self.illuminaOldWithSuffix
                elif self.illuminaOldColon.match ( self.deflineString ) :
                    m = self.illuminaOldColon.match ( self.deflineString )
                    self.foundRE = self.illuminaOldColon
                elif self.illuminaOldUnderscore.match ( self.deflineString ):
                    m = self.illuminaOldUnderscore.match ( self.deflineString )
                    self.foundRE = self.illuminaOldUnderscore
                elif self.illuminaOldWithSuffix2.match ( self.deflineString ):
                    m = self.illuminaOldWithSuffix2.match ( self.deflineString )
                    self.foundRE = self.illuminaOldWithSuffix2
                elif self.illuminaOldNoPrefix.match ( self.deflineString ):
                    m = self.illuminaOldNoPrefix.match ( self.deflineString )
                    self.foundRE = self.illuminaOldNoPrefix

                # Set here because last check may not be executed

                if m and self.saveDeflineType:
                    self.illuminaOldSelected = self.foundRE
                    self.platform = "ILLUMINA"

            # Confirm regular expression succeeded

            if m is None :
                self.isValid = False
                return self.isValid

            # Collect values from regular expression pattern

            (self.prefix, sep1, self.lane, sep2, self.tile, sep3, self.x, sep4, self.y,
             self.spotGroup, self.readNum, endSep) = m.groups()
            if self.readNum:
                self.readNum = self.readNum[1:]

            # Check for suffix

            # s2 = self.illuminaOldSuffix2.match ( self.y )
            # if s2:
            #     (self.y, self.filterRead, reserved, self.spotGroup) = s2.groups()
            # else:
            s = self.illuminaOldSuffix.match ( self.y )
            if s:
                ( self.y, self.suffix ) = s.groups()
                if len(self.suffix) < 3:
                    self.suffix = None
                elif self.readNumberInSuffix.search(self.suffix):
                    self.suffix = self.suffix[2:]

            # Determine number of discards first time through (retain float values in suffix if present)

            if ( not self.deflineType and
                 self.prefix) :
                self.numDiscards = self.countExtraNumbersInIllumina(sep4)

            if ( self.numDiscards == 2 and
                 self.floatInDiscards.search(self.x) ):
                if self.suffix:
                    self.suffix = f'{sep3}{self.x}{sep4}{self.y}{self.suffix}'
                else:
                    self.suffix = f'{sep3}{self.x}{sep4}{self.y}'
            elif ( self.numDiscards == 1 and
                   self.floatInDiscards.search(self.y) ):
                if self.suffix:
                    self.suffix = sep4 + self.y + self.suffix
                else:
                    self.suffix = sep4 + self.y

            # Change hyphens to colons because Illumina name processing fails on this

            if sep1 == '-':
                sep1 = ':'
            if sep2 == '-':
                sep2 = ':'
            if sep3 == '-':
                sep3 = ':'
            if sep4 == '-':
                sep4 = ':'

            # Discard extra numbers (which can cause a difference between pair deflines)
            # and set defline name

            if self.numDiscards > 0:
                if self.numDiscards == 1:
                    self.y = self.x
                    self.x = self.tile
                    self.tile = self.lane
                    m2 = re.match("([!-~]*?)(" + sep4 + ")(\d+)$", self.prefix)
                    if m2:
                        (self.prefix,sep0,self.lane) = m2.groups()
                        self.name = f'{self.prefix}{sep0}{self.lane}{sep1}{self.tile}{sep2}{self.x}{sep3}{self.y}'
                    else:
                        self.lane = self.prefix
                        self.prefix = ""
                        self.name = f'{self.lane}{sep1}{self.tile}{sep2}{self.x}{sep3}{self.y}'

                elif self.numDiscards == 2:
                    self.y = self.tile
                    self.x = self.lane
                    m2 = re.match("([!-~]*?)(" + sep4 + ")(\d+)(" + sep4 + ")(\d+)(\s+|$)",self.prefix)
                    if m2:
                        (self.prefix,sep_1,self.lane,sep0,self.tile,sepUnused) = m2.groups()
                        self.name = f'{self.prefix}{sep_1}{self.lane}{sep0}{self.tile}{sep1}{self.x}{sep2}{self.y}'
                    else:
                        m2 = re.match("(\d+)(" + sep4 + ")(\d+)(\s+|$)",self.prefix)
                        if m2:
                            (self.lane,sep0,self.tile,sepUnused) = m2.groups()
                            self.prefix = ""
                            self.name = f'{self.lane}{sep0}{self.tile}{sep1}{self.x}{sep2}{self.y}'

            # If no discards and prefix exists, set name including prefix

            elif self.prefix:
                self.name = f'{self.prefix}{sep1}{self.lane}{sep2}{self.tile}{sep3}{self.x}{sep4}{self.y}'

            # If no discards and no prefix, set name without prefix

            else:
                self.name = f'{self.lane}{sep2}{self.tile}{sep3}{self.x}{sep4}{self.y}'

            # Remove # from front of spot group if present
            # Value of 0 for spot group is equivalent to no spot group

            if self.spotGroup:
                self.spotGroup = self.spotGroup[1:]

            if self.spotGroup == "0":
                self.spotGroup = ''

            # Save defline type and regular expression (must occur after determination of numDiscards)

            if ( not self.deflineType and
                 self.saveDeflineType ):
                self.deflineType = self.ILLUMINA_OLD

        ############################################################
        # qiime with new illumina (only first spot qiime name retained)
        #
        # @B11.13210.SIV.Barouch.Stool.250.06.8.13.12_5644 M00181:229:000000000-AAPUA:1:1101:9433:3327 1:N:0:1 orig_bc=TGACCTCCTAGA new_bc=TGACCTCCAAGA bc_diffs=1 (YYY908068)
        # @2wkRT.79_123 M00176:18:000000000-A0DK4:1:1:13923:1732 1:N:0:0 orig_bc=ATGCTAACCACG new_bc=ATGCTAACCACG bc_diffs=0 (XXX776282)
        # @ HWI-M01929:28:000000000-A6VG4:1:2109:8848:7133 1:N:0:GTGTT  orig_bc=GCTTA   new_bc=GCTTA    bc_diffs=0 (WWW000005)
        # @cp2  :1:1101:17436:1559:1:N:0:5/1_:1:1101:17436:1559:2:N:0:5/2   124 124 (WWW000020)
        # @10_194156 HWI-M02808:46:AAHRM:1:1101:17744:1823 1:N:0:ATGAGACTCCAC orig_bc=ATGAGACTCCAC new_bc=ATGAGACTCCAC bc_diffs=0 (WWW000013)
        # @AM-B-CON M02233:62:000000000-A9GLW:1:1101:15425:1859 1:N:0:111^M (WWW000024)
        # >26.04.2015.WO.Comp.S55_1 M00596.112.000000000.AHGJM.1.1101.20901.1309 1.N.0.55 (SRR3112744)
        # >PSF.1d.20_0 HISEQ:128:160215_SNL128_0128_AHJT2MBCXX:1:1101:3422:2184 1:N:0: orig_bc=TATAGCGACTACTATA new_bc=TATAGCGACTACTATA bc_diffs=0 (SRR3992252)
        # @2-796964 M01929:5:000000000-A46YE:1:1108:16489:18207 1:N:0:2 (XXX1778155)
        ############################################################

        elif ( self.deflineType in (self.QIIME_ILLUMINA_NEW,
                                    self.QIIME_ILLUMINA_NEW_BC,
                                    self.QIIME_ILLUMINA_NEW_DBL,
                                    self.QIIME_ILLUMINA_NEW_DBL_BC) or
               ( self.deflineType is None and
                 ( self.qiimeIlluminaNew.match(self.deflineString) or
                   self.qiimeIlluminaNew2.match(self.deflineString) or
                   self.qiimeIlluminaNewPeriods.match(self.deflineString) or
                   self.qiimeIlluminaNewUnderscores.match(self.deflineString) ) ) ):

            m = None
            if ( self.deflineType and
                 self.qiimeSelected ):
                m = self.qiimeSelected.match ( self.deflineString )

            # Allowing for some variation (i.e. check for m is None)

            if ( self.deflineType is None or
                 m is None ):
                if self.qiimeIlluminaNew.match ( self.deflineString ) :
                    m = self.qiimeIlluminaNew.match ( self.deflineString )
                    self.foundRE = self.qiimeIlluminaNew
                elif self.qiimeIlluminaNew2.match ( self.deflineString ) :
                    m = self.qiimeIlluminaNew2.match ( self.deflineString )
                    self.foundRE = self.qiimeIlluminaNew2
                elif self.qiimeIlluminaNewPeriods.match ( self.deflineString ):
                    m = self.qiimeIlluminaNewPeriods.match ( self.deflineString )
                    self.foundRE = self.qiimeIlluminaNewPeriods
                elif self.qiimeIlluminaNewUnderscores.match ( self.deflineString ):
                    m = self.qiimeIlluminaNewUnderscores.match ( self.deflineString )
                    self.foundRE = self.qiimeIlluminaNewUnderscores

                # Set here because last check may not be executed

                if m and self.saveDeflineType:
                    self.qiimeSelected = self.foundRE
                    self.platform = "ILLUMINA"

            # Confirm regular expression succeeded

            if m is None :
                self.isValid = False
                return self.isValid

            # Get match values

            if self.qiimeSelected == self.qiimeIlluminaNew2:
                (self.qiimeName, qiimeName2, self.prefix, sep1, self.lane, sep2, self.tile, sep3, self.x, sep4, self.y, sep5,
                self.readNum, self.filterRead, reserved, self.spotGroup, endSep ) = m.groups()
                self.qiimeName = self.qiimeName + '_' + qiimeName2
            else:
                (self.qiimeName, self.prefix, sep1, self.lane, sep2, self.tile, sep3, self.x, sep4, self.y, sep5,
                self.readNum, self.filterRead, reserved, self.spotGroup, endSep ) = m.groups()

            # Check for the presence of barcode corrections

            m_bc = None
            if ( self.deflineType == self.QIIME_ILLUMINA_NEW_BC or
                 ( not self.deflineType and
                   self.qiimeBc.match(self.deflineString) ) ) :
                m_bc = self.qiimeBc.match(self.deflineString)
                if ( self.deflineType == self.QIIME_ILLUMINA_NEW_BC and
                     m_bc is None ):
                    self.isValid = False
                    return self.isValid
                else:
                    self.spotGroup = m_bc.group(2)

            # Set defline name

            self.name = f'{self.prefix}{sep1}{self.lane}{sep2}{self.tile}{sep3}{self.x}{sep4}{self.y}'

            # Check for doubled-up defline for both reads
            # Potential for mixed double-up and fragment deflines, too.
            # So setting to self.ILLUMINA_NEW_DOUBLE on first occurrence
            # Assuming single character separator between the two names.

            if ( self.deflineType in (self.QIIME_ILLUMINA_NEW_DBL,
                                      self.QIIME_ILLUMINA_NEW_DBL_BC) or
                 ( not self.deflineType and
                   len(self.spotGroup) > len(self.name) and
                   re.search( re.escape(self.name), self.spotGroup) ) ):
                start = str.find(self.spotGroup,self.name)
                if start != -1:
                    self.spotGroup = self.spotGroup[0:start-1]
                    if self.readNumInSpotGroup.search(self.spotGroup):
                        self.spotGroup = self.spotGroup[:len(self.spotGroup) - 2]
                    if ( not self.deflineType and
                         self.saveDeflineType ):
                        if self.deflineType == self.QIIME_ILLUMINA_NEW_BC:
                            self.deflineType = self.QIIME_ILLUMINA_NEW_DBL_BC
                        else:
                            self.deflineType = self.QIIME_ILLUMINA_NEW_DBL

            # Set filter value better

            if self.filterRead == 'Y':
                self.filterRead = 1
            else:
                self.filterRead = 0

            # Retain defline type if desired (must stay here)

            if ( not self.deflineType and
                 self.saveDeflineType ):
                if m_bc is None:
                    self.deflineType = self.QIIME_ILLUMINA_NEW
                else:
                    self.deflineType = self.QIIME_ILLUMINA_NEW_BC

        ############################################################
        # qiime with old illumina (only first spot qiime name retained)
        #
        # @B11.13210.SIV.Barouch.Stool.250.06.8.13.12_378 M00181:229:000000000-AAPUA:1:1101:19450:2192#0/1 orig_bc=TGACCTCCAAGA new_bc=TGACCTCCAAGA bc_diffs=0 (ZZZ908068)
        # @B11.13210.SIV.Barouch.Stool.250.06.8.13.12_158 M00181:229:000000000-AAPUA:1:1101:19450:2192#0/2 orig_bc=TGACCTCCAAGA new_bc=TGACCTCCAAGA bc_diffs=0 (ZZZ908068)
        ############################################################

        elif ( self.deflineType in (self.QIIME_ILLUMINA_OLD,
                                    self.QIIME_ILLUMINA_OLD_BC) or
               ( self.deflineType is None and
                 self.qiimeIlluminaOld.match(self.deflineString) ) ):

            m = self.qiimeIlluminaOld.match(self.deflineString)

            # Confirm regular expression succeeded

            if m is None :
                self.isValid = False
                return self.isValid

            # Get match values

            (self.qiimeName, self.prefix, sep1, self.lane, sep2, self.tile, sep3, self.x, sep4, self.y,
             self.spotGroup, self.readNum, endSep) = m.groups()
#            self.qiimeName = self.qiimeName.replace(" ","_")
            if self.readNum:
                self.readNum = self.readNum[1:]
            if self.spotGroup:
                self.spotGroup = self.spotGroup[1:]

            # Check for the presence of barcode corrections

            m_bc = None
            if ( self.deflineType == self.QIIME_ILLUMINA_OLD_BC or
                 ( not self.deflineType and
                   self.qiimeBc.match(self.deflineString) ) ) :
                m_bc = self.qiimeBc.match(self.deflineString)
                if ( self.deflineType == self.QIIME_ILLUMINA_OLD_BC and
                     m_bc is None ):
                    self.isValid = False
                    return self.isValid
                else:
                    self.spotGroup = m_bc.group(2)

            # Set defline name

            self.name = f'{self.prefix}{sep1}{self.lane}{sep2}{self.tile}{sep3}{self.x}{sep4}{self.y}'

            # Retain defline type if desired

            if ( not self.deflineType and
                 self.saveDeflineType ):
                if m_bc is None:
                    self.deflineType = self.QIIME_ILLUMINA_OLD
                else:
                    self.deflineType = self.QIIME_ILLUMINA_OLD_BC
                self.platform = "ILLUMINA"

        ############################################################
        # 454 defline
        #
        # @GG3IVWD03F5DLB length=97 xy=2404_1917 region=3 run=R_2010_05_11_11_15_22_ (XXX529889)
        # @EV5T11R03G54ZJ (YYY307780)
        # @GKW2OSF01D55D9 (ERR016499)
        # @OS-230b_GLZVSPV04JTNWT (ERR039808)
        # @HUT5UCF07H984F (SRR2035362)
        # @EM7LVYS02FOYNU/1 (WWW000042 or WWW000043)
        ############################################################

        elif ( self.deflineType == self.LS454 or
               ( self.deflineType is None and
                 self.ls454.match( self.deflineString ) ) ):

            # Capture 454 values

            m = self.ls454.match( self.deflineString )

            # Confirm regular expression succeeded

            if m is None :
                self.isValid = False
                return self.isValid

            # Get match values

            (self.prefix,self.dateAndHash454,self.region454,self.xy454,self.readNum,endSep) = m.groups()
            if self.readNum:
                self.readNum = self.readNum[1:]

            # Set name

            self.name = f'{self.prefix}{self.dateAndHash454}{self.region454}{self.xy454}'

            # Retain defline type if desired

            if ( not self.deflineType and
                 self.saveDeflineType ):
                self.deflineType = self.LS454
                self.platform = "LS454"

        ############################################################
        # qiime with 454
        #
        # @T562_7000012 H29C5KU01AZBDB orig_bc=AGCTCACGTA new_bc=AGCTCACGTA bc_diffs=0 (WWW000018)
        ############################################################

        elif ( self.deflineType in (self.QIIME_454,
                                    self.QIIME_454_BC) or
               ( self.deflineType is None and
                 self.qiime454.match( self.deflineString ) ) ):

            # Capture 454 values

            m = self.qiime454.match( self.deflineString )

            # Confirm regular expression succeeded

            if m is None :
                self.isValid = False
                return self.isValid

            # Get match values

            (self.qiimeName,self.dateAndHash454,self.region454,self.xy454,self.readNum,endSep) = m.groups()
            if self.readNum:
                self.readNum = self.readNum[1:]

            # Extract barcode if present

            m_bc = None
            if ( self.deflineType == self.QIIME_454_BC or
                 ( not self.deflineType and
                   self.qiimeBc.match(self.deflineString) ) ) :
                m_bc = self.qiimeBc.match(self.deflineString)
                if ( self.deflineType == self.QIIME_454_BC and
                     m_bc is None ):
                    self.isValid = False
                    return self.isValid
                else:
                    self.spotGroup = m_bc.group(2)

            # Set name

            self.name = f'{self.dateAndHash454}{self.region454}{self.xy454}'

            # Retain defline type if desired

            if ( not self.deflineType and
                 self.saveDeflineType ):
                if m_bc is None:
                    self.deflineType = self.QIIME_454
                else:
                    self.deflineType = self.QIIME_454_BC
                self.platform = "LS454"

        ############################################################
        # Pacbio CCS/RoIs reads or subreads
        #
        # See https://www.biostars.org/p/146048/ for field descriptions
        # Also see https://speakerdeck.com/pacbio/specifics-of-smrt-sequencing-data
        #
        # @m120525_202528_42132_c100323432550000001523017609061234_s1_p0/43 (XXX941211)
        # @m101111_134728_richard_c000027022550000000115022502211150_s1_p0/1 (YYY075011)
        # @m110115_082846_Uni_c000000000000000000000012706400001_s3_p0/1/0_508 (SRR497981)
        # @m130727_043304_42150_c100538232550000001823086511101337_s1_p0/16/0_5273 (XXX989791)
        # @m120328_022709_00128_c100311312550000001523011808061260_s1_p0/129/0_4701 (WWW000010)
        # @m120204_011539_00128_c100220982555400000315052304111230_s2_p0/8/0_1446 (WWW000009)
        # @m54261_181223_050738_4194378 (SRR9066822)
        # @m64049_200827_171349/2/ccs (SRR13103091)
        # @m54336U_191028_160945/1/120989_128438 (SRR14585543)
        ############################################################

        elif ( self.deflineType == self.PACBIO or
               ( self.deflineType is None and
                 ( self.pacbio.match ( self.deflineString ) or
                   self.pacbio2.match ( self.deflineString ) or
                   self.pacbio3.match ( self.deflineString ) or
                   self.pacbio4.match ( self.deflineString ) ) ) ):

            # Capture name after '@' sign

            m = None
            if ( self.deflineType and
                 self.pacbioSelected ) :
                m = self.pacbioSelected.match ( self.deflineString )

            if ( self.deflineType is None or
                 m is None ):
                if self.pacbio.match ( self.deflineString ) :
                    m = self.pacbio.match( self.deflineString )
                    self.foundRE = self.pacbio
                elif self.pacbio2.match ( self.deflineString ) :
                    m = self.pacbio2.match( self.deflineString )
                    self.foundRE = self.pacbio2
                elif self.pacbio3.match ( self.deflineString ) :
                    m = self.pacbio3.match( self.deflineString )
                    self.foundRE = self.pacbio3
                elif self.pacbio4.match ( self.deflineString ) :
                    m = self.pacbio4.match( self.deflineString )
                    self.foundRE = self.pacbio4

                # Set here because last check may not be executed

                if m and self.saveDeflineType:
                    self.pacbioSelected = self.foundRE
                    self.platform = "PACBIO"

            # Confirm regular expression succeeded

            if m is None :
                self.isValid = False
                return self.isValid

            # Get match value

            self.name = m.group(1)

            # Retain defline type if desired

            if ( not self.deflineType and
                 self.saveDeflineType ):
                self.deflineType = self.PACBIO
                self.platform = "PACBIO"

        ############################################################
        # Old illumina bar code and/or read number only
        # Must occur after pacbio defline check
        #
        # @_2_#GATCAGAT/1 (WWW000004)
        # @Read_190546#BC005 length=1419 (XXX1616052)
        ############################################################

        elif ( self.deflineType == self.ILLUMINA_OLD_BC_RN or
               ( self.deflineType is None and
                 self.retainAltNum is False and
                 ( self.illuminaOldBcRnOnly.match ( self.deflineString ) or
                   self.illuminaOldBcOnly.match ( self.deflineString ) or
                   self.illuminaOldRnOnly.match ( self.deflineString ) ) ) ):

            m = None
            if ( self.deflineType and
                 self.illuminaOldBcRnSelected ):
                m = self.illuminaOldBcRnSelected.match ( self.deflineString )

            if ( self.deflineType is None or
                 m is None ):
                if self.illuminaOldBcRnOnly.match ( self.deflineString ) :
                    m = self.illuminaOldBcRnOnly.match ( self.deflineString )
                    self.foundRE = self.illuminaOldBcRnOnly
                elif self.illuminaOldBcOnly.match ( self.deflineString ):
                    m = self.illuminaOldBcOnly.match ( self.deflineString )
                    self.foundRE = self.illuminaOldBcOnly
                elif self.illuminaOldRnOnly.match ( self.deflineString ):
                    m = self.illuminaOldRnOnly.match ( self.deflineString )
                    self.foundRE = self.illuminaOldRnOnly

                if m and self.saveDeflineType:
                    self.illuminaOldBcRnSelected = self.foundRE

            # Confirm regular expression succeeded

            if m is None :
                self.isValid = False
                return self.isValid

            # Assign values

            (self.name, self.spotGroup, self.readNum, endSep) = m.groups()
            if self.spotGroup[0] == "#":
                self.spotGroup = self.spotGroup[1:]
                if ( self.readNum[0] == "/" or
                     self.readNum[0:1] == "\\" ):
                    self.readNum = self.readNum[1:]
                else:
                    self.readNum = None
            else:
                self.readNum = self.spotGroup[1:]
                self.spotGroup = None

            if self.spotGroup == "0":
                self.spotGroup = None

            # Save defline type

            if ( not self.deflineType and
                 self.saveDeflineType) :
                self.deflineType = self.ILLUMINA_OLD_BC_RN

            # Switch to ALTNUM if readNum values are getting large

            elif self.readNum == '5':
                if self.saveDeflineType:
                    self.deflineType = self.ALTNUM
                self.altNum = self.readNum
                self.retainAltNum = True
                self.readNum = ''

        ############################################################
        # New illumina data group
        #
        # @V300033286L2C001R0010000001 1:N:0:TCTGCACG
        ############################################################

        elif ( self.deflineType == self.ILLUMINA_NEW_DATA_GRP or
               ( self.deflineType is None and
                 self.illuminaNewDataGroup.match(self.deflineString) ) ):

            m = self.illuminaNewDataGroup.match(self.deflineString)

            # Confirm regular expression succeeded

            if m is None :
                self.isValid = False
                return self.isValid

            # Get match values

            (self.name, sep1, self.readNum, self.filterRead, reserved, self.spotGroup, endSep) = m.groups()

            # Set filter value better

            if self.filterRead == 'Y':
                self.filterRead = 1
            else:
                self.filterRead = 0

            # Spot group of zero changed to ''

            if self.spotGroup == '0':
                self.spotGroup = ''

            # Retain defline type if desired

            if ( not self.deflineType and
                 self.saveDeflineType ):
                self.deflineType = self.ILLUMINA_NEW_DATA_GRP
                self.sharqCompatible = True

        ############################################################
        # ion torrent
        #
        # @A313D:7:49 (XXX486160)
        # @RD4FE:00027:00172 (XXX2925654)
        # @ONBWR:00329:02356/1 (WWW000044)
        # >311CX:3560:2667   length=347 (SRR547526)
        #
        ############################################################

        elif ( self.deflineType == self.ION_TORRENT or
               ( self.deflineType is None and
                 ( self.ionTorrent.match ( self.deflineString ) or
                   self.ionTorrent2.match ( self.deflineString ) ) ) ):

            # Capture name after '@' sign

            m = None
            if ( self.deflineType and
                 self.ionTorrentSelected ) :
                m = self.ionTorrentSelected.match ( self.deflineString )
                self.foundRE = self.ionTorrentSelected

            if ( self.deflineType is None or
                 m is None ):
                if self.ionTorrent2.match ( self.deflineString ) :
                    m = self.ionTorrent2.match( self.deflineString )
                    self.foundRE = self.ionTorrent2
                elif self.ionTorrent.match ( self.deflineString ) :
                    m = self.ionTorrent.match( self.deflineString )
                    self.foundRE = self.ionTorrent

                # Set here because last check may not be executed

                if m and self.saveDeflineType:
                    self.ionTorrentSelected = self.foundRE
                    self.platform = "ION_TORRENT"

            # Confirm regular expression succeeded

            if m is None :
                self.isValid = False
                return self.isValid

            # Get match values

            if self.foundRE == self.ionTorrent:
                (self.runId, sep1, self.row, sep2, self.column, self.readNum, endSep) = m.groups()

                # Interpret readNum, if found

                if self.readNum:
                    if (self.readNum[0] == "/" or
                            self.readNum[0:1] == "\\"):
                        self.readNum = self.readNum[1:]
                    elif self.readNum == "L":
                        self.readNum = "1"
                    elif self.readNum == "R":
                        self.readNum = "2"
            else:
                (self.runId, sep1, self.row, sep2, self.column, sep3,
                 self.readNum, self.filterRead, reserved, self.spotGroup, endSep) = m.groups()

            # Set defline name

            self.name = f'{self.runId}{sep1}{self.row}{sep2}{self.column}'

            if (not self.deflineType and
                self.saveDeflineType):
                self.deflineType = self.ION_TORRENT
                self.platform = "ION_TORRENT"

        ############################################################
        # generic qiime
        #
        # @5SP.pna.R_1030 METII:00037:00101 orig_bc=TTCCTACCAGTC new_bc=TTCCTACCAGTC bc_diffs=0 (ERR4604646)
        # @10317.000016458_0 orig_bc=TGCACCTCTGTC new_bc=TGCACCTCTGTC bc_diffs=0 (WWW000022)
        ############################################################

        elif ( self.deflineType == self.QIIME_GENERIC or
               ( self.deflineType is None and
                 ( self.qiimeBc2.match( self.deflineString ) or
                   self.qiimeBc.match( self.deflineString ) ) ) ):

            # Capture generic qiime values

            m = None
            if ( self.deflineType and
                 self.qiimeSelected ):
                m = self.qiimeSelected.match ( self.deflineString )
                self.foundRE = self.qiimeSelected

            # Allowing for some variation (i.e. check for m is None)

            if ( self.deflineType is None or
                 m is None ):
                if self.qiimeBc2.match( self.deflineString ):
                    m = self.qiimeBc2.match( self.deflineString )
                    self.foundRE = self.qiimeBc2
                elif self.qiimeBc.match( self.deflineString ):
                    m = self.qiimeBc.match( self.deflineString )
                    self.foundRE = self.qiimeBc

                # Set here because last check may not be executed

                if m and self.saveDeflineType:
                    self.qiimeSelected = self.foundRE
                    self.platform = "UNDEFINED"

            # Confirm regular expression succeeded

            if m is None :
                self.isValid = False
                return self.isValid

            # Get match values

            if self.foundRE == self.qiimeBc:
                (self.name,self.spotGroup) = m.groups()
            else:
                (self.qiimeName,self.name,self.spotGroup) = m.groups()

            # Retain defline type if desired

            if ( not self.deflineType and
                 self.saveDeflineType ):
                self.deflineType = self.QIIME_GENERIC
                self.platform = "UNDEFINED"

        ############################################################
        # Nanopore/MinION fastq
        #
        # @77_2_1650_1_ch100_file0_strand_twodirections (XXX2761339)
        # @77_2_1650_1_ch100_file16_strand_twodirections:pass\77_2_1650_1_ch100_file16_strand.fast5 (SRR2761339)
        # @channel_108_read_11_twodirections:flowcell_17/LomanLabz_PC_E.coli_MG1655_ONI_3058_1_ch108_file21_strand.fast5 (XXX637417 or YYY637417)
        # @channel_108_read_11_complement:flowcell_17/LomanLabz_PC_E.coli_MG1655_ONI_3058_1_ch108_file21_strand.fast5 (XXX637417 or YYY637417)
        # @channel_108_read_11_template:flowcell_17/LomanLabz_PC_E.coli_MG1655_ONI_3058_1_ch108_file21_strand.fast5 (XXX637417 or YYY637417)
        # @channel_346_read_183-1D (SRR1747417)
        # @channel_346_read_183-complement (SRR1747417)
        # @channel_346_read_183-2D (SRR1747417)
        # @ch120_file13-1D (SRR1980822)
        # @ch120_file13-2D (SRR1980822)
        # @channel_108_read_8:LomanLabz_PC_E.coli_MG1655_ONI_3058_1_ch108_file18_strand.fast5 (R-based poRe fastq requires filename, too) (WWW000025)
        # @1dc51069-f61f-45db-b624-56c857c4e2a8_Basecall_2D_000_2d oxford_PC_HG02.3attempt_0633_1_ch96_file81_strand_twodirections:CORNELL_Oxford_Nanopore/oxford_PC_HG02.3attempt_0633_1_ch96_file81_strand.fast5 (SRR2848544 - nanopore3)
        # @ae74c4fb-2c1d-4176-9584-3dfcc6dce41e_Basecall_2D_2d UT317077_20160808_FNFAD22478_MN19846_sequencing_run_FHV_Barcoded_TakeII_88358_ch93_read2620_strand NB06\UT317077_20160808_FNFAD22478_MN19846_sequencing_run_FHV_Barcoded_TakeII_88358_ch93_read2620_strand.fast5 (SRR5085901 - nanopore3)
        # @ddb7d987-73c0-4d9a-8ac0-ac0dbc462ab5_Basecall_2D_2d UT317077_20160808_FNFAD22478_MN19846_sequencing_run_FHV_Barcoded_TakeII_88358_ch100_read4767_strand1 NB06\UT317077_20160808_FNFAD22478_MN19846_sequencing_run_FHV_Barcoded_TakeII_88358_ch100_read4767_strand1.fast5 (SRR5085901 - nanopore3)
        # @f286a4e1-fb27-4ee7-adb8-60c863e55dbb_Basecall_Alignment_template MINICOL235_20170120_FN__MN16250_sequencing_throughput_ONLL3135_25304_ch143_read16010_strand
        # @channel_101_read_1.1C|1T|2D (ERR1121618 bam converted to fastq)
        # @channel_100_read_20_twodirections:/Users/blbrown/Documents/DATA/Biology Stuff/New Building/Oxford Nanopore/Pan-MAP/VCUEGLequalFAA23773/reads/downloads/pass/Bonnie_PC_VCUEGLequalFAA23773_2353_1_ch100_file20_strand.fast5 (SRR3473970)
        # @5f8415e3-46ae-48fc-9092-a291b8b6a9b9 run_id=47b8d024d71eef532d676f4aa32d8867a259fc1b read=279 mux=3 ch=87 start_time=2017-01-20T16:26:27Z (SRR5621803 - nanopore4)
        # >85d5c1f1-fbc7-4dbf-bf17-215992eb7a08_Basecall_Barcoding_2d MinION2_MinION2_PoreCamp_FAA81710_GrpB2_1140_1_ch103_file54_strand /mnt/data/bioinformatics/Projects/MINion/Aureus_2998-174/Data/minion_run2(demultiplexed)/pass//MinION2_MinION2_PoreCamp_FAA81710_GrpB2_1140_1_ch103_file54_strand.fast5 (ERR1424936 - nanopore3)
        # @2ba6978b-759c-47a5-a3b9-77f03ca3ba66_Basecall_2D_template UNC_MJ0127RV_20160907_FN028_MN17984_sequencing_run_hansen_pool_32604_ch468_read5754_strand (SRR5817721 - nanopore3)
        # >9f328492-d6be-43b1-b2c9-3bf524508841_Basecall_Alignment_template minion_20160730_FNFAD16610_MN17211_sequencing_run_rhodotorula_36861_ch234_read302508_strand pass/minion_20160730_FNFAD16610_MN17211_sequencing_run_rhodotorula_36861_ch234_read302508_strand.fast5 (SRR5821557 - nanopore3)
        # >9f328492-d6be-43b1-b2c9-3bf524508841_Basecall_2D_2d minion_20160730_FNFAD16610_MN17211_sequencing_run_rhodotorula_36861_ch234_read302508_strand pass/minion_20160730_FNFAD16610_MN17211_sequencing_run_rhodotorula_36861_ch234_read302508_strand.fast5 (SRR5821557 - nanopore3)
        # @dd6189de-1023-4092-b1e9-76b291c07723_Basecall_1D_template Kelvin_20170412_FNFAF12973_MN19810_sequencing_run_170412_fungus_scedo_29052_ch239_read26_strand (SRR5812844 - nanopore3)
        # @channel_181_b44bbc58-3753-46d6-882c-0021c0697b55_template pass/6/Athena_20170324_FNFAF13858_MN19255_sequencing_run_fc2_real1_0_53723_ch181_read15475_strand.fast5 (SRR6329415 - nanopore3)
        # @channel_95_2663511f-6459-4e8b-8201-2360d199b9d8_template (SRR6329415 - nanopore)
        # @08441923-cb1a-490c-89b4-d209e234eb30_Basecall_1D_template GPBE6_F39_20170913_FAH26527_MN19835_sequencing_run_R265_r1_FAH26527_96320_read_2345_ch_440_strand fast5/GPBE6_F39_20170913_FAH26527_MN19835_sequencing_run_R265_r1_FAH26527_96320_read_2345_ch_440_strand.fast5 (SRR6377102 - nanopore3_1)
        # @d1aa0fa2-0e49-4030-b9c3-2785d2efd8ed_Basecall_1D_template ifik_cm401_20170420_FNFAE22716_MN16142_sequencing_run_CFsamplesB2bis_15133_ch74_read12210_strand
        # @6bbe187c-50a2-457e-997e-6be564f5980a_Basecall_2D 9B93VG2_20170410_FNFAF20574_MN19395_sequencing_run_AML_001_run2_82880_ch166_read9174_strand (WWW000050 - nanopore3 with no poreRead specified)
        # @aba5dfd4-af02-46d1-9bce-3b62557aa8c1 runid=91c917caaf7b201766339e506ba26eddaf8c06d9 read=29 ch=350 start_time=2018-03-02T16:12:39Z barcode=barcode01(SRR8695851 - nanopore4)
        # @72ad9b11-af72-4a2f-b943-650c9d88962f protocol_group_id=NASA_WCR_PCR_BC_083019 ch=503 barcode=BC08 read=17599 start_time=2019-08-30T21:26:18Z flow_cell_id=FAK67070 runid=4976f978bd6496df7a0ff30873137235afba9834 sample_id=NASA_WCR_083019 (SRR10303645 - nanopore4)
        ############################################################

        elif ( self.deflineType == self.NANOPORE or
               ( self.deflineType is None and
                 ( self.nanopore4.match(self.deflineString) or
                   self.nanopore.match( self.deflineString ) or
                   self.nanopore2.match( self.deflineString ) or
                   self.nanopore3.match( self.deflineString ) or
                   self.nanopore3_1.match( self.deflineString ) or
                   self.nanopore5.match( self.deflineString ) ) ) ):

            m = None
            if ( self.deflineType and
                 self.nanoporeSelected ):
                m = self.nanoporeSelected.match ( self.deflineString )
                self.foundRE = self.nanoporeSelected

            if ( self.deflineType is None or
                 m is None ):
                if self.nanopore.match ( self.deflineString ) :
                    m = self.nanopore.match ( self.deflineString )
                    self.foundRE = self.nanopore
                elif self.nanopore2.match ( self.deflineString ) :
                    m = self.nanopore2.match ( self.deflineString )
                    self.foundRE = self.nanopore2
                elif self.nanopore3.match ( self.deflineString ) :
                    m = self.nanopore3.match ( self.deflineString )
                    self.foundRE = self.nanopore3
                elif self.nanopore3_1.match ( self.deflineString ) :
                    m = self.nanopore3_1.match ( self.deflineString )
                    self.foundRE = self.nanopore3_1
                elif self.nanopore5.match ( self.deflineString ) :
                    m = self.nanopore5.match ( self.deflineString )
                    self.foundRE = self.nanopore5
                else:
                    m = self.nanopore4.match ( self.deflineString )
                    self.foundRE = self.nanopore4

                if m and self.saveDeflineType:
                    self.nanoporeSelected = self.foundRE
                    self.platform = "NANOPORE"

            # Confirm regular expression succeeded

            if m is None :
                self.isValid = False
                return self.isValid

            # Assign values (note that readNum may actually be a file number which is not the same)

            poreMid = None
            if self.foundRE == self.nanopore4:
                self.name = m.group(1)
                r = self.getPoreReadNo.search(self.deflineString)
                if r:
                    self.readNo = r.group(1)
                c = self.getPoreChannel.search(self.deflineString)
                if c:
                    self.channel = c.group(1)
                b = self.getPoreBarcode.search(self.deflineString)
                if ( b and
                     b.group(1) != 'unclassified'):
                    self.spotGroup = b.group(1)
            elif self.foundRE == self.nanopore3:
                ( prefix, self.name, self.suffix, self.poreRead, discard, poreStart, self.channel, poreMid, self.readNo, poreEnd, endSep ) = m.groups()
                self.poreFile = poreStart + self.channel + poreMid + self.readNo + poreEnd
            elif self.foundRE == self.nanopore3_1:
                ( prefix, self.name, self.suffix, self.poreRead, discard, poreStart, self.readNo, poreMid, self.channel, poreEnd, endSep ) = m.groups()
                self.poreFile = poreStart + self.readNo + poreMid + self.channel + poreEnd
            elif self.foundRE == self.nanopore5:
                ( self.name, self.suffix, self.poreRead, endSep ) = m.groups()
            else:
                ( poreStart, self.channel, poreMid, self.readNo, poreEnd, self.poreRead, self.poreFile, endSep ) = m.groups()
                if self.readNo:
                    self.name = poreStart + self.channel + poreMid + self.readNo + poreEnd
                else:
                    self.name = poreStart + self.channel + poreEnd
                    r = self.getPoreReadNo2.search(self.deflineString)
                    if r:
                        self.readNo = r.group(1)

            # Set readNum to None if actually a file number

            if ( ( poreMid == "_file" ) or
                 ( not self.readNo ) ):
                self.readNo = 0

            # Process poreFile if present

            if self.poreFile:

                # Check for 'pass' or 'fail'

                if self.getPorePass.search(self.poreFile):
                    self.filterRead = 0
                elif self.getPoreFail.search(self.poreFile):
                    self.filterRead = 1

                # Check for barcode

                b = self.getPoreBarcode2.search(self.poreFile)
                if b:
                    (self.spotGroup,delimiter) = b.groups()
                    if self.barcodeStringPresent.search(self.spotGroup):
                        self.spotGroup = re.sub(r'barcode(\d+)$',r'BC\1',self.spotGroup,1)

                # Split poreFile on '/' or '\' if present

                poreFileChunks = re.split(r'[/\\]',self.poreFile)
                if len ( poreFileChunks) > 1:
                    self.poreFile = poreFileChunks.pop()

            if self.poreRead and self.appendPoreReadToName:
                self.name += self.suffix + self.poreRead

            # Check for missing poreRead (from R-based poRe fastq dump) and normalize read type

            if ( self.poreRead is None or
                 self.poreRead == '' ):
                if self.filename:
                    if self.pore2Dpresent.search(self.filename):
                        self.poreRead = "2D"
                    elif self.poreTemplatePresent.search(self.filename):
                        self.poreRead = "template"
                    elif self.poreComplementPresent.search(self.filename):
                        self.poreRead = "complement"
                    else:
                        self.poreRead = ""

                    if self.poreRead and self.appendPoreReadToName:
                        self.name += self.suffix + "_" + self.poreRead

                else:
                    self.poreRead = ""
            elif ( self.poreRead == "_twodirections" or
                   self.poreRead[1:] == "2D" or
                   self.poreRead[0:3] == "_2d" ):
                self.poreRead = "2D"
            elif ( self.poreRead[0:9] == "_template" or
                   self.poreRead == "-1D" or
                   self.poreRead == ".1T" ):
                self.poreRead = "template"
            else:
                self.poreRead = "complement"

            if ( not self.deflineType and
                 self.saveDeflineType ):
                self.deflineType = self.NANOPORE

        ############################################################
        # read_id and barcode
        #
        # @12-Dfasci_84178 read_id=12-Dfasci::G2J4TZQ02D3VUU barcode=AAAAAATT (WWW000023)
        # @32_L3_60077 read_id=24_PPC4::HDOFHVG03GOW52 barcode=AAAAAACT (WWW000023)
        # @PF01_76 read_id=P2034:00008:00038 barcode=CTATACACT (SRR2420289)
        ############################################################

        elif ( self.deflineType == self.READID_BARCODE or
               ( self.deflineType is None and
                 self.readIdBarcode.match( self.deflineString ) ) ) :

            # Capture 'qiimeName', prefix, read_id, and barcode

            m = self.readIdBarcode.match ( self.deflineString )

            # Confirm regular expression succeeded

            if m is None:
                self.isValid = False
                return self.isValid

            # Get match values

            (self.qiimeName,self.prefix,self.name,self.spotGroup,endSep) = m.groups()

            # Prepend self.prefix onto self.name if it exists and not at start of self.qiimeName

            if ( self.prefix and
                 not re.search( re.escape(self.prefix[0:len(self.prefix)-2]), self.qiimeName ) ):
                self.name = self.prefix + self.name

            # Retain defline type if desired

            if ( not self.deflineType and
                 self.saveDeflineType ):
                self.deflineType = self.READID_BARCODE
                self.platform = "UNDEFINED"

        ############################################################
        # Capillary/Sanger fastq with template & dir for input to newbler
        # Template is used for read/spot name; first name is discarded
        #
        # @Msex-P09-F_A01 template=Msex-P09-A01 dir=fwd library=BAC_end (SRR2762665)
        # @Msex-P09-R_A01 template=Msex-P09-A01 dir=rev library=BAC_end
        # >bac-190o01.f template=190o01 dir=f library=BACends (from https://contig.wordpress.com/2011/01/21/newbler-input-ii-sequencing-reads-from-other-platforms)
        # >bac-190o01.r template=190o01 dir=r library=BACends
        # >originalreadname_1 template=originalreadname dir=F library=somename (from https://contig.wordpress.com/2010/06/10/running-newbler-de-novo-assembly/
        # >originalreadname_2 template=originalreadname dir=R library=somename
        # >DJS045A03F template=DJS054A03 dir=F library=DJS045 trim=12-543 (from http://454.com/downloads/my454/documentation/gs-flx-plus/454SeqSys_SWManual-v2.6_PartC_May2011.pdf)
        ############################################################

        elif ( self.deflineType == self.SANGER_NEWBLER or
               ( self.deflineType is None and
                 self.sangerNewbler.match ( self.deflineString ) ) ) :

            # Capture 'qiimeName', prefix, read_id, and barcode

            m = self.sangerNewbler.match( self.deflineString )

            # Confirm regular expression succeeded

            if m is None :
                self.isValid = False
                return self.isValid

            # Get match values (template becomes the name)

            (localName,self.name,self.dir,endSep) = m.groups()

            # Set readNum (expected to be char)

            if self.dir[0] in ('f', 'F'):
                self.readNum = '1'
            elif self.dir[0] in ('r', 'R'):
                self.readNum = '2'
            else:
                self.statusWriter.outputErrorAndExit( "Unexpected sanger read dir value ... {}".format(self.deflineString) )

            # Retain defline type if desired

            if ( not self.deflineType and
                 self.saveDeflineType ):
                self.deflineType = self.SANGER_NEWBLER
                self.platform = "CAPILLARY"

        ############################################################
        # ab solid without tag or suffix
        #
        # @1_8_179 (ERR2134219)
        ############################################################

        elif ( self.deflineType == self.ABSOLID2 or
               ( self.deflineType is None and
                 self.abSolid2.match(self.deflineString) ) ):

            m = self.abSolid2.match(self.deflineString)

            # Confirm regular expression succeeded

            if m is None :
                self.isValid = False
                return self.isValid

            # Get match values

            (self.panel, self.x, self.y, endSep) = m.groups()

            # Set defline name

            if self.abiTitle:
                self.name = f'{self.abiTitle}_{self.panel}_{self.x}_{self.y}'
            else:
                self.name = f'{self.panel}_{self.x}_{self.y}'

            # Retain defline type if desired

            if ( not self.deflineType and
                 self.saveDeflineType ):
                self.deflineType = self.ABSOLID2
                self.platform = "ABSOLID"

        ############################################################
        # Spot numbering or read numbering offset by '-', '/', or '_'
        # @HWI-ST765_r518_[12] (SRR2169708)
        # >Read_1 (SRR5350446)
        # @G1907151113-00-2_YSS-0715-F (SRR10053979)
        # @S3_332 (ERR1288564)
        ############################################################

        elif ( self.deflineType == self.ALTNUM or
               ( self.deflineType is None and
                 ( self.altNumOnly.match ( self.deflineString ) or
                   self.altNumOnly2.match ( self.deflineString ) ) ) ):

            m = None
            if ( self.deflineType and
                 self.altNumSelected ):
                m = self.altNumSelected.match(self.deflineString)

            # Not allowing switch back and forth between two altnum varieties

            if ( self.deflineType is None or
                 m is None ):
                if self.altNumOnly.match( self.deflineString ):
                    m = self.altNumOnly.match( self.deflineString )
                    self.foundRE = self.altNumOnly
                elif self.altNumOnly2.match( self.deflineString ):
                    m = self.altNumOnly2.match( self.deflineString )
                    self.foundRE = self.altNumOnly2

                if m and self.saveDeflineType:
                    self.altNumSelected = self.foundRE
                    self.platform = "UNDEFINED"

            # Confirm regular expression succeeded

            if m is None :
                self.isValid = False
                return self.isValid

            # Extract components

            (self.name, sep, self.altNum, endSep) = m.groups()

            # Convert F or R read identifier to 1 or 2.

            if self.altNum == 'F':
                self.altNum = '1'
            elif self.altNum == 'R':
                self.altNum = '2'

            # If sequential, add 'sep' and 'self.altNum' to name

            if ( self.retainAltNum or
                 int ( self.altNum ) > 4 ):
                self.name = self.name + sep + self.altNum
                self.retainAltNum = True

            # Retain type if desired

            if ( not self.deflineType and
                 self.saveDeflineType ):
                self.deflineType = self.ALTNUM
                self.platform = "UNDEFINED"

        ############################################################
        # Default just use non-whitespace characters after > or @
        #
        # @MB03~zSEQ034LingCncrtPool1~0000282 (SRR869399)
        # @Lab1.3.ab1 1249 10 1068 (WWW000012)
        # @CH_BAC1_C05.ab1_extraction_2 (WWW000011)
        # @G15-D_3_1_903_603_0.81 (XXX006565)
        # >No_name^M (WWW000021)
        # @SN7001204_0288_BH97LHADXX_R_SRi_L5503_L5508:11106:1433:74165 (WWW000035)
        # @HWI-ST170:292:8:1101:1239-2176 (WWW000036)
        # @Read_1-Barcode=BC001-PIPELINE=V41 length=6487 (WWW000038)
        # @contig_2_to_3_R_inner_clone_1_9589_A11_BJ-674528_048.ab1 (SRR2064214)
        # @sim_CFSAN001140-756880/1 (SRR3020730)
        ############################################################

        elif ( self.deflineType == self.UNDEFINED or
               ( self.deflineType is None and
                 self.generic.match ( self.deflineString ) ) ):

            # Capture name after '@' sign

            m = self.generic.match( self.deflineString )

            # Confirm regular expression succeeded

            if m is None :
                self.isValid = False
                return self.isValid

            # Get match value

            self.name = m.group(1)

            # Retain defline type if desired

            if ( not self.deflineType and
                 self.saveDeflineType ):
                self.deflineType = self.UNDEFINED
                self.platform = "UNDEFINED"

        else:
            self.isValid = False

        if self.isValid:

            # Ensure name contains at least one alphanumeric

            if ( self.name and
                 not self.alphanumericPresent.search(self.name) ):
                self.isValid = False
                self.deflineType = None

            # Make desired adjustments

            else:
                if self.ignLeadCharsNum:
                    self.ignoredLeadChars = self.name[0:self.ignLeadCharsNum]
                    self.name = self.name[self.ignLeadCharsNum:]
                if self.name and self.ignTailCharsNum:
                    self.ignoredTailChars = self.name[len(self.name)-self.ignTailCharsNum]
                    self.name = self.name[0:len(self.name)-self.ignTailCharsNum]
                if self.addFileToName and self.filenameTrunc:
                    self.name = self.filenameTrunc + '_' + self.name

        return self.isValid

        # SRA fastq output with read numbers
        # SRA fastq output without read numbers

    ############################################################
    # Count extra numbers in Illumina prefix (at most 2)
    ############################################################

    def countExtraNumbersInIllumina (self, sep):

        # Determine how many numbers at the end of self.prefix
        # separated by colons

        prefixChunks = re.split(sep,self.prefix)
        prefixChunksLen = len(prefixChunks)
        numCount = 0

        # Determine if prefix ends in one or two digits
        # delineated with 'sep'

        if prefixChunks[prefixChunksLen-1].isdigit():
            numCount += 1
        if ( prefixChunksLen > 1 and
             prefixChunks[prefixChunksLen-2].isdigit() ):
            numCount+=1

        # Determine how many numbers to discard (capping at 2 based on
        # what I have seen in the data)

        discardCount = 0
        if numCount > 0:
            numberChunks = re.split('\.',self.y)
            if int(numberChunks[0]) < 4:
                discardCount += 1
                numberChunks = re.split('\.',self.x)
                if ( numCount == 2 and
                     int(numberChunks[0]) < 4 ):
                    discardCount += 1

        return discardCount

############################################################
# Seq Class
############################################################

class Seq:
    """ Parses/validates sequence string """

    leadLowerBases = re.compile("(^[-actgnu.?]*)")
    tailLowerBases = re.compile("([-actgnu.?]*)$")

    def __init__(self, seqString, getClips=False):
        self.seqOrig = None
        self.seq = None
        self.length = 0
        self.isValid = False
        self.isColorSpace = False
        self.isBaseSpace = False
        self.clipLeft = 0
        self.clipRight = 0
        self.csKey = None
        self.transBase1 = str.maketrans('', '', "-ACTGNWSBVDHKMRYUX.?")
        self.transBase2 = str.maketrans('', '', "ACTGNWSBVDHKMRYUX")
        self.transColor = str.maketrans('', '', "0123.")
        if seqString:
            self.parseSeq(seqString,getClips)

    ############################################################
    # Determine number of bases to clip from left of sequence
    ############################################################

    @classmethod
    def getClipLeft (cls, seqOrig, seq):
        if seqOrig[0] == seq[0]:
            return 0
        else:
            m = Seq.leadLowerBases.search(seqOrig)
            clipString = m.group(1)
            return len(clipString)

    ############################################################
    # Determine number of bases to clip from right of sequence
    ############################################################

    @classmethod
    def getClipRight (cls, seqOrig, seq):
        seqLen = len(seq)
        if seqOrig[seqLen-1] == seq[seqLen-1]:
            return 0
        else:
            m = Seq.tailLowerBases.search(seqOrig)
            clipString = m.group(1)
            return len(clipString)

    ############################################################
    # Determine if provided seqString matches only sequence characters
    ############################################################

    def parseSeq ( self, seqString, getClips=False ):
        self.seqOrig = seqString.strip()
        self.seq = self.seqOrig.upper()
        self.length = len(self.seq)
        self.isValid = False
        self.clipLeft = 0
        self.clipRight = 0
        self.csKey = None

        if self.length == 0:
            pass

        elif self.isBaseSpace:
            empty = self.seq.translate(self.transBase1)
            if not empty:
                self.isValid = True
                if getClips:
                    self.clipLeft = Seq.getClipLeft(self.seqOrig,self.seq)
                    if self.clipLeft != self.length:
                        self.clipRight = Seq.getClipRight(self.seqOrig,self.seq)

        elif self.isColorSpace:
            empty2 = self.seq[1:].translate(self.transColor)
            if ( not empty2 and
                 self.seq[0] in "ACTG" ):
                self.isValid = True
                self.csKey = self.seq[0]
                self.seq = self.seq[1:]
                self.length -= 1

        else:
            empty = self.seq.translate(self.transBase2)
            if not empty:
                self.isValid = True
                self.isBaseSpace = True
                if getClips:
                    self.clipLeft = Seq.getClipLeft(self.seqOrig,self.seq)
                    if self.clipLeft != self.length:
                        self.clipRight = Seq.getClipRight(self.seqOrig,self.seq)

            else:
                empty2 = self.seq[1:].translate(self.transColor)
                if ( not empty2 and
                     self.seq[0] in "ACTG" ):
                    self.isValid = True
                    self.isColorSpace = True
                    self.csKey = self.seq[0]
                    self.seq = self.seq[1:]
                    self.length -= 1

                else:

                    # Check for non-colorspace seq with dots
                    # (2nd check here to properly handle colorspace seq consisting of all dots)

                    empty3 = self.seq.translate(self.transBase1)
                    if not empty3:
                        self.isValid = True
                        self.isBaseSpace = True
                        if getClips:
                            self.clipLeft = Seq.getClipLeft(self.seqOrig,self.seq)
                            if self.clipLeft != self.length:
                                self.clipRight = Seq.getClipRight(self.seqOrig,self.seq)

        return self.isValid

############################################################
# Qual Class
############################################################

class Qual:
    """ Parses/validates/characterizes quality string """

    def __init__(self, qualString, seqLen):
        self.qual = None
        self.length = 0
        self.isValid = False
        self.isNumQual = False
        self.isAscQual = False
        self.offset = 33
        self.minQual = 1000
        self.maxQual = 0
        self.minOrd = 1000
        self.maxOrd = 0
        self.seqLen = 0
        if seqLen:
            self.seqLen = seqLen
        if qualString:
            self.parseQual( qualString, seqLen )

    ############################################################
    @staticmethod
    def fixBadFlatNumQual( qualString ):
        qualFixed = ''
        for qualStr in qualString.strip().split():
            qLen = len ( qualStr )
            if qLen == 4:
                qualStr = qualStr[0:2] + ' ' + qualStr[2:]
            elif qLen == 3:
                q1 = qualStr[0:2]
                q2 = qualStr[0]
                q3 = qualStr[1:]
                q4 = qualStr[2:]
                if ( q2 == '0' or
                     ( int(q1) > 45 > int(q3) ) ):
                    qualStr = q2 + ' ' + q3
                else:
                    qualStr = q1 + ' ' + q4
            elif ( int(qualStr) > 45 or
                   qualStr[0] == '0' ):
                qualStr = qualStr[0] + ' ' + qualStr[1:]
            qualFixed += qualStr + ' '
        return qualFixed.strip()

    ############################################################
    # Determine if provided qualString matches only qual characters
    # in range or if numerical quality
    ############################################################

    def parseQual( self, qualString, seqLen ):
        self.qual = qualString.strip()
        self.length = 0
        self.isValid = False
        self.minQual = 1000
        self.maxQual = 0
        self.minOrd = 1000
        self.maxOrd = 0

        # isValid stays False if zero length

        if not self.qual:
            return self.isValid

        # Check for numerical quality

        singleIntQual = False
        if ( ( not seqLen or seqLen == 1 ) and
             self.isInt(self.qual) and
             int(self.qual) <= 100 ):
            singleIntQual = True

        if ( self.isNumQual or
             ( not self.isAscQual and
               ( " " in self.qual or
                 singleIntQual ) ) ):
            for qualSubString in self.qual.split():
                if self.isInt(qualSubString):
                    if int(qualSubString) < self.minQual:
                        self.minQual = int(qualSubString)
                    if int(qualSubString) > self.maxQual:
                        self.maxQual = int(qualSubString)
                    self.length += 1
                    if self.maxQual > 100:
                        sys.exit( "Numerical quality is too high (perhaps try '--badFlatNumQual')  ... {}".format(self.qual) )
                else:
                    self.length = 0
                    break

            if self.length != 0:
                self.isValid = True
                if not singleIntQual: # Not concluding numerical qual based on single integer quality
                    self.isNumQual = True

        # Check for ascii quality.
        # Defaulting to 33 offset for now
        # Subsequent processing will be used to better establish the offset
        # A space in ascii quality invalids the quality

        else:
            self.minOrd, self.maxOrd, spacePresent = cParseQual(self.qual)

            if self.maxOrd <= 126 and self.minOrd >= -40 and not spacePresent: # no special characters present
                self.length = len(self.qual)
                self.minQual = self.minOrd - self.offset
                self.maxQual = self.maxOrd - self.offset
                self.isValid = True
                self.isAscQual = True

        return self.isValid

    ############################################################
    # Check for integer
    ############################################################

    @staticmethod
    def isInt ( qStr ):
        try:
            _ = int(qStr)
            return True
        except:
            return False

############################################################
# FastqReader Class
############################################################

class FastqReader:
    """ Retains information read/parsed from fastq file """

    def __init__(self,filename,handle):
        self.filename = filename
        self.handle = handle
        self.qualHandle = None
        self.deflineCharSeq = "@"
        self.deflineCharQual = "+"
        self.deflineStringSeq = ''
        self.deflineStringQual = ''
        self.defline = Defline( None )
        self.defline.filename = self.filename
        self.deflineQual = None
        self.deflineCheck = None
        self.seqOrig = ''
        self.seq = ''
        self.qualOrig = ''
        self.qual = ''
        self.length = 0
        self.lengthQual = 0
        self.spotCount = 0
        self.lineCount = 0
        self.lineCountQual = 0
        self.prevLineCount = 0
        self.eof = False
        self.csKey = ''
        self.isColorSpace = False
        self.isNumQual = False
        self.isMultiLine = False
        self.isFasta = False
        self.offset = 33
        self.defaultQual = "?" # 30 using 33 offset
        self.isSplitSeqQual = False
        self.seqQualDeflineMismatch = False
        self.setClips = False
        self.clipLeft = 0
        self.clipRight = 0
        self.extraQualForCSkey = False
        self.checkSeqQual = True
        self.seqParser = Seq( '' )
        self.qualParser = Qual( '', 0 )
        self.seqLineCount = None
        self.qualLineCount = None
        self.savedDeflineString = ''
        self.savedDeflineStringSeq = ''
        self.savedDeflineStringQual = ''
        self.writingToArchive = False
        self.outputStatus = False
        self.statusWriter = None
        self.ignoreNames = False
        self.ignAndDiscardNames = False
        self.concatPairFiles = False
        self.firstSpotNames = {}
        self.firstSpotsCount = 0
        self.prevFilePos = None
        self.concatOne = None
        self.concatTwo = None
        self.concatPos = None
        self.badCtrlChars = re.compile(r'[\x00-\x08\x0B\x0C\x0E-\x1F\x7F-\xFF]')
        self.badCtrlCharsPlusCR = re.compile(r'[\x00-\x08\x0B\x0C\x0D\x0E-\x1F\x7F-\xFF]+') # Found CR in middle of sequences
        self.removeBadChars = False
        self.removeSeqSpaces = False
        self.discardCount = 0
        self.maxDiscardCount = 100
        self.maxSearchCount = 100
        self.maxDeflineLen = 1000
        self.deflineDiscardDueToLength = 0
        self.maxSeqLineCountAllowed = 20000
        self.badFlatNumQual = False
        self.line = ''
        self.prevLine = ''
        self.convertEmptyDeflines = False
        self.emptyDeflineCount = 0
        self.filenameTrunc = ''
        self.prevEmptyLine = 0
        self.prevEmptyLineQual = 0
        self.keepEmptySeqs = False
        self.emptySeq = False
        self.spotNumAtNameStart = False
        self.seqOnlyFile = False
        self.genbankFile = False
        self.fastqWithGtrThan = False
        self.fastqWithoutAtSign = False
        self.removeLineQuotes = False
        self.removeLastChar = False
        self.allowEarlyFileEnd = False
        self.abiLastPrimerBase = None
        self.badLineCounts = {}
        self.badLineFirstOccurs = {}
        self.corruptionStart = None
        self.corruptionLine = None
        self.truncateLine = None
        self.transDashUridineX = str.maketrans('-uUX?', 'NtTNN')
#       self.transDashUridineX = str.maketrans('-uUXWSBVDHKMRY', 'NtTNNNNNNNNNNN')
        handle.seek(0)
        self.headerLineCount = self.processHeader(handle,self.defline,self.seqOnlyFile,self.removeLastChar)
        self.deflineCheck = copy.copy(self.defline) # Must come after processHeader (for capturing ABI title if present)

    ############################################################
    # Initialize fastq reader to reflect current understanding
    # of the data. Some sw values may not be set yet
    # In retrospect, should just retain 'status' object.
    ############################################################

    def setStatus (self, status):
        self.isNumQual = status.isNumQual
        self.isColorSpace = status.isColorSpace
        self.extraQualForCSkey = status.extraQualForCSkey
        self.checkSeqQual = status.checkSeqQual
        self.setClips = status.setClips
        self.statusWriter = status.statusWriter
        self.ignoreNames = status.ignoreNames
        self.ignAndDiscardNames = status.ignAndDiscardNames
        self.concatPairFiles = status.concatPairFiles
        self.removeBadChars = status.removeBadChars
        self.removeSeqSpaces = status.removeSeqSpaces
        self.maxDiscardCount = status.maxErrorCount
        self.maxSearchCount = status.maxSearchCount
        self.maxDeflineLen = status.maxDeflineLen
        self.maxSeqLineCountAllowed = status.maxSeqLineCountAllowed
        self.badFlatNumQual = status.badFlatNumQual
        self.convertEmptyDeflines = status.convertEmptyDeflines
        self.fastqWithGtrThan = status.fastqWithGtrThan
        self.fastqWithoutAtSign = status.fastqWithoutAtSign
        self.removeLineQuotes = status.removeLineQuotes
        self.removeLastChar = status.removeLastChar
        self.abiLastPrimerBase = status.abiLastPrimerBase
        if self.convertEmptyDeflines:
            self.filenameTrunc = self.getFilenameTrunc(self.filename)
        self.keepEmptySeqs = status.keepEmptySeqs
        self.spotNumAtNameStart = status.spotNumAtNameStart # Only used for single line fastq in columns
        self.allowEarlyFileEnd = status.allowEarlyFileEnd
        self.seqOnlyFile = status.seqOnlyFile
        self.genbankFile = status.genbankFile
        if status.isMultiLine: # Can only turn on multiline - assuming if some files are multiline, then all could be
            self.isMultiLine = True
        self.offset = status.offset
        if self.offset == 64:
            self.defaultQual = "^"
        self.writingToArchive = status.writingToArchive
        self.defline.setStatus ( status )
        self.deflineCheck.setStatus ( status )
        if self.deflineQual:
            self.deflineQual.setStatus ( status )
        self.restart()

    ############################################################
    # Read line from fastq handle
    ############################################################

    def readFastqLine (self,handle):
        safe = False
        line = ''
        attemptCount = 0
        while True:
            attemptCount += 1
            try:
                line = handle.readline()
            except IOError as err:
                self.statusWriter.outputWarning("I/O error at line {} in file {}: {}"
                                                .format(self.lineCount+1,self.filename,err))

            except UnicodeDecodeError:
                if not self.allowEarlyFileEnd:
                    self.statusWriter.outputErrorAndExit( "UnicodeDecodeError occurred. Likely UTF-16LE file {} ended early at line {}. Use '--allowEarlyFileEnd' to allow load to proceed."
                                                          .format(self.filename,self.lineCount) )
                elif self.outputStatus:
                    self.statusWriter.outputWarning("UnicodeDecodeError at line {} in file {} - treating as EOF (typically truncated utf-16-le)"
                                                    .format(self.lineCount+1,self.filename))
            except ValueError:
                self.addToBadLines(line, self.lineCount)
                line = '\n'
                self.discardCount += 1
                if ( not self.maxDiscardCount or
                     self.discardCount < self.maxDiscardCount ):
                    self.statusWriter.outputWarning("ValueError at line {} in file {} - treating as empty line\n{}"
                                                    .format(self.lineCount+1,self.filename,sys.exc_info()))
                elif ( self.maxDiscardCount and
                       self.discardCount >= self.maxDiscardCount ):
                    self.processExcessiveDiscards()
            except EOFError:
                if not self.allowEarlyFileEnd:
                    self.statusWriter.outputErrorAndExit( "EOFError occurred. Likely gzipped file {} ended early at line {}. Use '--allowEarlyFileEnd' to allow load to proceed."
                                                          .format(self.filename,self.lineCount) )
                self.statusWriter.outputWarning("Failed to read at line {} in file {} - treating as EOF (typically truncated .gz)"
                                                .format(self.lineCount,self.filename) )
            except:
                if ( re.search ( "zlib.error", str(sys.exc_info()[0] ) ) and
                     not self.allowEarlyFileEnd ):
                    self.statusWriter.outputErrorAndExit( "zlib.error occurred. Likely gzipped file {} ended early at line {}. Use '--allowEarlyFileEnd' to allow load to proceed."
                                                          .format(self.filename,self.lineCount) )
                elif re.search ( "zlib.error", str (sys.exc_info()[0] ) ):
                    self.statusWriter.outputWarning("Failed to read at line {} in file {} - treating as EOF (typically truncated .gz)"
                                                    .format(self.lineCount,self.filename) )
                else:
                    self.statusWriter.outputErrorAndExit("Unexpected error while reading at line {} in file {}: {}"
                                                         .format(self.lineCount+1,self.filename,sys.exc_info()[0]) )

            if ( self.removeLineQuotes and
                 line.startswith('"') and
                 line[0:-1].endswith('"') ):
                line = line[1:-2] + line[-1]

            if line == '':
                break
            elif ( self.qualHandle and
                   handle == self.qualHandle):
                self.lineCountQual += 1
            else:
                self.lineCount += 1
                if ( self.truncateLine and
                     self.lineCount >= self.truncateLine ):
                    line = ''
                    break

            if self.removeBadChars:
                line = re.sub(self.badCtrlCharsPlusCR,'',line)
                if ord(line[0]) == 65279: # Check for and remove Byte Order Mark or BOM character
                    line = line[1:]

            line = line.lstrip(' \t') # Retention of linefeed is desired for whitespace only lines

            # Very limited solution for empty deflines - requires otherwise perfect 4-line fastq to work

            if self.convertEmptyDeflines:
                line_strip = line.strip()
                if ( (line_strip in ('', '0', '@', '@_', '@+') and self.lineCount % 4 == 1) or
                     (line_strip in ('>', '>_') and self.lineCount % 2 == 1) ):
                    self.emptyDeflineCount += 1
                    if line_strip in ('>', '>_'):
                        line = f'>{self.filenameTrunc}_{self.emptyDeflineCount}'
                    else:
                        line = f'@{self.filenameTrunc}_{self.emptyDeflineCount}'

            if ( line.strip() == '' and
                 ( not self.keepEmptySeqs ) ):
                if ( self.qualHandle and
                     handle == self.qualHandle):
                    self.prevEmptyLineQual = self.lineCountQual
                else:
                    self.prevEmptyLine = self.lineCount

            elif ( not self.removeBadChars and
                   self.badCtrlChars.search ( line ) ):
                self.addToBadLines(line, self.lineCount)
                self.discardCount += 1
                if ( not self.maxDiscardCount or
                     self.discardCount < self.maxDiscardCount ):
                    self.statusWriter.outputWarning("Discarding line {} containing unexpected characters in file {}"
                                                    .format(self.lineCount,self.filename) )
                elif ( self.maxDiscardCount and
                       self.discardCount >= self.maxDiscardCount ):
                    self.processExcessiveDiscards()

            else:
                #self.statusWriter.outputInfo("Safe line ... {}\t{}".format(self.lineCount,line) )
                safe = True
                break

        if ( self.removeLastChar and
             line.strip() != ''):
            line = line.strip()
            line = line[0:-1] + '\n'

        self.prevLine = self.line
        self.line = line

#        self.statusWriter.outputInfo("file {} line {}:\t{}".format(self.filename,self.lineCount,line) )

        if safe:
#            sys.stderr.write("{}".format(line))
            return line
        return ''

    ############################################################
    # Read from fastq handle
    ############################################################

    def read (self):

        while True:

            self.resetFastqValues()

            if self.concatPairFiles:
                self.prevFilePos = self.handle.tell()
                self.prevLineCount = self.lineCount
                if self.savedDeflineString:
                    self.prevFilePos -= ( len(self.savedDeflineString) + 1 )
                elif ( ( self.isFasta or self.isSplitSeqQual) and
                         self.savedDeflineStringSeq ):
                    self.prevFilePos -= ( len(self.savedDeflineStringSeq) + 1 )
                elif ( not ( self.isFasta or self.isSplitSeqQual) and
                       self.savedDeflineStringQual ):
                    self.prevFilePos -= ( len(self.savedDeflineStringQual) + 1 )

            # Process seq defline

            if self.seqOnlyFile:
                self.createSeqDefline()
            else:
                self.processSeqDefline()

            # Read seq

            self.length = self.readSeq()

            # Set qual values

            self.processQuality()

            # Break if not empty seq or eof; otherwise discard and try to read another read/spot

            if ( ( self.keepEmptySeqs and
                   self.emptySeq ) or
                 self.seq or
                 self.isEof() ):
                break

        # Check for excessive line discards

        if ( self.outputStatus and
             self.maxDiscardCount and
             self.discardCount > self.maxDiscardCount ):
            self.processExcessiveDiscards()

        elif ( ( self.keepEmptySeqs and
                 self.emptySeq ) or
               self.seq ):

            # Retain spot name if one of the first names and processing concatenated files

            if ( self.concatPairFiles and
                 ( self.maxSearchCount == 0 or self.firstSpotsCount < self.maxSearchCount ) and
                 self.firstSpotsCount < 1000000 ):
                self.firstSpotNames[self.defline.name] = 1
                self.firstSpotsCount += 1

            # Adjust for color space

            if ( self.isColorSpace or
                 self.seqParser.isColorSpace ):
                self.adjustForCS()

        # Increment spot count and check for eof
        # (or reached concat position when concatOne)

        self.updateSpotCountOrEof()

    ############################################################
    # Reset values before reading spot
    ############################################################

    def resetFastqValues (self):
        self.seqOrig = ''
        self.seq = ''
        self.qualOrig = ''
        self.qual = ''
        self.length = 0
        self.lengthQual = 0
        self.csKey = ''
        self.clipLeft = 0
        self.clipRight = 0

    ############################################################
    # Process quality related lines
    ############################################################

    def processQuality (self):

        # Process qual defline

        self.processQualDefline(self.handle)

        # Read qual

        if self.deflineStringQual:
            self.lengthQual = self.readQual( self.handle )

    ############################################################
    # Process seq defline
    ############################################################

    def processSeqDefline (self):

        if self.savedDeflineString:
            self.deflineStringSeq = self.savedDeflineString
            self.savedDeflineString = ''
        elif self.savedDeflineStringSeq:
            self.deflineStringSeq = self.savedDeflineStringSeq
            self.savedDeflineStringSeq = ''
        else:
            self.deflineStringSeq = self.readFastqLine(self.handle)
            if self.deflineStringSeq == '': # eof
                self.defline.reset()
                return

        self.deflineStringSeq = self.deflineStringSeq.strip()

        if self.genbankFile:
            m = re.match("LOCUS\s+(.+)",self.deflineStringSeq)
            if m:
                self.deflineStringSeq = '>' + m.group(1)
            else:
                self.statusWriter.outputErrorAndExit("Unable to find LOCUS at line {} in GenBank file {}"
                                                     .format(self.lineCount,self.filename) )
        elif ( self.fastqWithGtrThan and
               self.deflineStringSeq != '' and
               self.deflineStringSeq[0] == '>' ):
            self.deflineStringSeq = '@' + self.deflineStringSeq[1:]
        elif ( self.fastqWithoutAtSign and
               not self.deflineStringSeq[0] in "@+"):
            self.deflineStringSeq = '@' + self.deflineStringSeq

        self.deflineDiscardDueToLength = 0

        if ( self.deflineStringSeq and
             ( len ( self.deflineStringSeq ) > self.maxDeflineLen or
               self.deflineStringSeq[0] != self.deflineCharSeq or
               not self.defline.parseDeflineString( self.deflineStringSeq ) ) ):
            self.discardCount += 1
            self.addToBadLines(self.deflineStringSeq,self.lineCount)
            if self.outputStatus:
                if ( not self.maxDiscardCount or
                     self.discardCount < self.maxDiscardCount ):
                    if ( len ( self.deflineStringSeq ) > self.maxDeflineLen and
                         self.deflineStringSeq[0] == self.deflineCharSeq ):
                        self.statusWriter.outputWarning("Discarding line {} in processSeqDefline due to excessive length while looking for a valid defline in file {}"
                                                        .format(self.lineCount,self.filename) )
                        self.deflineDiscardDueToLength += 1
                    else:
                        self.statusWriter.outputWarning("Discarding line {} in processSeqDefline while looking for a valid defline in file {}"
                                                        .format(self.lineCount,self.filename) )

                elif ( self.maxDiscardCount and
                       self.discardCount > self.maxDiscardCount ):
                    self.processExcessiveDiscards()
            elif ( len ( self.deflineStringSeq ) > self.maxDeflineLen and
                   self.deflineStringSeq[0] == self.deflineCharSeq ):
                self.statusWriter.outputWarning("Discarding line {} in processSeqDefline due to excessive length while looking for a valid defline in file {}"
                                                .format(self.lineCount,self.filename) )
                self.deflineDiscardDueToLength += 1

            self.findValidDefline(True)

    ############################################################
    # Create defline with filename when seq-only file
    ############################################################

    def createSeqDefline (self):
        filenameTrunc = self.getFilenameTrunc(self.filename)
        self.deflineStringSeq = '>' + filenameTrunc
        self.defline.parseDeflineString( self.deflineStringSeq )

    ############################################################
    # Retain bad lines to look for patterns
    ############################################################
    def addToBadLines (self, line, lineCount):
        if ( self.maxDiscardCount and
             self.discardCount <= 1000 < lineCount and
             len(line) >= 5):
            if not ( line in self.badLineCounts ):
                self.badLineCounts[line] = 0
                self.badLineFirstOccurs[line] = lineCount
            self.badLineCounts[line] += 1

    ############################################################
    # Process retained bad lines
    ############################################################
    def processBadLines (self):
        sortedDict = OrderedDict( sorted(self.badLineCounts.items(), key = itemgetter(1), reverse = True) )
        for badLine in sortedDict:
            if sortedDict[badLine] >= 25:
                if ( not self.corruptionStart or
                     self.corruptionStart > self.badLineFirstOccurs[badLine] ):
                    self.corruptionStart = self.badLineFirstOccurs[badLine]
                    self.corruptionLine = badLine

    ############################################################
    # Exit due to excessive line discards
    ############################################################
    def processExcessiveDiscards (self,message='Excessive lines were discarded'):
        self.processBadLines()
        if self.corruptionStart:
            self.statusWriter.outputErrorAndExit(
                "{} due to repetitive corruption starting in {} at line {}"
                .format(message,self.filename, self.corruptionStart))
        else:
            if message == 'Excessive lines were discarded':
                badCount = self.discardCount
            else:
                badCount = self.maxSearchCount
            self.statusWriter.outputErrorAndExit(
                "{} for {} cases/lines in {} at line {}"
                .format(message,badCount, self.filename, self.lineCount))

    ############################################################
    # Get truncated file name without suffixes if possible
    ############################################################
    filenameMatch = re.compile(r"^(.*)[.](fastq|fq|txt|qual|fasta|csfasta|fa|fna|out|seq)(.gz|.bz2|)$")

    @staticmethod
    def getFilenameTrunc(filename):
        if FastqReader.filenameMatch.match(filename):
            m = FastqReader.filenameMatch.match(filename)
            (filenameTrunc, suffix1, suffix2 ) = m.groups()
            return filenameTrunc
        return filename

    ############################################################
    # Read seq and convert to uppercase
    ############################################################

    def readSeq (self):

        attemptCount = 0
        self.emptySeq = False

        if not self.getSeq(): # Populates self.seq and self.seqOrig and savedDeflines (multline only)
            return 0

        if self.checkSeqQual: # Skip this for early characterization of files
            while True:
                if ( ( self.keepEmptySeqs and
                       self.seq == '' ) or
                     self.seqParser.parseSeq( self.seq, self.setClips ) ):
                    if self.setClips:
                        self.clipLeft = self.seqParser.clipLeft
                        self.clipRight = self.seqParser.clipRight
                    break

                else:
                    if self.seq != '': # Ignore empty seqs rather than counting as an attempt
                        attemptCount += 1
                    self.seq = self.seqOrig.strip() # In case spaces were removed due to self.removeSeqSpaces

                    # Max attempt count (i.e. maxSearchCount) reached

                    if self.maxSearchCount != 0 and attemptCount > self.maxSearchCount:
                        self.processExcessiveDiscards('Unable to find a valid seq')

                    # Encountered seq is actually a defline

                    elif ( not self.isMultiLine and
                           ( self.seq[0] in (self.deflineCharQual, self.deflineCharSeq) or
                             ( self.fastqWithGtrThan and self.seq[0] == '>' ) or
                             self.fastqWithoutAtSign ) and
                           self.isDeflineString ( self.seq ) ):

                        # If seq defline, restart spot read

                        if ( self.seq[0] == self.deflineCharSeq or
                             ( self.fastqWithGtrThan and
                               self.seq[0] == '>' ) or
                             ( self.fastqWithoutAtSign and
                               self.seq[0] != '+' ) ):
                            self.discardCount += 1
                            self.addToBadLines(self.savedDeflineStringSeq,(self.lineCount-1))
                            if self.outputStatus:
                                if ( not self.maxDiscardCount or
                                     self.discardCount < self.maxDiscardCount ):
                                    self.statusWriter.outputWarning("Discarding line {} in readSeq due to absence of sequence between seq deflines in file {}"
                                                                    .format((self.lineCount-1),self.filename) )
                                elif ( self.maxDiscardCount and
                                       self.discardCount > self.maxDiscardCount ):
                                    self.processExcessiveDiscards()

                            self.savedDeflineStringSeq = self.seq
                            self.seq = ''
                            self.processSeqDefline()

                        # If qual defline, return with empty seq

                        elif not self.isSplitSeqQual:
                            # Assuming empty spot if previous line was empty - treating as not an error
                            if not (self.prevEmptyLine == (self.lineCount - 1) ):
                                self.discardCount += 1
                                self.addToBadLines(self.seq,self.lineCount)
                                if self.outputStatus:
                                    if ( not self.maxDiscardCount or
                                         self.discardCount < self.maxDiscardCount ):
                                        self.statusWriter.outputWarning("Discarding line {} in readSeq due to absence of sequence between seq and qual deflines in file {}"
                                                                        .format(self.lineCount,self.filename) )
                                    elif ( self.maxDiscardCount and
                                           self.discardCount > self.maxDiscardCount ):
                                        self.processExcessiveDiscards()

                            self.savedDeflineStringQual = self.seq
                            self.seq = ''
                            break

                        # Just skip this line (unlikely but just in case)

                        else:
                            self.seq = ''
                            break

                    # Not valid seq and not seq string is defline (non-multiline case)

                    else:
                        self.discardCount += 1
                        self.addToBadLines(self.seq,self.lineCount)
                        if self.outputStatus:
                            if ( not self.maxDiscardCount or
                                 self.discardCount < self.maxDiscardCount ):
                                self.statusWriter.outputWarning("Discarding line {} in readSeq while looking for valid sequence in file {}"
                                                                .format(self.lineCount,self.filename) )
                            elif ( self.maxDiscardCount and
                                   self.discardCount > self.maxDiscardCount ):
                                self.processExcessiveDiscards()

                        if self.savedDeflineStringSeq: # from multiline seq read that resulted in bad seq but valid seq defline
                            self.processSeqDefline()
                        else:
                            self.findValidDefline(True) # Start spot read over

                    if not self.getSeq(): # Populates self.seq and self.seqOrig and savedDeflines (multline only)
                        return 0

        self.seq = self.seq.upper()
        if ( self.keepEmptySeqs and
             self.seq == '' ):
            self.emptySeq = True
        elif self.writingToArchive: # If I do this too early it messes up other checks (i.e. -1 in quality is changed to N1)
            self.seq = self.seq.translate(self.transDashUridineX)

        return len ( self.seq )

    ############################################################
    # Get seq string
    ############################################################
    def getSeq (self):
        if self.isMultiLine:
            savedDeflineString = self.readMultiLineSeq()

            # Should be a qual defline with a '+' but checking for seq defline just in case

            if ( self.fastqWithGtrThan and
                 savedDeflineString and
                 savedDeflineString[0] == '>' ):
                savedDeflineString = '@' + savedDeflineString[1:]
            elif ( self.fastqWithoutAtSign and
                   savedDeflineString and
                   savedDeflineString[0] not in '+@' ):
                self.deflineStringSeq = '@' + self.deflineStringSeq

            if not savedDeflineString and not self.seqOrig:
                self.seq = ''
                return 0

            if savedDeflineString and savedDeflineString[0] == self.deflineCharSeq:
                self.savedDeflineStringSeq = savedDeflineString
            elif savedDeflineString and not self.isSplitSeqQual:
                self.savedDeflineStringQual = savedDeflineString
            self.seq = self.seqOrig

        else:
            self.seqOrig = self.readFastqLine(self.handle)
            if not self.seqOrig: # eof
                self.seq = ''
                return 0

            self.seq = self.seqOrig.strip()

        if self.removeSeqSpaces:
            self.seq = self.seq.replace(' ','')

        if ( self.abiLastPrimerBase and
             not self.isDeflineString ( self.seq ) ):
            self.seq = self.abiLastPrimerBase + self.seq
            self.seqOrig = self.abiLastPrimerBase + self.seqOrig

        return 1

    ############################################################
    # Process qual defline string
    ############################################################

    def processQualDefline (self,handle):
        if self.savedDeflineStringQual:
            self.deflineStringQual = self.savedDeflineStringQual
            self.savedDeflineStringQual = ''

        # Handle cases where quality defline and quality are absent (multiline seq)

        elif ( not self.isSplitSeqQual and
                not self.savedDeflineStringSeq == '' ):
            self.deflineStringQual = self.savedDeflineStringSeq

        else:
            self.deflineStringQual = self.readFastqLine(handle)

        attemptCount = 0
        while True:
            attemptCount += 1
            if not self.deflineStringQual: # eof
                if self.isSplitSeqQual:
                    self.deflineQual.reset()
                break
            else:
                self.deflineStringQual = self.deflineStringQual.strip()

                if ( self.deflineStringQual[0] == self.deflineCharQual and
                     self.isDeflineString( self.deflineStringQual ) ):
                    break
                elif ( not self.isSplitSeqQual and
                       ( self.deflineStringQual[0] == '@' or
                         ( self.fastqWithGtrThan and
                           self.deflineStringQual[0] == '>' ) or
                         ( self.fastqWithoutAtSign and
                           self.deflineStringQual[0] != '+' ) ) and
                       self.isDeflineString(self.deflineStringQual ) ):
                    self.savedDeflineStringSeq = self.deflineStringQual
                    self.deflineStringQual = ''
                    break
                elif self.maxSearchCount != 0 and attemptCount > self.maxSearchCount:
                    self.processExcessiveDiscards('Unable to find a valid qual defline')
                else:
                    self.discardCount += 1
                    self.addToBadLines(self.deflineStringQual,self.lineCount)
                    if self.outputStatus:
                        if ( not self.maxDiscardCount or
                             self.discardCount < self.maxDiscardCount ):
                            self.statusWriter.outputWarning("Discarding line {} in processQualDefline while looking for a valid qual defline in file {}"
                                                        .format(self.lineCount,self.filename) )
                        elif ( self.maxDiscardCount and
                               self.discardCount > self.maxDiscardCount ):
                            self.processExcessiveDiscards()

                self.deflineStringQual = self.readFastqLine(handle)

    ############################################################
    # Read qual (possibly from separate file)
    ############################################################

    def readQual ( self, handle ):

        if self.isMultiLine:
            savedDeflineString = self.readMultiLineQual( handle )
            if not savedDeflineString and not self.qualOrig:
                self.qual = ''
                return 0

            if self.isSplitSeqQual:
                self.savedDeflineStringQual = savedDeflineString
            elif not self.isSplitSeqQual and savedDeflineString:
                self.savedDeflineStringSeq = savedDeflineString
            if not self.qualOrig.rstrip(): # Went direct to seq defline - qual missing
                return 0
        else:
            self.qualOrig = self.readFastqLine(handle)
            if not self.qualOrig: # eof
                self.qual = ''
                return 0

        self.qual = self.qualOrig.rstrip() # leave leading space here; just strip trailing spaces
        if ( self.qual and
             self.badFlatNumQual and
             " " in self.qual):
            self.qual = self.qualParser.fixBadFlatNumQual(self.qual)

        if ( self.qual and
             self.isNumQual ):
            lengthQual = len ( self.qual.split() )
        else:
            lengthQual = len ( self.qual )

        if self.checkSeqQual: # Skipped for early characterization of files
            if ( lengthQual != self.length and
                 ( self.qual[0] in (self.deflineCharQual, self.deflineCharSeq) or
                   ( self.fastqWithGtrThan and self.qual[0] == '>' ) or
                   self.fastqWithoutAtSign ) and
                 self.isDeflineStringQual(self.qual) ):
                if ( self.isSplitSeqQual and
                     self.qual[0] == '>' ):
                    self.savedDeflineStringQual = self.qual
                elif self.qual[0] == '+': # This qual defline will not be used
                    self.discardCount += 1
                    self.addToBadLines(self.qual,self.lineCount)
                    if self.outputStatus:
                        if ( not self.maxDiscardCount or
                             self.discardCount < self.maxDiscardCount ):
                            self.statusWriter.outputWarning("Discarding line {} in readQual (unexpected qual defline) while looking for valid quality in file {}"
                                                            .format(self.lineCount,self.filename) )
                        elif ( self.maxDiscardCount and
                               self.discardCount > self.maxDiscardCount ):
                            self.processExcessiveDiscards()
                else:
                    self.savedDeflineStringSeq = self.qual
                self.qual = ''
            elif ( ( self.keepEmptySeqs and
                     not self.qual and
                     self.emptySeq ) or
                   ( self.qual and
                     self.qualParser.parseQual( self.qual, self.length) ) ):
                pass
            else:
                self.discardCount += 1
                self.addToBadLines(self.qual,self.lineCount)
                if self.outputStatus:
                    if ( not self.maxDiscardCount or
                         self.discardCount < self.maxDiscardCount ):
                        self.statusWriter.outputWarning("Discarding line {} in readQual while looking for valid quality in file {}"
                                                        .format(self.lineCount,self.filename) )
                    elif ( self.maxDiscardCount and
                           self.discardCount > self.maxDiscardCount ):
                        self.processExcessiveDiscards()
                self.qual = ''

        if ( self.qual and
             self.isNumQual ):
            lengthQual = len ( self.qual.split() )
        else:
            lengthQual = len ( self.qual )

        if ( self.qual and
             lengthQual != self.length and
             self.qual[0] == '"' and
             not self.qual[1:2] == '"' and
             self.qual[lengthQual-1] == '"' ):
            self.qual = self.qual[1:lengthQual-1]
            lengthQual -= 2

##        if ( ( not self.keepEmptySeqs ) and
##             self.length == 0 ): # Seq was not found but going through quality to get to next seq defline
##            if ( not self.maxDiscardCount or
##                     self.discardCount < self.maxDiscardCount ):
##                self.statusWriter.outputWarning("Line {} in file {} (and perhaps prior lines) will be discarded due to being preceded by an empty sequence"
##                                                .format(self.lineCount,self.filename) )

        return lengthQual

    ############################################################
    # Read seq occurring on multiple lines
    ############################################################

    def readMultiLineSeq (self):
        self.seqLineCount = 0
        self.seqOrig = ''
        seqList = []

        while True:
            fileStringOrig = self.readFastqLine(self.handle)
            if not fileStringOrig: # eof
                self.seqOrig = ''.join(seqList)
                return ''

            fileString = fileStringOrig.strip()

            if self.genbankFile:
                if fileString == "//":
                    self.seqOrig = ''.join(seqList)
                    return ''
                else:
                    m = re.match("^\s*\d+\s+(.+)",fileString)
                    if m:
                        fileString = m.group(1)
                    else:
                        fileString = ''

            if fileString and self.removeSeqSpaces:
                fileString = fileString.replace(' ','')

            # Check for next seq or qual defline

            if ( ( fileString and fileString[0] in '@+>' ) or
                   ( self.fastqWithoutAtSign and self.isDeflineString ( fileString ) ) ):
                fileString = fileStringOrig.strip() # In case spaces were removed
                self.seqOrig = ''.join(seqList)
                return fileString

            # Collect seq

            elif fileString and self.seqParser.parseSeq( fileString, self.setClips ):
                seqList.append(fileString)
                
                # If self.maxSeqLineCountAllowed lines between deflines

                self.seqLineCount += 1
                if self.seqLineCount > self.maxSeqLineCountAllowed:
                    self.statusWriter.outputErrorAndExit("Exceeded {} sequence lines between deflines near line {} in file {}"
                                                         .format(self.maxSeqLineCountAllowed,self.lineCount,self.filename) )

    ############################################################
    # Read multi-line qual (handle can be same or separate file)
    ############################################################

    def readMultiLineQual ( self, handle ):
        self.qualLineCount = 0
        self.qualOrig = ''

        # When validating found qual defline normally use self.deflineStringSeq

        deflineStringSeq = self.deflineStringSeq

        # For split seq and qual, use self.savedDeflineStringSeq for validation of defline

        if self.isSplitSeqQual:
            deflineStringSeq = self.savedDeflineStringSeq

        while True:
            fileStringOrig = self.readFastqLine(handle)
            fileString = fileStringOrig.rstrip()

            # At EOF/self.maxSeqLineCountAllowed lines and did not find next read defline

            if not fileStringOrig: # eof
                return ''

            elif self.qualLineCount > self.maxSeqLineCountAllowed:
                self.statusWriter.outputErrorAndExit("Exceeded {} quality lines between deflines near line {} in file {}\n"
                                                     .format(self.maxSeqLineCountAllowed,self.lineCount,self.filename))

            # Check for next seq defline. '+' and '>' are valid qualities so
            # have to depend on length check below for qual defline

            elif ( ( ( fileString and
                       fileString[0] == self.deflineCharSeq ) or
                     ( self.fastqWithGtrThan and
                       fileString[0] == '>' ) or
                     self.fastqWithoutAtSign ) and
                   len(fileString) <= self.maxDeflineLen and
                   self.isDeflineStringQual(fileString,deflineStringSeq) ):
                return fileString

            # Account for wrapped numerical quality

            if ( self.isNumQual or
                 ( " " in fileString and
                 re.search ( "^-?\d+\s+-?\d+\s+", fileString ) ) ):
                self.qualOrig += " "
                self.isNumQual = True

            # Collect qual

            self.qualLineCount += 1
            self.qualOrig += fileString

            if self.isNumQual: # Must be in the loop
                lengthQual = len ( self.qualOrig.split() ) - 1
            else:
                lengthQual = len ( self.qualOrig )

            if lengthQual >= self.length: # Also must be in the loop
                if not self.isNumQual: # Sometimes ascii qual longer than seq length contains a valid defline in corrupted files
                    self.qualOrig = self.qualOrig[0:self.length]
                return '' # No retained defline for next spot

    ############################################################
    # See if provided string represents a valid defline string when reading quality
    ############################################################
    def isDeflineStringQual ( self, fileString, deflineStringSeq=None ):

        if deflineStringSeq is None:
            deflineStringSeq = self.deflineStringSeq

        if self.fastqWithGtrThan:
            if fileString[0] == '>':
                fileString = '@' + fileString[1:]
            if deflineStringSeq[0] == '>':
                deflineStringSeq = '@' + deflineStringSeq[1:]
        if self.fastqWithoutAtSign:
            if fileString[0] not in '@+':
                fileString = '@' + fileString
            if deflineStringSeq[0] not in '@+':
                deflineStringSeq = '@' + deflineStringSeq

        if ( fileString == '+' and
             not self.isSplitSeqQual ):
            return True
        elif ( ( ( ( not self.defline.saveDeflineType or
                     self.defline.deflineType != Defline.UNDEFINED ) and # just match @|> for determined platform or not saving defline type
                   ( fileString[0:3] == deflineStringSeq[0:3] or
                     ( fileString[0] == deflineStringSeq[0] and
                       self.defline.deflineType in
                       ( Defline.NANOPORE, Defline.ABSOLID, Defline.READID_BARCODE, Defline.QIIME_GENERIC,
                         Defline.QIIME_454, Defline.QIIME_454_BC, Defline.QIIME_ILLUMINA_OLD,
                         Defline.QIIME_ILLUMINA_OLD_BC, Defline.QIIME_ILLUMINA_NEW, Defline.QIIME_ILLUMINA_NEW_BC,
                         Defline.QIIME_ILLUMINA_NEW_DBL, Defline.QIIME_ILLUMINA_NEW_DBL_BC ) ) ) ) or
                 ( self.defline.deflineType == Defline.UNDEFINED and # match @|> and next character for undefined platform
                   ( ( fileString[0:3] == deflineStringSeq[0:3] or
                       self.qualHandle and fileString[0:3] == self.deflineStringQual[0:3] ) or
                     ( re.search ( "^\d+", fileString[1:] ) and # match @|> and both start with digits for undefined platform
                       ( fileString[0] == '@' and
                         re.search ( "^\d+", deflineStringSeq[1:] ) ) or
                       ( self.qualHandle and
                         fileString[0] == '>' and
                         re.search ( "^\d+", self.deflineStringQual[1:] ) ) ) ) )  or
                 ( ( self.ignoreNames or self.ignAndDiscardNames) and
                   self.qualHandle and
                   fileString[0:3] == deflineStringSeq[0:3] ) ) and
               self.isDeflineString(fileString) ):
            return True
        else:
            return False

    ############################################################
    # See if provided string represents a valid defline string
    ############################################################
    def isDeflineString (self, fileString):
        if self.defline.saveDeflineType:
            # Optimization: save on name resolution
            deflineCheck = self.deflineCheck
            defline = self.defline

            deflineCheck.deflineType = defline.deflineType
            deflineCheck.illuminaNewSelected = defline.illuminaNewSelected
            deflineCheck.illuminaOldSelected = defline.illuminaOldSelected
            deflineCheck.illuminaOldBcRnSelected = defline.illuminaOldBcRnSelected
            deflineCheck.qiimeSelected = defline.qiimeSelected
            deflineCheck.pacbioSelected = defline.pacbioSelected
            deflineCheck.nanoporeSelected = defline.nanoporeSelected
            deflineCheck.altNumSelected = defline.altNumSelected
            deflineCheck.ionTorrentSelected = defline.ionTorrentSelected
            deflineCheck.bgiSelected = defline.bgiSelected
            deflineCheck.platform = defline.platform
        return self.deflineCheck.parseDeflineString(fileString)

    ############################################################
    # Find valid seq defline (ignore qual deflines)
    # Unlikely to succeed for split seq/qual
    ############################################################

    def findValidDefline (self,callProcessSeqDefline):
        count = 0
        while True:
            count += 1
            deflineStringSeqOrig = self.readFastqLine(self.handle)

            if not deflineStringSeqOrig: # eof
                self.deflineStringSeq = ''
                self.defline.reset()
                break
            else:
                self.deflineStringSeq = deflineStringSeqOrig.strip()

                if ( self.deflineStringSeq and
                     len(self.deflineStringSeq) <= self.maxDeflineLen and
                     ( self.deflineStringSeq[0] == self.deflineCharSeq or
                       ( self.fastqWithGtrThan and
                         self.deflineStringSeq[0] == '>' ) or
                       ( self.fastqWithoutAtSign and
                           self.deflineStringSeq[0] != '+' ) ) and
                     self.isDeflineString ( self.deflineStringSeq ) ):
                    self.savedDeflineStringSeq = self.deflineStringSeq
                    if callProcessSeqDefline:
                        self.processSeqDefline()
                    break

                elif self.maxSearchCount != 0 and count == self.maxSearchCount:
                    if self.deflineDiscardDueToLength > 0:
                        self.processExcessiveDiscards ('Unable to find a valid defline (maybe due to unexpected lengths)')
                    else:
                        self.processExcessiveDiscards ('Unable to find a valid defline')

                else:
                    self.discardCount += 1
                    self.addToBadLines(self.deflineStringSeq,self.lineCount)
                    if self.outputStatus:
                        if ( not self.maxDiscardCount or
                             self.discardCount < self.maxDiscardCount ):
                            if ( len ( self.deflineStringSeq ) > self.maxDeflineLen and
                                 ( self.deflineStringSeq[0] == self.deflineCharSeq or
                                   ( self.fastqWithGtrThan and
                                     self.deflineStringSeq[0] == '>' ) or
                                   ( self.fastqWithoutAtSign and
                                     self.deflineStringSeq[0] != '+' and
                                     self.isDeflineString ( self.deflineStringSeq ) ) ) ):
                                self.statusWriter.outputWarning("Discarding line {} in processSeqDefline due to excessive length while looking for a valid defline in file {}"
                                                                .format(self.lineCount,self.filename) )
                                self.deflineDiscardDueToLength += 1
                            else:
                                self.statusWriter.outputWarning("Discarding line {} in findValidDefline while looking for a valid defline in file {}"
                                                                .format(self.lineCount,self.filename) )
                        elif ( self.maxDiscardCount and
                               self.discardCount > self.maxDiscardCount ):
                            self.processExcessiveDiscards()

                    elif ( len ( self.deflineStringSeq ) > self.maxDeflineLen and
                           ( self.deflineStringSeq[0] == self.deflineCharSeq or
                             ( self.fastqWithGtrThan and
                               self.deflineStringSeq[0] == '>' ) or
                             ( self.fastqWithoutAtSign and
                               self.deflineStringSeq[0] != '+' and
                               self.isDeflineString ( self.deflineStringSeq ) ) ) ):
                        self.statusWriter.outputWarning("Discarding line {} in processSeqDefline due to excessive length while looking for a valid defline in file {}"
                                                        .format(self.lineCount,self.filename) )
                        self.deflineDiscardDueToLength += 1

    ############################################################
    # Update spotcount or eof based on quality and defline name
    ############################################################

    def updateSpotCountOrEof (self):
        if self.isEof():
            self.eof = True
        else:
            self.spotCount += 1

    ############################################################
    # EOF check (orig strings should still have linefeed)
    # Possible to have partial spot without quality at the end
    # of a file - thus the check for seqOrig, too
    ############################################################
    def isEof (self):
        return ( self.lineCount > 0 and
                 ( ( not self.qualOrig and
                     not self.seqOrig and
                     not self.savedDeflineStringSeq and
                     not self.savedDeflineStringQual ) or
                   ( self.concatOne and
                     self.prevFilePos >= self.concatPos ) ) )

    ############################################################
    # Adjust for CS
    ############################################################

    def adjustForCS (self):
        self.csKey = self.seq[0]
        self.seq = self.seq[1:]
        self.length -= 1
        space = self.qual.find(' ')
        if self.extraQualForCSkey:
            if space != -1:
                self.qual = self.qual[space+1:]
            else:
                self.qual = self.qual[1:]
            self.lengthQual -= 1
        if ( space != -1 and
             self.qual.find("-1") != -1 ):
            self.qual = self.qual.replace("-1","0")

    ############################################################
    # Start reading from front of fastq again
    ############################################################

    def restart (self):
        self.handle.seek(0)
        self.processHeader(self.handle,self.defline,self.seqOnlyFile,self.removeLastChar)
        if self.qualHandle:
            self.qualHandle.seek(0)
            self.processHeader (self.qualHandle,self.deflineQual,self.seqOnlyFile,self.removeLastChar)
            if ( self.defline.abiTitle and
                 not self.deflineQual.abiTitle ):
                self.defline.abiTitle = ''
            elif ( self.deflineQual.abiTitle
                   and not self.defline.abiTitle ):
                self.deflineQual.abiTitle = ''
        self.eof = False
        self.spotCount = 0
        self.lineCount = 0
        self.deflineStringSeq = ''
        self.deflineStringQual = ''
        self.seqQualDeflineMismatch = False
        self.savedDeflineString = ''
        self.savedDeflineStringSeq = ''
        self.savedDeflineStringQual = ''
        self.checkSeqQual = True
        self.firstSpotNames = {}
        self.firstSpotsCount = 0
        self.prevFilePos = None
        self.discardCount = 0
        self.emptyDeflineCount = 0
        self.line = ''
        self.badLineCounts = {}
        self.badLineFirstOccurs = {}
        self.corruptionStart = None
        self.corruptionLine = None

    ############################################################
    # Skip lines at start of file beginning with '#'
    ############################################################

    @staticmethod
    def processHeader (handle,defline,seqOnlyFile,removeLastChar):
        headerLineCount = 0
        while True:
            prevPos = handle.tell()
            fileString = handle.readline()
            try:
                fileString = fileString.decode("utf-8") # In case reading from bz2 file
            except:
                pass

            if removeLastChar and fileString.strip():
                fileString = fileString.strip()
                fileString = fileString[0:-1] + '\n'

            if seqOnlyFile:
                parsedSeq = Seq(fileString)
                if fileString.strip() and parsedSeq.isValid:
                    handle.seek(prevPos)
                    break
                else:
                    headerLineCount+=1
            else:
                if fileString.strip() and fileString[0] != "#":
                    handle.seek(prevPos)
                    break
                else:
                    headerLineCount+=1

                # Check for abi title

                if ( not defline.abiTitle and
                     re.match("# Title: ([!-~]+)",fileString) ):
                    m = re.match("# Title: ([!-~]+)",fileString)
                    defline.abiTitle = m.group(1)
                    if re.search ( "_$", defline.abiTitle ):
                        defline.abiTitle = defline.abiTitle[0:len(defline.abiTitle)-1]

        return headerLineCount

    ############################################################
    # Check if fastq file contains multi-line fastq
    ############################################################

    def isMultiLineFastq ( self ):
        self.seqLineCount = 0
        self.qualLineCount = 0

        self.restart()
        maxSeqLineCount = 0
        maxQualLineCount = 0
        while ( self.spotCount < 1001 and
                not ( ( self.lineCount > 0 and
                        not self.savedDeflineStringSeq and
                        ( not self.qualOrig or
                          not self.seqOrig ) ) or
                      self.seqLineCount > self.maxSeqLineCountAllowed or
                      self.qualLineCount > self.maxSeqLineCountAllowed ) ):
            self.read()
            if self.seqLineCount > maxSeqLineCount:
                maxSeqLineCount = self.seqLineCount
            if self.qualLineCount > maxQualLineCount:
                maxQualLineCount = self.qualLineCount

        self.restart()

        if ( 1 < maxSeqLineCount <= self.maxSeqLineCountAllowed or
             1 < maxQualLineCount <= self.maxSeqLineCountAllowed ):
            return True

        elif ( maxSeqLineCount > self.maxSeqLineCountAllowed or
               maxQualLineCount > self.maxSeqLineCountAllowed ):
            self.statusWriter.outputErrorAndExit( "Distance between fastq deflines exceeds {} lines in file ... {}"
                                                  .format(self.maxSeqLineCountAllowed,self.filename) )

        else:
            return False

    ############################################################
    # Find second read start in concat file
    ############################################################

    def findConcatSecReadStart (self):
        while ( ( ( self.maxSearchCount and self.spotCount <= self.maxSearchCount ) or
                  not self.defline.name in self.firstSpotNames ) and
                not self.eof ):
            self.read()
        if self.eof:
            return 0
        return self.prevFilePos

############################################################
# Split seq/qual Fastq Reader (i.e. separate seq and qual files)
############################################################

class SeqQualFastqReader (FastqReader):
    """ Split seq/qual variation on FastqReader """

    def __init__(self,filename, handle, qualFilename, qualHandle):
        FastqReader.__init__(self, filename, handle)
        qualHandle.seek(0)
        self.isSplitSeqQual = True
        self.isSeqQualCheck = False
        self.deflineCharSeq = ">"
        self.deflineCharQual = ">"
        self.qualFilename = qualFilename
        self.qualHandle = qualHandle
        self.deflineQual = Defline( '' )
        self.deflineQual.filename = self.qualFilename
        self.headerLineCountQual = self.processHeader(qualHandle,self.deflineQual,self.seqOnlyFile,self.removeLastChar)

    ############################################################
    # Process seq defline
    ############################################################

    def processSeqDefline (self):

        # Determine which string to used for defline

        if self.savedDeflineStringSeq:
            self.deflineStringSeq = self.savedDeflineStringSeq
            self.savedDeflineStringSeq = ''
        else:
            self.deflineStringSeq = self.readFastqLine(self.handle)

        if not self.deflineStringSeq:
            self.defline.reset()
        else:
            self.deflineStringSeq = self.deflineStringSeq.strip()

            # Parse defline and confirm validity

            self.defline.parseDeflineString ( self.deflineStringSeq )

            if not self.defline.isValid:
                attempts = 0
                while True:
                    attempts += 1
                    self.deflineStringSeq = self.readFastqLine(self.handle)
                    if not self.deflineStringSeq: # Seq file can end early
                        self.defline.reset()
                        return # End of truncated seq file, ignore rest of qual file?
                    else:
                        self.deflineStringSeq = self.deflineStringSeq.strip()
                        if (self.convertEmptyDeflines and
                            self.deflineStringSeq in ('>', '>_')):
                            self.emptyDeflineCount += 1
                            self.deflineStringSeq = f'>{self.filenameTrunc}_{self.emptyDeflineCount}'
                        self.defline.parseDeflineString( self.deflineStringSeq )
                        if self.defline.isValid:
                            break
                        self.addToBadLines(self.deflineStringSeq,self.lineCount)
                        if ( ( self.maxSearchCount != 0 and attempts == self.maxSearchCount ) or
                               self.ignoreNames or
                               self.ignAndDiscardNames ):
                            self.processExcessiveDiscards('Unable to find a valid seq defline')

            if not self.defline.isValid:
                if self.outputStatus:
                    self.statusWriter.outputErrorAndExit("Unable to parse defline in filename {} at line {}"
                                                         .format(self.filename,self.lineCount))

    ############################################################
    # Process qual defline
    ############################################################

    def processQualDefline ( self, handle ):
        if self.savedDeflineStringQual:
            self.deflineStringQual= self.savedDeflineStringQual
            self.savedDeflineStringQual = ''
        else:
            self.deflineStringQual = self.readFastqLine(handle)

        if not self.deflineStringQual: # Qual file can end early
            self.deflineQual.reset()
        else:
            self.deflineStringQual = self.deflineStringQual.strip()

            # Parse defline and confirm validity
            # (needs more work to match up with seq defline or wait for matching seq defline)

            self.deflineQual.parseDeflineString( self.deflineStringQual )

            if ( ( not self.deflineQual.isValid ) or
                 ( not self.ignoreNames and
                   not self.ignAndDiscardNames and
                   self.defline.name != self.deflineQual.name ) ):
                if self.outputStatus:
                    self.statusWriter.outputWarning("Discarding line {} in processQualDefline while looking for valid and matching qual defline in file {}"
                                                    .format(self.lineCountQual,self.qualFilename) )
                attempts = 0
                while True:
                    attempts += 1
                    self.deflineStringQual = self.readFastqLine(handle)
                    if not self.deflineStringQual: # Qual file can end early
                        self.deflineQual.reset()
                        return # End of truncated qual file, continue with seq
                    else:
                        self.deflineStringQual = self.deflineStringQual.strip()
                        self.deflineQual.parseDeflineString( self.deflineStringQual )
                        if ( self.deflineQual.isValid and
                             self.defline.name == self.deflineQual.name ): # Check for equal deflines in case seq was skipped for some reason
                            break
                        elif ( ( self.maxSearchCount and attempts == self.maxSearchCount ) or
                               self.ignoreNames or
                               self.ignAndDiscardNames ):
                            if not self.isSeqQualCheck:
                                self.statusWriter.outputErrorAndExit( "Unable to find valid or matching quality defline after {} lines in filename {} at lineCount {}"
                                                                      .format(attempts,self.qualFilename,self.lineCountQual) )
                            break
                        if self.outputStatus:
                            self.statusWriter.outputWarning("Discarding line {} in processQualDefline while looking for valid and matching qual defline in file {}"
                                                            .format(self.lineCountQual,self.qualFilename) )

            if ( not self.ignoreNames and
                 not self.ignAndDiscardNames and
                 self.deflineStringQual and # Allowing qual file to end early
                 self.defline.name != self.deflineQual.name ):
                self.seqQualDeflineMismatch = True

    ############################################################
    # Process quality related lines
    ############################################################
    def processQuality(self):
        self.processQualDefline(self.qualHandle)
        if self.deflineStringQual:
            self.lengthQual = self.readQual( self.qualHandle )

    ############################################################
    # Check if two handles represent seq/qual fastq pair
    ############################################################

    def isSeqQualFastq (self):
        self.restart()
        self.isSeqQualCheck = True
        self.checkSeqQual = False
        self.read()
        self.restart()

        parsedSeq1 = Seq(self.seqOrig)
        parsedSeq2 = Seq(self.qual)
        if ( self.seq and   # In case qual file comes before seq file and self.badFlatNumQual is True
             self.badFlatNumQual and
             " " in self.seq):
            self.seq = self.qualParser.fixBadFlatNumQual(self.seq)
        parsedQual1 = Qual(self.seq,parsedSeq2.length)
        parsedQual2 = Qual(self.qual,parsedSeq1.length)

        if Defline.isPairedDeflines ( self.defline, self.deflineQual, True ):
            if parsedSeq1.isValid and parsedSeq2.isValid: # 2 fasta files with the same name
                return False
            if parsedSeq1.isValid and parsedQual2.isValid:
                return 1
            elif parsedSeq2.isValid and parsedQual1.isValid:
                return 2

        return False

    ############################################################
    # EOF check (orig strings include linefeeds)
    ############################################################
    def isEof (self):
        return not self.seqOrig

############################################################
# Fasta only Fastq Reader with generated qual
############################################################

class FastaFastqReader (FastqReader):
    """ Fasta-only variation loaded as Fastq """

    def __init__(self,filename,handle):
        FastqReader.__init__(self, filename, handle)
        self.deflineCharSeq = ">"
        self.isFasta = True

    ############################################################
    # Generate quality for fasta submission
    ############################################################

    def processQuality(self):
        if self.length != 0:
            self.qualOrig = self.defaultQual * self.length
            self.qual = self.qualOrig
        self.lengthQual = self.length

    ############################################################
    # Check if file is fasta
    ############################################################

    def isFastaFastq(self):
        if self.genbankFile:
            return True
        self.restart()
        self.isMultiLine = True # I don't know if it is multiline fasta at this point but it doesn't hurt and could be required
        self.checkSeqQual = False
        self.read()
        parsedSeq = Seq(self.seqOrig) # Use seqOrig here because cskey is stripped on just 'seq', violating validity check
        self.read()
        if self.eof: # Special case single sequence fasta file
            self.restart()
            return parsedSeq.isValid
        else:
            parsedSeq2 = Seq(self.seqOrig)
            self.restart()
            return parsedSeq.isValid and parsedSeq2.isValid

    ############################################################
    # EOF check
    ############################################################
    def isEof (self):
        return ( ( not self.seqOrig and
                   not self.savedDeflineStringSeq ) or # Potential for bad seq between deflines leading to empty seqOrig
                 ( self.concatOne and
                   self.prevFilePos >= self.concatPos ) )

    ############################################################
    # Check if fasta file contains multi-line sequence
    ############################################################

    def isMultiLineFasta ( self ):
        self.seqLineCount = 0
        maxSeqLineCount = 0

        if self.genbankFile:
            return True

        self.restart()
        while ( self.spotCount < 1001 and
                not ( ( self.spotCount > 0 and
                        not self.seqOrig ) or
                       self.seqLineCount > self.maxSeqLineCountAllowed ) ):
            self.read()
            if self.seqLineCount > maxSeqLineCount:
                maxSeqLineCount = self.seqLineCount
        self.restart()

        if 1 < maxSeqLineCount <= self.maxSeqLineCountAllowed:
            return True

        if maxSeqLineCount > self.maxSeqLineCountAllowed:
            self.statusWriter.outputErrorAndExit( "Distance between fasta deflines exceeds {} lines in file ... {}"
                                                  .format(self.maxSeqLineCountAllowed,self.filename) )

        else:
            return False

############################################################
# Single Line Fastq Reader Class
############################################################

class SingleLineFastqReader (FastqReader):
    """ Single line variation on FastqReader (colon-separated fields) """

    def __init__(self,filename,handle):
        FastqReader.__init__(self, filename, handle)
        self.deflineCharSeq = ''
        self.deflineChunks = []
        self.fileString = ''
        self.delim = ":"
        self.filter = 0
        self.singleLineCheck = False
        self.nameColumns = []
        self.seqColumns = []
        self.qualColumns = []
        self.groupColumns = []
        self.spotGroupString = ''
        self.isQseq = False
        self.isExport = False
        self.handle = handle

    ############################################################
    # Read/parse single line fastq from file
    ############################################################

    def read (self):

        attemptCount = 0

        while True:

            self.resetFastqValues()
            attemptCount += 1

            # Read line from file

            fileStringOrig = self.readFastqLine(self.handle)
            self.fileString = fileStringOrig.strip()

            # Check for EOF

            if not fileStringOrig:
                self.resetFastqValues()

            # Split file line into chunks

            else:
                self.processLineChunks()
                if ( not self.singleLineCheck and
                     ( not self.deflineStringSeq or
                       not self.defline.parseDeflineString( self.deflineStringSeq ) ) ):
                    self.findValidDefline(False)
                if ( self.spotGroupString and
                     not self.defline.spotGroup ):
                    self.defline.spotGroup = self.spotGroupString

            # Break if not empty seq or eof; otherwise discard and try to read another read/spot

            if self.singleLineCheck or self.seq or self.isEof():
                break
            else:
                if ( ( self.maxSearchCount != 0 and attemptCount > self.maxSearchCount ) or
                     not fileStringOrig ):
                    self.statusWriter.outputErrorAndExit( "Unable to find valid spot in {} lines or before eof ... \n\tfilename ... {}\n\tlineCount ... {}"
                                                          .format(self.maxSearchCount,self.filename,self.lineCount) )
                self.discardCount += 1
                if ( self.outputStatus and
                     ( not self.maxDiscardCount or
                       self.discardCount < self.maxDiscardCount ) ):
                    self.statusWriter.outputWarning("Discarding line {} in SingleLineFastqReader spot read due to an invalid spot in file {}"
                                                    .format(self.lineCount,self.filename) )
                elif ( self.maxDiscardCount and
                       self.discardCount > self.maxDiscardCount ):
                    self.statusWriter.outputErrorAndExit( "Excessive errors: {} cases occurred where lines were discarded from {} at line {}"
                                                          .format(self.discardCount,self.filename,self.lineCount) )

        # Increment spot count and check for eol

        self.updateSpotCountOrEof()

    ############################################################
    # Process single line chunks
    ############################################################

    def processLineChunks (self):

        self.deflineChunks = self.fileString.split(self.delim)

        # Get filter value if qseq or export format

        if ( not self.nameColumns and
             self.delim == '\t' ):
            filterVal = self.deflineChunks.pop()
            if filterVal in 'NY':
                self.isExport = True
            else:
                self.isQseq = True

            if filterVal in '0N':     # Qseq format indicating did not pass filtering (0) OR
                                      # Export format indicating did not pass filtering (N) (http://allaboutbioinfo.blogspot.com/2011/08/qseq-and-export-file-format-of-illumina.html)
                self.filter = 'Y'     # 'Y' indicates 'yes', the read should be filtered
            else:
                self.filter = 'N'

        # Process qual chunk

        self.lengthQual = self.readQual( None )

        # Process seq chunk

        self.length = self.readSeq()

        # Rebuild/parse defline string, confirming validity

        self.buildDeflineString()

    ############################################################
    # Process seq defline
    ############################################################

    def buildDeflineString (self):
        if self.seqColumns:
            self.deflineStringSeq = '@'
            if self.nameColumns == -1:
                self.deflineStringSeq += str(self.spotCount + 1)
            else:
                if self.spotNumAtNameStart:
                    self.deflineStringSeq += str(self.spotCount + 1) + '_'
                colNum = len ( self.nameColumns )
                count = 0
                for nameCol in self.nameColumns:
                    if len(self.deflineChunks) < nameCol:
                        self.deflineStringSeq = ''
                        break
                    else:
                        count += 1
                        if ( self.deflineStringSeq == '@' and
                             (self.deflineChunks[nameCol-1][0] in '@>') ):
                            self.deflineChunks[nameCol-1] = self.deflineChunks[nameCol-1][1:]
                        if count != colNum:
                            self.deflineStringSeq += str(self.deflineChunks[nameCol-1]) + '_'
                        else:
                            self.deflineStringSeq += self.deflineChunks[nameCol-1]

            if (self.deflineStringSeq and
                self.groupColumns):
                groupLen = self.readGroup()
                if groupLen > 0:
                    self.deflineStringSeq += ' ' + self.spotGroupString

        elif self.delim == '\t': #qseq format
            if len(self.deflineChunks) < 3:
                self.resetFastqValues()
            else:
                readNum = self.deflineChunks.pop()
                index = self.deflineChunks.pop()
                self.deflineStringSeq = '@'
                self.deflineStringSeq += ':'.join(self.deflineChunks)
                self.deflineStringSeq += " {}:{}:0:{}".format(readNum,self.filter,index)
        else:
            if not self.deflineChunks:
                self.resetFastqValues()
            else:
                self.deflineStringSeq = '@' + self.delim.join(self.deflineChunks)

    ############################################################
    # Read seq and convert to uppercase
    ############################################################

    def readSeq (self):
        if self.seqColumns:
            self.seqOrig = ''
            for seqCol in self.seqColumns:
                if len(self.deflineChunks) < seqCol:
                    return 0
                else:
                    self.seqOrig += self.deflineChunks[seqCol-1]
            self.seq = self.seqOrig.upper()
            if not self.qualColumns:
                self.qualOrig = self.defaultQual * len(self.seq)
                self.qual = self.qualOrig
                self.lengthQual = len(self.seq)
        elif not self.deflineChunks:
            return 0
        else:
            self.seqOrig = self.deflineChunks.pop()
            self.seq = self.seqOrig.upper()

        if self.outputStatus: # If I do this too early it messes up other checks (i.e. -1 in quality is changed to N1)
            self.seq = self.seq.translate(self.transDashUridineX)

        if self.seqParser.parseSeq( self.seq, self.setClips ):
            return len(self.seq)
        else:
            self.seq = ""
            return 0

    ############################################################
    # Read qual (handle,defline not used - retained to meet arg needs)
    ############################################################

    def readQual (self,handle):
        if self.seqColumns:
            if not self.qualColumns:
                return 0
            else:
                self.qualOrig = ''
                for qualCol in self.qualColumns:
                    if len(self.deflineChunks) < qualCol:
                        return 0
                    else:
                        self.qualOrig += self.deflineChunks[qualCol-1]
                self.qual = self.qualOrig
                return len ( self.qual )
        elif not self.deflineChunks:
            return 0
        else:
            if self.isExport:
                for i in range(0,11):
                    self.deflineChunks.pop()
            self.qualOrig = self.deflineChunks.pop()
            self.qual = self.qualOrig
            if ( self.qual and
                 self.qualParser.parseQual( self.qual, len (self.qual ) ) ):
                if self.isNumQual:
                    return len ( self.qual.split() )
                else:
                    return len ( self.qual )
            else:
                self.qual = ""
                return 0

    ############################################################
    # Read group columns
    ############################################################

    def readGroup (self):
        if self.groupColumns:
            self.spotGroupString = ''
            for groupCol in self.groupColumns:
                if len(self.deflineChunks) < groupCol:
                    return len ( self.spotGroupString )
                elif len ( self.deflineChunks[groupCol-1] ) == 0:
                    continue
                elif not self.spotGroupString:
                    self.spotGroupString = self.deflineChunks[groupCol-1]
                else:
                    self.spotGroupString += '_' + self.deflineChunks[groupCol-1]
        return len ( self.spotGroupString )

    ############################################################
    def resetFastqValues ( self ):
        self.seqOrig = ''
        self.seq = ''
        self.qualOrig = ''
        self.qual = ''
        self.deflineStringSeq = ''
        self.spotGroupString = ''
        self.defline.reset()

    ############################################################
    # Find valid defline
    ############################################################

    def findValidDefline (self,callProcessSeqDefline):
        findAttemptCount = 0
        while True:
            findAttemptCount += 1
            if self.fileString:
                self.discardCount += 1
                if self.outputStatus:
                    if ( not self.maxDiscardCount or
                         self.discardCount < self.maxDiscardCount ):
                        self.statusWriter.outputWarning("Discarding line {} in SingleLineFastqReader findValidDefline while looking for a valid defline in file {}"
                                                        .format(self.lineCount,self.filename) )
                    elif ( self.maxDiscardCount and
                           self.discardCount > self.maxDiscardCount ):
                        if self.singleLineCheck:
                            return
                        else:
                            self.statusWriter.outputErrorAndExit( "Excessive errors: {} cases occurred where lines were discarded from {} at line {}"
                                                                  .format(self.discardCount,self.filename,self.lineCount) )

            # Read line from file

            fileStringOrig = self.readFastqLine(self.handle)
            self.fileString = fileStringOrig.strip()

            # Check for EOF

            if not fileStringOrig:
                self.resetFastqValues()
                break

            # Split file line into chunks

            else:
                self.processLineChunks()
                if ( self.deflineStringSeq and
                     self.defline.parseDeflineString( self.deflineStringSeq ) ):
                    break

            if ( ( self.maxSearchCount != 0 and findAttemptCount > self.maxSearchCount ) or
                 not fileStringOrig):
                if self.singleLineCheck:
                    return
                else:
                    self.statusWriter.outputErrorAndExit( "Unable to parse defline in {} lines or before eof ... \n\tfilename ... {}\n\tlineCount ... {}"
                                                          .format(self.maxSearchCount,self.filename,self.lineCount) )

    ############################################################
    # Check if fastq file contains single line fastq
    ############################################################

    def isSingleLineFastq (self,skipLines):
        self.singleLineCheck = True
        self.restart()
        self.checkSeqQual = False
        for i in range(0, skipLines):
            self.handle.readline()
        self.read()
        parsedSeq = Seq(self.seq)
        parsedQual = Qual( self.qual, parsedSeq.length )
        self.restart()
        if ( parsedSeq.isValid and
             parsedQual.isValid and
             parsedSeq.length == parsedQual.length and
             parsedSeq.length > 1 and
             parsedQual.length > 1):
            self.singleLineCheck = False
            return True
        else:
            self.delim = '\t'
            self.restart()
            self.checkSeqQual = False
            for i in range(0, skipLines):
                self.handle.readline()
            self.read()
            self.singleLineCheck = False
            parsedSeq = Seq(self.seq)
            parsedQual = Qual( self.qual, parsedSeq.length )
            return bool( parsedSeq.isValid and
                         parsedQual.isValid and
                         parsedSeq.length == parsedQual.length and
                         parsedSeq.length > 1 and
                         parsedQual.length > 1 )

############################################################
# Class for writing fastq-based spots to SRA archive
############################################################

class FastqSpotWriter:
    """ Container for collecting/writing fasta spot information to gw """

    READ_TYPE_TECHNICAL                 = 0
    READ_TYPE_BIOLOGICAL                = 1
    READ_TYPE_GROUP                     = 2

    def __init__(self):

        self.readCount = 0
        self.readCountNonGroup = 0
        self.readLengths = []
        self.readStarts = []
        self.readTypes = array.array('B', [1, 1])
        self.platform = array.array('B', [0])
        self.platformString = ''
        self.fileType = ''
        self.readTypeString = ''

        self.clipQualityLeft = array.array('I', [0])
        self.clipQualityRight = array.array('I', [0])
        self.labels = ''
        self.firstLabelRE = ''
        self.labelStarts = array.array('I', [ 0, 0 ])
        self.labelLengths = array.array('I', [ 0, 0 ])
        self.offset = 33
        self.minQual = 1000
        self.maxQual = 0
        self.isNumQual = False
        self.orphanReads = False
        self.dumpOrphans = False
        self.dumpOrphans2D = False
        self.pairedRead1 = {}
        self.pairedRead2 = {}
        self.foundTemplate = False
        self.foundComplement = False
        self.found2D = False
        self.nanopore2Dread = {}
        self.nanopore2Donly = False
        self.nanoporeTemplateOnly = False
        self.nanoporeTemplateComplementOnly = False
        self.nanoporeNoLabels = False
        self.nanopore_one_type_only = False
        self.appendPoreReadToName = False
        self.read1PairFiles = None
        self.read2PairFiles = None
        self.read3PairFiles = None
        self.read4PairFiles = None
        self.read5PairFiles = None
        self.read6PairFiles = None
        self.read1QualFiles = None
        self.read2QualFiles = None
        self.filesToIgnore = None
        self.filesToTruncate = None
        self.fileSpotGroupsProvided = None
        self.useFilenameForSG = False
        self.nameColumns = None
        self.seqColumns = None
        self.qualColumns = None
        self.groupColumns = None
        self.spotGroup = ''
        self.spotGroupsFound = {}
        self.discardBarcodes = False
        self.discardBCchars = 0
        self.useFileOrder = 0
        self.appendBCtoName = False
        self.bcSpaceOffset = False
        self.concatPairFiles = False
        self.removeBadChars = False
        self.removeSeqSpaces = False
        self.ignoreNames = False
        self.ignAndDiscardNames = False
        self.useAndDiscardNames = False
        self.nameCount = 0
        self.isEightLine = False
        self.mixedDeflines = False
        self.genericDeflines = False
        self.mixedTypes = False
        self.setClips = False
        self.retainAltNum = False
        self.badFlatNumQual = False
        self.convertEmptyDeflines = False
        self.isMultiLine = False
        self.duplicateReads = False
        self.convertDeflineSpaces = False
        self.addFileToName = False
        self.filenameTrunc = None
        self.keepEmptySeqs = False
        self.spotNumAtNameStart = False # Only used for single line fastq in columns
        self.allowEarlyFileEnd = False # Only used for 3 or more files
        self.seqOnlyFile = False
        self.genbankFile = False
        self.noPairFileCheck = False
        self.ignoreSGforPairing = False
        self.fastqWithGtrThan = False
        self.fastqWithoutAtSign = False
        self.removeLineQuotes = False
        self.removeLastChar = False
        self.abiLastPrimerBase = None
        self.foundIndex = None
        self.foundI1 = None
        self.foundR1 = None
        self.foundR2 = None

        self.dst = None
        self.dst2D = None

        self.lengthsProvided = False
        self.variableLengthSpecified = False
        self.labelsProvided = False
        self.typesProvided = False
        self.spotGroupProvided = False
        self.offsetProvided = False
        self.pairFilesProvided = False
        self.qualFilesProvided = False
        self.spacesInFilenames = False
        self.dashFileNames = False

        self.useSharq = False

        self.workDir = None
        self.outdir = "out.sra"
        self.finalDest = ''
        self.spotCount = 0
        self.totalSize = 0
        self.totalSize2 = 0
        self.totalSize3 = 0
        self.totalSize4 = 0
        self.isCompressed = 0
        self.cumulativeBytesRead = 0
        self.cumulativeBytesRead2 = 0

        self.isColorSpace = False       # ABI fastq rather than base-based fastq -
                                        # handled by abi-load currently; comment lines at the start
        self.extraQualForCSkey = False  # Set to True if extra qual value for cs key
        self.fabricateNames = False     # Cases of no name or all the same name will result in the spot number as the name for CS
        self.readColumn = 'READ'

        self.logOdds = False            # Indicated by presence of negative qualities
        self.changeNegOneQual = False   # '-1' only is likely used for dot or N qualities
        self.transNegOne = str.maketrans('?', '@')
        self.readNums = []

        self.gw = None
        self.db = 'NCBI:SRA:GenericFastq:db'
        self.schema = 'sra/generic-fastq.vschema'
        self.maxPairCheckCount = 250000
        self.maxErrorCount = 100
        self.maxSearchCount = 100
        self.maxErrorOutput = 100
        self.maxSeqLineCountAllowed = 20000
        self.maxDeflineLen = 1000
        self.errorCount = 0
        self.infoSizeForSpotCount = 1000000 # Determines frequency of progress notifications
        self.profile = False
        self.checkSeqQual = True
        self.statusWriter = None
        self.ignLeadCharsNum = None
        self.ignTailCharsNum = None
        self.ignCharsAfterPlus = False
        self.forceLoad = False
        self.writingToArchive = False

        ############################################################
        # Define SEQUENCE table components for general loader
        # and create instance of General Writer
        ############################################################

        self.tbl = {
            'SEQUENCE': {
                'READ': {
                    'expression': '(INSDC:dna:text)READ',
                    'elem_bits': 8
                    },
                'CSREAD': {
                    'expression': '(INSDC:color:text)CSREAD',
                    'elem_bits': 8
                    },
                'CS_KEY': {
                    'expression': '(INSDC:dna:text)CS_KEY',
                    'elem_bits': 8
                    },
                'READ_START': {
                    'expression': '(INSDC:coord:zero)READ_START',
                    'elem_bits': 32
                    },
                'READ_LENGTH': {
                    'expression': '(INSDC:coord:len)READ_LEN',
                    'elem_bits': 32
                    },
                'READ_TYPE': {
                    'expression': '(U8)READ_TYPE',
                    'elem_bits': 8,
                    'default': array.array('B', [ 1, 1 ])
                    },
                'READ_FILTER': {
                    'expression': '(U8)READ_FILTER',
                    'elem_bits': 8,
                    'default': array.array('B', [ 0, 0 ])
                    },
                'QUALITY': {
                    'expression': '(INSDC:quality:text:phred_33)QUALITY',
                    'elem_bits': 8
                    },
                'NAME': {
                    'expression': '(ascii)NAME',
                    'elem_bits': 8,
                    'default': ' '.encode('ascii')
                    },
                'SPOT_GROUP': {
                    'elem_bits': 8,
                    'default': ''.encode('ascii')
                    },
                'CLIP_QUALITY_LEFT': {
                    'expression': '(INSDC:coord:one)CLIP_QUALITY_LEFT',
                    'elem_bits': 32,
                    'default': array.array('I', [0])
                    },
                'CLIP_QUALITY_RIGHT': {
                    'expression': '(INSDC:coord:one)CLIP_QUALITY_RIGHT',
                    'elem_bits': 32,
                    'default': array.array('I', [0])
                    },
                'LABEL': {
                    'elem_bits': 8,
                    'default': ''.encode('ascii')
                    },
                'LABEL_START': {
                    'elem_bits': 32,
                    'default': array.array('I', [ 0, 0 ])
                    },
                'LABEL_LEN': {
                    'elem_bits': 32,
                    'default': array.array('I', [ 0, 0 ])
                    },
                'PLATFORM': {
                    'expression': '(U8)PLATFORM',
                    'elem_bits': 8,
                    'default': array.array('B', [0])
                    },
                'CHANNEL': {
                    'elem_bits': 32,
                    'default': array.array('I', [0])
                    },
                'READ_NO': {
                    'expression': 'READ_NUMBER',
                    'elem_bits': 32,
                    'default': array.array('I', [0])
                    },
                },
            'CONSENSUS': {
                'READ': {
                    'expression': '(INSDC:dna:text)READ',
                    'elem_bits': 8,
                    },
                'QUALITY': {
                    'expression': '(INSDC:quality:text:phred_33)QUALITY',
                    'elem_bits': 8,
                    },
                'NAME': {
                    'expression': '(ascii)NAME',
                    'elem_bits': 8,
                    'default': ' '.encode('ascii')
                    },
                'SPOT_GROUP': {
                    'elem_bits': 8
                    },
                'CHANNEL': {
                    'elem_bits': 32
                    },
                'READ_NO': {
                    'expression': 'READ_NUMBER',
                    'elem_bits': 32,
                    'default': array.array('I', [0])
                    },
                'READ_START': {
                    'expression': '(INSDC:coord:zero)READ_START',
                    'elem_bits': 32,
                    'default': array.array('I', [0])
                    },
                'READ_LENGTH': {
                    'expression': '(INSDC:coord:len)READ_LEN',
                    'elem_bits': 32,
                    },
                'READ_TYPE': {
                    'expression': '(U8)READ_TYPE',
                    'elem_bits': 8,
                    'default': array.array('B', [1])
                    },
                'READ_FILTER': {
                    'expression': '(U8)READ_FILTER',
                    'elem_bits': 8,
                    'default': array.array('B', [0])
                    },
                'LABEL': {
                    'elem_bits': 8,
                    'default': '2D'.encode('ascii')
                    },
                'LABEL_START': {
                    'elem_bits': 32,
                    'default': array.array('I', [0])
                    },
                'LABEL_LEN': {
                    'elem_bits': 32,
                    'default': array.array('I', [2])
                    },
                }
            }

    ############################################################
    # Open general writer for output to general loader
    ############################################################

    def openGeneralWriter(self):

        # Address nanopore possibilities

        if self.platformString == "NANOPORE":
            self.db = 'NCBI:SRA:GenericFastqNanopore:db'
            if ( self.nanopore2Donly or
                 self.nanoporeNoLabels or
                 self.nanoporeTemplateOnly ):
                self.nanopore_one_type_only = True
                del self.tbl['CONSENSUS']
            elif self.nanoporeTemplateComplementOnly:
                self.setReadLabels("template,complement")
                self.orphanReads = True
                del self.tbl['CONSENSUS']
            else:
                self.setReadLabels("template,complement")
                self.orphanReads = True
        else:
            self.removeNanoporeColumns()
            del self.tbl['CONSENSUS']

        # Address color space possibilities

        if self.isColorSpace:
            self.db = 'NCBI:SRA:GenericFastqAbsolid:db'
            self.removeNonColorSpaceColumns()
            self.readColumn = 'CSREAD'

        elif 'SEQUENCE' in self.tbl:
            self.removeColorSpaceColumns()

        # Address other special cases (not overlapping with nanopore and color space)

        if self.logOdds:
            self.db = 'NCBI:SRA:GenericFastqLogOdds:db'

        elif ( self.ignAndDiscardNames or
               self.useAndDiscardNames ):
            self.db = 'NCBI:SRA:GenericFastqNoNames:db'

        # Remove label and clip columns if necessary

        if 'SEQUENCE' in self.tbl:

            # Remove labels from SEQUENCE table if not provided

            if not self.labels:
                self.removeLabelColumns()

            # Remove clips from SEQUENCE table if not needed

            if not self.setClips:
                self.removeClipColumns()

        # Set destinations

        if 'CONSENSUS' in self.tbl:
            self.dst2D = self.tbl['CONSENSUS']

        if 'SEQUENCE' in self.tbl:
            self.dst = self.tbl['SEQUENCE']

        # Open writer

        self.gw = GeneralWriter.GeneralWriter( self.outdir,
                                               self.schema,
                                               self.db,
                                               'fastq-load.py',
                                               '1.1.2',
                                               self.tbl )

    ############################################################
    # Set values that are hopefully consistent across a file
    # Reduces the write time
    ############################################################

    def setUnchangingSpotValues(self,fastq2=None,fastq3=None):

        if self.platformString:
            self.dst['PLATFORM']['data'] = self.platform

        # Handle recognized nanopore data somewhat differently

        if self.platformString == "NANOPORE" :
            if self.spotCount == 0:
                if self.nanopore_one_type_only:
                    self.dst['READ_START']['data'] = array.array( 'I', [ 0 ] )
                elif not self.nanoporeTemplateComplementOnly:
                    self.dst2D['READ_START']['data'] = array.array( 'I', [ 0 ] )

        # Orphan reads trump setting constant values

        elif self.orphanReads or self.duplicateReads:
            pass

        # Fragment file

        elif ( not fastq2 and
               not self.lengthsProvided ):
            self.dst['READ_START']['data'] = array.array( 'I', [ 0 ] )
            self.dst['READ_TYPE']['data'] = array.array( 'B', [ self.READ_TYPE_BIOLOGICAL ] )

        # Multiple reads

        else:
            if ( fastq3 and
                 len ( self.readTypes ) < 3 ):
                self.statusWriter.outputErrorAndExit("If 3 or more reads are provided for the same spots for other than nanopore data, then the 'readTypes' option must be utilized.")

            self.dst['READ_TYPE']['data'] = array.array( 'B', [] )
            if ( self.lengthsProvided and
                 not self.variableLengthSpecified ):
                self.dst['READ_START']['data'] = array.array( 'I', [] )
                self.dst['READ_LENGTH']['data'] = array.array( 'I', [] )

            readIndex = -1
            readStart = 0
            for readType in self.readTypes:
                readIndex += 1
                if readType != self.READ_TYPE_GROUP:
                    self.dst['READ_TYPE']['data'].append( readType )
                    if ( self.lengthsProvided and
                         not self.variableLengthSpecified ):
                        self.dst['READ_START']['data'].append( readStart )
                        self.dst['READ_LENGTH']['data'].append( self.readLengths[readIndex] )
                        readStart += self.readLengths[readIndex]

        # Can be called twice if fastq1 or fastq2 ends early
        # Execute the following lines only once

        if self.spotCount == 0:

            if self.labels :
                self.dst['LABEL']['data'] = self.labels
                self.dst['LABEL_START']['data'] = self.labelStarts
                self.dst['LABEL_LEN']['data'] = self.labelLengths

    ############################################################
    # Selects write function for handling fastq data
    ############################################################

    def writeSpot ( self, fastq1, fastq2, fastq3, fastq4, fastq5, fastq6 ):

        fastqs = tuple ( [_f for _f in [ fastq1, fastq2, fastq3, fastq4, fastq5, fastq6 ] if _f] )

        # Very occasionally submitted file contains empty seq/qual for a spot.

        seqFound = False
        for fastq in fastqs:
            if fastq and fastq.seq:
                seqFound = True
                break

        if not seqFound:
            pass

        # Process nanopore spots

        elif self.platformString == "NANOPORE":
            if self.nanopore_one_type_only:
                self.write2Dread ( fastq1, None, self.dst )
            elif ( self.dumpOrphans or
                   self.dumpOrphans2D ) :
                self.processOrphanNanoporeReads ( fastq1 )
            else:
                self.processNanoporeReads ( fastq1, fastq2, fastq3 )

        # Pre-split 454 would be handled here. Assumes no 454 adapter
        # is present (i.e. assuming no sub-division of sequence in each
        # pair file)

        elif self.orphanReads:
            if self.dumpOrphans:
                self.processOrphanReads ( fastq1 )
            else:
                self.processMixedReads ( fastq1, fastq2 )

        else:

            if ( len(fastqs) > 1 or
                 not self.lengthsProvided ):
                self.processFastqSpot ( fastqs )

            else:
                self.processMultiReadSingleFastqSpot ( fastq1 )

            # Set common columns and track spotCount

            self.spotCount += 1
            self.setDstName ( fastq1, self.dst )
            self.setDstSpotGroup ( fastqs, self.dst )

            # Output spotcount

            if fastq1 is not None and ( fastq1.spotCount % self.infoSizeForSpotCount ) == 0:
                self.outputCount(False,fastq1,fastq2,fastq3,fastq4,fastq5,fastq6)

            # Edge case where one seq in a pair is an empty string (very sparse)

            if not self.typesProvided and not self.duplicateReads:
                if len(fastqs) == 2:
                    if not fastq1.seq:
                        self.dst['READ_TYPE']['data'] = array.array( 'B', [ self.READ_TYPE_TECHNICAL, self.READ_TYPE_BIOLOGICAL ] )
                    elif not fastq2.seq:
                        self.dst['READ_TYPE']['data'] = array.array( 'B', [ self.READ_TYPE_BIOLOGICAL, self.READ_TYPE_TECHNICAL ] )

            # Finally write to general writer

            self.gw.write(self.dst)

            # Recover from edge case

            if not self.typesProvided and not self.duplicateReads:
                if ( len(fastqs) == 2 and
                     ( not fastq1.seq or
                       not fastq2.seq ) ):
                    self.dst['READ_TYPE']['data'] = array.array( 'B', [ self.READ_TYPE_BIOLOGICAL, self.READ_TYPE_BIOLOGICAL ] )

    ############################################################
    # Process multiple fastq files
    ############################################################

    def processFastqSpot ( self, fastqs ):
        seq = ''
        qual = ''
        vdbStart = 0
        clipLeft = 0
        clipRight = 0
        csKey = ''
        readIndex = -1
        prevName = None
        prevNum = None
        prevSeq = None
        isFiltered = False
        nonGroupReadCount = 0

        if len(fastqs) == 2:
            if ('data' in self.dst['READ_TYPE'] and
                not len(self.dst['READ_TYPE']['data']) == 2): # Don't change types if already present
                self.dst['READ_TYPE']['data'] = array.array( 'B', [ self.READ_TYPE_BIOLOGICAL, self.READ_TYPE_BIOLOGICAL ] )
        elif len(fastqs) == 1:
            self.dst['READ_TYPE']['data'] = array.array( 'B', [ self.READ_TYPE_BIOLOGICAL ] )

        collectLengths = False
        if ( not self.lengthsProvided or
             self.variableLengthSpecified ):
            self.dst['READ_START']['data'] = array.array( 'I', [] )
            self.dst['READ_LENGTH']['data'] = array.array( 'I', [] )
            collectLengths = True

        for fastq in fastqs:

            if ( not ( self.ignoreNames or self.ignAndDiscardNames ) and
                 prevName and
                 prevName != fastq.defline.name ):
                self.statusWriter.outputErrorAndExit("Expected matching deflines but encountered name mismatch at line {} in {} ... {}\tvs\t{}"
                                                     .format(fastq.lineCount,fastq.filename,prevName,fastq.defline.name) )

            # Skip identical reads - leads to fragments in the middle but hashing the whole file is too big

            if ( self.duplicateReads and
                 prevNum and
                 prevNum == fastq.defline.readNum and
                 prevSeq == fastq.seq ):
                self.dst['READ_TYPE']['data'] = array.array( 'B', [ self.READ_TYPE_BIOLOGICAL ] )
                continue

            prevName = fastq.defline.name
            prevNum = fastq.defline.readNum
            prevSeq = fastq.seq

            readIndex += 1
            if self.readTypes[readIndex] != self.READ_TYPE_GROUP:

                nonGroupReadCount += 1

                if collectLengths:
                    self.dst['READ_START']['data'].append( vdbStart )
                    self.dst['READ_LENGTH']['data'].append( len(fastq.seq) )
                    vdbStart += len(fastq.seq)

                if self.isColorSpace:
                    csKey += fastq.csKey

                if self.setClips:
                    clipLeft += fastq.clipLeft
                    clipRight += fastq.clipRight

                seq += fastq.seq

                if ( self.isNumQual and
                     qual and
                     fastq.qual ):
                    qual += ' '
                qual += fastq.qual

                if fastq.defline.filterRead:
                    isFiltered = True

        self.dst[self.readColumn]['data'] = seq.encode('ascii')
        self.setDstQual ( qual, self.dst )
        if isFiltered:
            self.dst['READ_FILTER']['data'] = array.array('B', [1] * nonGroupReadCount )
        else:
            self.dst['READ_FILTER']['data'] = array.array('B', [0] * nonGroupReadCount )

        if self.isColorSpace:
            self.dst['CS_KEY']['data'] = csKey.encode('ascii')

        if self.setClips:
            clipLeft += 1
            clipRight = len(seq) - clipRight
            self.dst['CLIP_QUALITY_LEFT']['data'] = array.array( 'I', [ int(clipLeft) ] )
            self.dst['CLIP_QUALITY_RIGHT']['data'] = array.array( 'I', [ int(clipRight) ] )

    ############################################################
    # Process multiple reads from single fastq file
    # Reads lengths provided and read starts calculated already
    ############################################################

    def processMultiReadSingleFastqSpot ( self, fastq ):
        seq = ""
        qual = ""
        start = 0
        vdbStart = 0
        readIndex = -1

        if self.variableLengthSpecified:
            self.dst['READ_START']['data'] = array.array( 'I', [] )
            self.dst['READ_LENGTH']['data'] = array.array( 'I', [] )

        qualValStrings = None
        if self.isNumQual:
            qualValStrings = fastq.qual.split()

        for readType in self.readTypes:
            readIndex += 1
            readLen = self.readLengths[readIndex]
            if readLen == 0:
                readLen = len ( fastq.seq[start:] ) # last read only currently

            if readType != self.READ_TYPE_GROUP:

                seq += fastq.seq[ start : start + readLen ]

                if self.variableLengthSpecified:
                    self.dst['READ_START']['data'].append( vdbStart )
                    self.dst['READ_LENGTH']['data'].append( readLen )

                if self.isNumQual:
                    qual += " ".join( qualValStrings[ start : start + readLen ] )
                    if readIndex != self.readCountNonGroup - 1:
                        qual += ' '
                else:
                    qual += fastq.qual[ start : start + readLen ]

                vdbStart += readLen

            start += readLen

        self.dst['READ']['data'] = seq.encode('ascii')
        self.setDstQual ( qual, self.dst )
        self.dst['READ_FILTER']['data'] = array.array('B', [fastq.defline.filterRead] * self.readCountNonGroup)
        if self.setClips:
            clipQualityLeft = fastq.clipLeft + 1
            clipQualityRight = len(seq) - fastq.clipRight
            self.dst['CLIP_QUALITY_LEFT']['data'] = array.array( 'I', [ int(clipQualityLeft) ] )
            self.dst['CLIP_QUALITY_RIGHT']['data'] = array.array( 'I', [ int(clipQualityRight) ] )

    ############################################################
    # Process spots from nanopore fastq files (potential orphan reads)
    # Note that if a single file contains all three read types
    # then I need to assume
    ############################################################

    def processNanoporeReads ( self, fastq1, fastq2, fastq3 ):

        if ( fastq2 and
             not fastq2.eof and
             fastq1.defline.name == fastq2.defline.name and
             fastq1.defline.poreRead != "2D" and
             fastq2.defline.poreRead != "2D"):
            self.setNanoporeColumns ( fastq1, self.dst )
            self.writeMixedPair ( fastq1, fastq2 )
            self.writeNanopore2Dread ( fastq1, fastq3 )

        elif ( fastq2 and
               not fastq2.eof and
               fastq3 and
               not fastq3.eof and
             fastq2.defline.name == fastq3.defline.name and
             fastq2.defline.poreRead != "2D" and
             fastq3.defline.poreRead != "2D"):
            if fastq1.defline.name != fastq2.defline.name:
                self.processUnpairedPoreRead ( fastq1, None )
            self.setNanoporeColumns ( fastq2, self.dst )
            self.writeMixedPair ( fastq2, fastq3 )
            if fastq1.defline.name == fastq2.defline.name:
                self.writeNanopore2Dread ( fastq2, fastq1 )
            else:
                self.writeNanopore2Dread ( fastq2 )

        elif ( fastq3 and # unlikely but just in case
               not fastq3.eof and
               fastq1.defline.name == fastq3.defline.name and
               fastq1.defline.poreRead != "2D" and
               fastq3.defline.poreRead != "2D"):
            self.setNanoporeColumns ( fastq1, self.dst )
            self.writeMixedPair ( fastq1, fastq3 )
            self.writeNanopore2Dread ( fastq1, fastq2 )

        else:
            self.processUnpairedPoreRead ( fastq1, fastq3 )
            if fastq2 :
                self.processUnpairedPoreRead ( fastq2, None ) # 2D read in fastq3 handled/cached by prior call

        if fastq1 is not None and ( fastq1.spotCount % self.infoSizeForSpotCount == 0 ):
            self.outputCount(True,fastq1,fastq2)

    ############################################################
    # If a read is unpaired then check for prior encounter with
    # the mate
    ############################################################

    def processUnpairedPoreRead (self, fastq, fastq2D ):

        # Occasionally seeing duplicates

        if ( fastq.defline.poreRead == "2D" or
             ( not fastq.defline.name in self.pairedRead1 and
               not fastq.defline.name in self.pairedRead2 ) or
             ( fastq.defline.name in self.pairedRead1 and
               fastq.defline.poreRead == "template" ) or
             ( fastq.defline.name in self.pairedRead2 and
               fastq.defline.poreRead == "complement" ) ):
            self.addToPoreReadHash( fastq )
            if ( fastq2D and
                 not fastq2D.eof ):
                self.addToPoreReadHash( fastq2D )

        else:
            self.setNanoporeColumns ( fastq, self.dst )
            if fastq.defline.name in self.pairedRead2:
                self.writeMixedRead1 ( fastq )
            else:
                self.writeMixedRead2 ( fastq )
            self.writeNanopore2Dread ( fastq, fastq2D )

    ############################################################
    # Process orphan reads from nanopore fastq files (probably
    # template only)
    ############################################################

    def processOrphanNanoporeReads ( self, fastq1 ):

        if ( self.dumpOrphans and
             not fastq1.defline.poreRead == "2D" and
             ( fastq1.defline.name in self.pairedRead1 or
               fastq1.defline.name in self.pairedRead2 ) ) :
            self.setNanoporeColumns ( fastq1, self.dst )
            self.writeMixedOrphan ( fastq1 )
            self.writeNanopore2Dread ( fastq1 )
            if fastq1.defline.name in self.pairedRead1:
                del self.pairedRead1 [ fastq1.defline.name ]
            elif fastq1.defline.name in self.pairedRead2:
                del self.pairedRead2 [ fastq1.defline.name ]
        elif ( self.dumpOrphans2D and
               fastq1.defline.poreRead == "2D" and
               fastq1.defline.name in self.nanopore2Dread ) :
            self.writeFakeTemplateRead ( fastq1 )
            self.writeNanopore2Dread ( fastq1 )

        if fastq1 is not None and ( fastq1.spotCount % self.infoSizeForSpotCount == 0 ):
            self.outputCount(True,fastq1)

    ############################################################
    # Write actual or fake 2D read
    ############################################################

    def writeNanopore2Dread ( self, fastq, fastq2D=None ):

        readWritten = False
        if ( fastq2D and
             not fastq2D.eof ):

            if ( fastq2D.defline.poreRead != "2D" or
                 fastq.defline.name != fastq2D.defline.name):
                self.processUnpairedPoreRead ( fastq2D, None )
            else:
                self.write2Dread ( fastq2D, None, self.dst2D )
                readWritten = True

        if not readWritten:
            read2D = None
            if fastq.defline.name in self.nanopore2Dread:
                read2D = self.nanopore2Dread [ fastq.defline.name ]

            if read2D:
                self.write2Dread ( fastq, read2D, self.dst2D )
                del self.nanopore2Dread [ fastq.defline.name ]
            elif not self.nanoporeTemplateComplementOnly:
                self.writeFake2Dread ( fastq )

    ############################################################
    # Write nanopore 2D reads to CONSENSUS table using either
    # values read only from 2D read file or combination from
    # template/complement file and save 2D read data
    ############################################################

    def write2Dread (self, fastq, read2D, dst):

        self.setNanoporeColumns ( fastq, dst )
        dst['READ_TYPE']['data'] = array.array( 'B', [ 1 ] )
        if read2D:
            dst['READ']['data'] = read2D['seq'].encode('ascii')
            dst['READ_LENGTH']['data'] = array.array( 'I', [ len(read2D['seq']) ] )
            self.setDstQual ( read2D['qual'], dst )
            dst['READ_FILTER']['data'] = array.array('B', [ read2D['filterRead'] ] )
            if read2D['suffix']:
                fastq.defline.suffix = read2D['suffix'] # Allowing for different suffices between template/complement and 2d reads
        else:
            dst['READ']['data'] = fastq.seq.encode('ascii')
            dst['READ_LENGTH']['data'] = array.array( 'I', [ len(fastq.seq) ] )
            self.setDstQual ( fastq.qual, dst )
            dst['READ_FILTER']['data'] = array.array('B', [ fastq.defline.filterRead ] )
        self.setDstName ( fastq, dst )
        self.setDstSpotGroup ( fastq, dst )
        self.gw.write( dst )
        if self.nanopore_one_type_only:
            self.spotCount += 1
            if fastq is not None and ( fastq.spotCount % self.infoSizeForSpotCount == 0 ):
                self.outputCount(False,fastq)

    ############################################################
    # Write technical nanopore template read to SEQUENCE table using
    # 2D fastq read values
    ############################################################

    def writeFakeTemplateRead (self, fastq):

        self.setNanoporeColumns ( fastq, self.dst )
        self.dst['READ_LENGTH']['data'] = array.array( 'I', [0] )
        self.dst['READ_FILTER']['data'] = array.array('B', [ fastq.defline.filterRead ] )
        self.dst['READ_TYPE']['data'] = array.array( 'B', [ 0 ] )
        self.dst['READ']['data'] = ''.encode('ascii')
        self.setDstQual ( '', self.dst )
        self.setDstName ( fastq, self.dst )
        self.setDstSpotGroup ( fastq, self.dst )
        self.gw.write( self.dst )

    ############################################################
    # Write technical nanopore 2D read to CONSENSUS table using
    # fastq read results from template or complement read
    ############################################################

    def writeFake2Dread (self, fastq):

        self.setNanoporeColumns ( fastq, self.dst2D )
        self.dst2D['READ_LENGTH']['data'] = array.array( 'I', [0] )
        self.dst2D['READ_FILTER']['data'] = array.array('B', [ fastq.defline.filterRead ] )
        self.dst2D['READ_TYPE']['data'] = array.array( 'B', [ 0 ] )
        self.dst2D['READ']['data'] = ''.encode('ascii')
        self.setDstQual ( '', self.dst2D )
        self.setDstName ( fastq, self.dst2D )
        self.setDstSpotGroup ( fastq, self.dst2D )
        self.gw.write( self.dst2D )

    ############################################################
    # Process spots from mixed fastq files (potential orphan reads)
    ############################################################

    def processMixedReads ( self, fastq1, fastq2 ):

        if ( fastq2 and
             fastq1.defline.name == fastq2.defline.name ):
            self.writeMixedPair ( fastq1, fastq2 )
        else:
            self.processUnpairedRead ( fastq1 )
            if fastq2 :
                self.processUnpairedRead ( fastq2 )

        if fastq1 is not None and ( fastq1.spotCount % self.infoSizeForSpotCount == 0 ):
            self.outputCount (True,fastq1,fastq2)

    ############################################################
    # Write spot when read matches up with previously encountered reads
    # Otherwise retain/hash the read
    ############################################################

    def processUnpairedRead ( self, fastq ):

        if fastq.defline.name in self.pairedRead2:
            self.writeMixedRead1 ( fastq )
        elif fastq.defline.name in self.pairedRead1:
            self.writeMixedRead2 ( fastq )
        else:
            self.addToPairedReadHash ( fastq )

    ############################################################
    # Write spot where read1 and read2 in fastq
    ############################################################

    def writeMixedPair ( self, fastq1, fastq2 ):
        if ( ( fastq1.defline.readNum and
               fastq1.defline.readNum != "1" and
               fastq2.defline.readNum == "1" ) or
             ( fastq1.defline.poreRead and
               fastq1.defline.poreRead == "complement" and
               fastq2.defline.poreRead == "template" ) ):
            save = fastq1
            fastq1 = fastq2
            fastq2 = save
        self.processFastqSpot ( ( fastq1, fastq2 ) ) # processFastqSpot expects a tuple
#        self.processPairFastqSpot ( fastq1, fastq2 )
#        self.dst['READ_TYPE']['data'] = array.array( 'B', [ 1, 1 ] )
        self.setDstName ( fastq1, self.dst )
        self.setDstSpotGroup ( fastq1, self.dst )
        self.gw.write(self.dst)

    ############################################################
    # Write spot where read1 in fastq and read2 in hash
    ############################################################

    def writeMixedRead1 ( self, fastq1 ):
        read2 = self.pairedRead2 [ fastq1.defline.name ]
        del self.pairedRead2 [ fastq1.defline.name ]
        self.dst[self.readColumn]['data'] = (fastq1.seq + read2['seq']).encode('ascii')
        if self.isColorSpace:
            self.dst['CS_KEY']['data'] = (fastq1.csKey + read2['csKey']).encode('ascii')
        self.dst['READ_START']['data'] = array.array( 'I', [ 0, len(fastq1.seq) ] )
        self.dst['READ_LENGTH']['data'] = array.array( 'I', [ len(fastq1.seq), len(read2['seq']) ] )
        if ( fastq1.defline.filterRead or
             read2['filterRead'] ):
            self.dst['READ_FILTER']['data'] = array.array('B', [ 1, 1 ] )
        else:
            self.dst['READ_FILTER']['data'] = array.array('B', [ 0, 0 ] )

        if ( self.isNumQual and
             fastq1.qual and
             read2['qual'] ):
            fastq1.qual += " "
        self.setDstQual ( fastq1.qual + read2['qual'], self.dst )

        self.dst['READ_TYPE']['data'] = array.array( 'B', [ 1, 1 ] )
        self.setDstName ( fastq1, self.dst )
        self.setDstSpotGroup ( fastq1, self.dst )
        self.gw.write(self.dst)

    ############################################################
    # Write spot where read2 in fastq and read1 in hash
    ############################################################

    def writeMixedRead2 ( self, fastq2 ):
        read1 = self.pairedRead1 [ fastq2.defline.name ]
        del self.pairedRead1 [ fastq2.defline.name ]
        self.dst[self.readColumn]['data'] = (read1['seq'] + fastq2.seq ).encode('ascii')
        if self.isColorSpace:
            self.dst['CS_KEY']['data'] = ( read1['csKey'] + fastq2.csKey ).encode('ascii')
        self.dst['READ_START']['data'] = array.array( 'I', [ 0, len(read1['seq']) ] )
        self.dst['READ_LENGTH']['data'] = array.array( 'I', [ len(read1['seq']), len(fastq2.seq) ] )
        if ( read1['filterRead'] or
             fastq2.defline.filterRead ):
            self.dst['READ_FILTER']['data'] = array.array('B', [ 1, 1 ] )
        else:
            self.dst['READ_FILTER']['data'] = array.array('B', [ 0, 0 ] )

        if ( self.isNumQual and
             read1['qual'] and
             fastq2.qual ):
            read1['qual'] += " "

        self.setDstQual ( read1['qual'] + fastq2.qual, self.dst )

        self.dst['READ_TYPE']['data'] = array.array( 'B', [ 1, 1 ] )

        # qiime name can vary between read1 and read2 (read1 name is used)

        if fastq2.defline.qiimeName:
            fastq2.defline.qiimeName = read1['qiimeName']

        if fastq2.defline.suffix:
            fastq2.defline.suffix = read1['suffix']

        self.setDstName ( fastq2, self.dst )
        self.setDstSpotGroup ( fastq2, self.dst )
        self.gw.write(self.dst)

    ############################################################
    # Write spot where read in fastq is an orphan
    ############################################################

    def writeMixedOrphan ( self, fastq ):
        self.processFastqSpot ( ( fastq, ) ) # processFastqSpot expects a tuple
#        self.dst['READ_TYPE']['data'] = array.array( 'B', [ 1 ] )
        self.setDstName ( fastq, self.dst )
        self.setDstSpotGroup ( fastq, self.dst )
        self.gw.write(self.dst)

    ############################################################
    # Populate paired read hashes
    ############################################################

    def addToPairedReadHash ( self, fastq ):
        if self.isColorSpace:
            readValues = { 'seq':fastq.seq, 'qual':fastq.qual, 'filterRead':fastq.defline.filterRead, 'csKey':fastq.csKey }
        else:
            readValues = { 'seq':fastq.seq, 'qual':fastq.qual, 'filterRead':fastq.defline.filterRead }

        # qiime name can vary between read1 and read2

        if fastq.defline.qiimeName:
            readValues['qiimeName'] = fastq.defline.qiimeName

        if fastq.defline.suffix:
            readValues['suffix'] = fastq.defline.suffix
        else:
            readValues['suffix'] = None

        if ( ( fastq.defline.readNum and
               fastq.defline.readNum == "1") or
             ( fastq.defline.poreRead and
               fastq.defline.poreRead == "template") or
             ( fastq.defline.tagType and
               fastq.defline.tagType == "F3") or
             ( fastq.defline.deflineType == fastq.defline.ABSOLID and
               fastq.defline.tagType == '' ) ):
            self.pairedRead1[fastq.defline.name] = readValues
        else:
            self.pairedRead2[fastq.defline.name] = readValues

    ############################################################
    # Populate paired read hashes
    ############################################################

    def addToPoreReadHash ( self, fastq ):
        if ( fastq and
             not fastq.eof):
            if fastq.defline.poreRead != "2D":
                self.addToPairedReadHash ( fastq )
            elif not fastq.defline.name in self.nanopore2Dread :
                readValues = { 'seq':fastq.seq, 'qual':fastq.qual, 'filterRead':fastq.defline.filterRead, 'suffix':fastq.defline.suffix }
                self.nanopore2Dread[fastq.defline.name] = readValues

    ############################################################
    # Dump leftover orphan reads
    # Not sure if it would be quicker to del entries in each hash
    # when found
    ############################################################

    def processOrphanReads ( self, fastq1 ):
        if ( fastq1.defline.name in self.pairedRead1 or
             fastq1.defline.name in self.pairedRead2 ):
            self.writeMixedOrphan ( fastq1 )
        if fastq1.spotCount % self.infoSizeForSpotCount == 0:
            self.outputCount(True,fastq1)

    ############################################################
    # Set name to be written in archive
    # Keeping qiimeName separate because readNum may be present
    # in the qiimeName and I have no way to remove it. But I can
    # remove it from what occurs after the qiimeName and use that
    # for comparing read names.
    ############################################################

    def setDstName ( self, fastq, dst ):
        self.nameCount += 1
        if ( self.ignAndDiscardNames or
             self.useAndDiscardNames ):
            pass
        else:
            if self.fabricateNames: # Used only for AB Solid submissions where names are repeated/identical/too many
                name = str(fastq.spotCount)
            else:
                name = fastq.defline.name
            if ( self.ignoreNames and
                 self.firstLabelRE ):
                m = re.search( self.firstLabelRE, name )
                if m:
                    name = name[:m.start()]
            if fastq.defline.qiimeName :
                name = fastq.defline.qiimeName + "_" + name
            if fastq.defline.ignoredLeadChars :
                name = fastq.defline.ignoredLeadChars + name
            if fastq.defline.suffix and (not self.appendPoreReadToName):
                name += fastq.defline.suffix
##            if fastq.defline.ignoredTailChars :
##                name = name + fastq.defline.ignoredTailChars
            if self.appendBCtoName:
                spotGroup = fastq.defline.spotGroup
                if self.discardBCchars:
                    if self.discardBCchars < 0:
                        spotGroup = spotGroup[0:len(spotGroup) + self.discardBCchars]
                    else:
                        spotGroup = spotGroup[self.discardBCchars:]
                if spotGroup:
                    name = name + "_" + spotGroup
            dst['NAME']['data'] = name.encode('ascii')

    ############################################################
    # Set spot group in destination array
    # The bar code may not provide the final spot group that is
    # indicated in the run xml data block.
    ############################################################

    def setDstSpotGroup ( self, fastqs, dst ):
        if self.appendBCtoName:
            return
        if not ( type ( fastqs ) is tuple ):
            fastqs = ( fastqs, )
        spotGroup = ''
        if not self.discardBarcodes:
            if self.spotGroupProvided:
                spotGroup = self.spotGroup
            elif ( self.typesProvided and
                   self.readCount != self.readCountNonGroup ):
                start = 0
                readIndex = -1
                for readType in self.readTypes:
                    readIndex += 1
                    if readType == self.READ_TYPE_GROUP:
                        if spotGroup:
                            spotGroup += "_"

                        if len(fastqs) > 1:
                            spotGroup += fastqs[readIndex].seq
                        else:
                            spotGroup += fastqs[0].seq[ start : start + self.readLengths[readIndex] ]
                            start += self.readLengths[readIndex]

            elif ( fastqs[0].defline.spotGroup and
                   not self.ignAndDiscardNames ):
                spotGroup = fastqs[0].defline.spotGroup
                if self.discardBCchars:
                    if self.discardBCchars < 0:
                        spotGroup = spotGroup[0:len(spotGroup) + self.discardBCchars]
                    else:
                        spotGroup = spotGroup[self.discardBCchars:]

            elif ( self.ignoreSGforPairing and  # For the 10x case where the Index read does not have a spot group but the R1 read does.
                   not fastqs[0].defline.spotGroup and
                   fastqs[1].defline.spotGroup and
                   not self.ignAndDiscardNames ):
                spotGroup = fastqs[1].defline.spotGroup
                if self.discardBCchars:
                    if self.discardBCchars < 0:
                        spotGroup = spotGroup[0:len(spotGroup) + self.discardBCchars]
                    else:
                        spotGroup = spotGroup[self.discardBCchars:]

        dst['SPOT_GROUP']['data'] = spotGroup.encode('ascii')
        if ( spotGroup and
             not self.spotGroupProvided and
             not spotGroup in self.spotGroupsFound ):
            self.spotGroupsFound[spotGroup] = 1

    ############################################################
    # Set quality in destination array
    ############################################################

    def setDstQual ( self, qualString, dst ):
        if self.isNumQual:
            qualVals = self.getNumQualArray ( qualString )
            dst['QUALITY']['data'] = qualVals
        else:
            if self.changeNegOneQual:
                qualString = qualString.translate(self.transNegOne)
            dst['QUALITY']['data'] = qualString.encode('ascii')

    ############################################################
    # Get numerical quality array
    ############################################################

    def getNumQualArray ( self, qualString ):

        qualValStrings = qualString.split()
        qualIndex = -1

        if self.logOdds:
            qualVals = array.array('b', [0] * len(qualValStrings) )
        else:
            qualVals = array.array('B', [0] * len(qualValStrings) )

        for qualValString in qualValStrings:
            qualIndex += 1
            qualVal = int(qualValString)

            # Handle case where dots or Ns are sometimes assigned '-1' for quality

            if ( self.changeNegOneQual and
                 qualVal == -1 ):
                qualVal = 0

            qualVals[qualIndex] = qualVal

        return qualVals

    ############################################################
    # Set quality offset
    ############################################################

    def setQualityOffset ( self, offset, provided ):

        if offset in ('PHRED_0', '0'):
            offset = 0
        elif offset in ('PHRED_33', '33'):
            offset = 33
        elif offset in ('PHRED_64', '64'):
            offset = 64

        if offset not in (0, 33, 64):
            self.statusWriter.outputErrorAndExit(f"Invalid offset determined or specified. Allowed values are 0 (numerical quality), 33, or 64 ... {offset}")
        else:
            self.offset = offset

            if self.offset == 64:
                if self.logOdds:
                    self.tbl['SEQUENCE']['QUALITY']['expression'] = '(INSDC:quality:text:log_odds_64)QUALITY'
                else:
                    self.tbl['SEQUENCE']['QUALITY']['expression'] = '(INSDC:quality:text:phred_64)QUALITY'

            elif self.offset == 0:
                if self.logOdds:
                    self.tbl['SEQUENCE']['QUALITY']['expression'] = '(INSDC:quality:log_odds)QUALITY'
                else:
                    self.tbl['SEQUENCE']['QUALITY']['expression'] = '(INSDC:quality:phred)QUALITY'
                self.isNumQual = True

            elif ( self.offset == 33 and
                   self.logOdds ):
                self.statusWriter.outputErrorAndExit( "Combining 33 offset with log odds is not yet supported" )

            if provided :
                self.offsetProvided = True

    ############################################################
    # Set read lengths from provided comma-separated list of lengths
    ############################################################

    def setReadLengths (self,readLengthString):
        readLengthString = readLengthString.strip()
        readLens = readLengthString.split(",")
        readLenCount = len(readLens)
        if self.readCount != 0:
            if readLenCount != self.readCount:
                self.statusWriter.outputErrorAndExit( "Read count disagreement between provided read lengths, {}, and prior argument read count, {}"
                                                      .format(str(readLenCount),str(self.readCount)))
        else:
            self.readCount = readLenCount

        readIndex = -1
        readStart = 0
        self.readLengths = array.array('I', [0] * self.readCount)
        self.readStarts = array.array('I', [0] * self.readCount)
        for readLen in readLens:
            readLen = readLen.strip()
            if readLen.isdigit():
                readIndex += 1
                self.readLengths[readIndex] = int(readLen)
                if self.readLengths[readIndex] == 0 :
                    self.variableLengthSpecified = True
                self.readStarts[readIndex] = readStart
                readStart += int(readLen)
            else:
                self.statusWriter.outputErrorAndExit( "Non-integer length specified ... {}".format(readLen))
        self.lengthsProvided = True

    ############################################################
    # Set read types from provided read type string
    ############################################################

    def setReadTypes (self,readTypeString):

        self.typesProvided = True

        readTypeString = readTypeString.strip()
        self.readTypeString = readTypeString
        transType = str.maketrans('', '', "BTG")
        empty = readTypeString.translate(transType)
        if len(empty) != 0:
            self.statusWriter.outputErrorAndExit( "Invalid read type specified (only B, T, or G allowed) ... {}".format(readTypeString) )

        else:
            typeCount = len(readTypeString)
            if self.readCount != 0:
                if typeCount != self.readCount:
                    self.statusWriter.outputErrorAndExit( "Read count disagreement between provided read types, {}, and prior argument read count, {}"
                                                          .format(typeCount,self.readCount))
            else:
                self.readCount = typeCount

            readIndex = -1
            self.readTypes = array.array('B', [0] * self.readCount)
            for readType in readTypeString:
                readIndex += 1
                if readType == "G":
                    self.readTypes[readIndex] = self.READ_TYPE_GROUP
                else:
                    self.readCountNonGroup += 1
                    if readType == "T":
                        self.readTypes[readIndex] = self.READ_TYPE_TECHNICAL
                    else:
                        self.readTypes[readIndex] = self.READ_TYPE_BIOLOGICAL

            if typeCount > 2:
                self.keepEmptySeqs = True

    ############################################################
    # Set read labels variables from provided read label string
    ############################################################

    def setReadLabels (self,readLabelString):
        readLabelString = readLabelString.strip()
        labels = readLabelString.split(",")
        labelCount = len(labels)
        if self.readCount != 0:
            if labelCount != self.readCount:
                self.statusWriter.outputErrorAndExit( "Read count disagreement between provided read labels, {}, and prior argument read count, {}"
                                                      .format(labelCount,self.readCount))
        else:
            self.readCount = labelCount

        labelStart = 0
        labelString = ''
        self.firstLabelRE = re.escape ( labels[0] )
        self.firstLabelRE += "$"
        self.labels = ''
        self.labelStarts = array.array('I', [] )
        self.labelLengths = array.array('I', [] )
        for label in labels:
            label = label.strip()
            if label != "": # Label is empty string for read designated for spot group (type 'G')
                labelString += label
                self.labelStarts.append(labelStart)

                label_len = len(label)
                self.labelLengths.append(label_len)
                labelStart += label_len
        self.labels = labelString.encode('ascii')
        self.labelsProvided = True

    ############################################################
    # Set spot group (need to add check for spaces/tabs)
    ############################################################

    def setSpotGroup (self,spotGroupString):
        self.spotGroup = spotGroupString.strip()
        self.spotGroupProvided = True

    ############################################################
    # Set platform
    ############################################################

    def setPlatform (self, platformString ):
        platformString.strip()
        platform = Platform.convertPlatformString(platformString)
        if platform is None :
            self.statusWriter.outputErrorAndExit( "Unrecognized platform was specified ... {}".format(platformString) )
        else:
            self.platform = array.array('B', [platform.value] )
            self.platformString = platformString
            if platform == Platform.SRA_PLATFORM_ILLUMINA:
                self.db = 'NCBI:SRA:Illumina:db'

    ############################################################
    # Ignore leading characters on deflines for pairing
    ############################################################

    def setIgnLeadChars ( self, ignLeadCharsNum ):
        self.ignLeadCharsNum = int(ignLeadCharsNum) - 1

    ############################################################
    # Ignore tailing characters on deflines for pairing
    ############################################################

    def setIgnTailChars ( self, ignTailCharsNum ):
        self.ignTailCharsNum = int(ignTailCharsNum)

    ############################################################
    # Ignore names and remove NAME column value
    ############################################################

    def setIgnoreNames ( self ):
        self.ignoreNames = True
        self.mixedDeflines = True

    ############################################################
    # Discard names and remove NAME column value
    ############################################################

    def setIgnAndDiscardNames ( self ):
        self.ignAndDiscardNames = True
        self.mixedDeflines = True
        self.db = 'NCBI:SRA:GenericFastqNoNames:db'
        del self.tbl['SEQUENCE']['NAME']

    ############################################################
    # Force load
    ############################################################

    def setForceLoad ( self ):
        self.forceLoad = True
        self.mixedDeflines = True

    ############################################################
    # Discard names and remove NAME column value
    ############################################################

    def setUseAndDiscardNames ( self ):
        self.useAndDiscardNames = True
        self.db = 'NCBI:SRA:GenericFastqNoNames:db'
        del self.tbl['SEQUENCE']['NAME']

    ############################################################
    # Discard barcodes because too many
    ############################################################

    def setDiscardBarcodes ( self ):
        self.discardBarcodes = True

    ############################################################
    # Discard leading barcode characters (not subtracting one
    # because the defline character does not come into consideration.
    ############################################################

    def setDiscardBCchars ( self, discardBCcharsNum ):
        self.discardBCchars = int(discardBCcharsNum)

    ############################################################
    # Set alternative schema to use
    ############################################################

    def setSchema ( self, schema ):
        self.schema = schema

    ############################################################
    # For special 10x related deflines containing BX: and BY:, set
    # these two options
    ############################################################

    def set10x ( self ):
        self.convertDeflineSpaces = True
        self.appendBCtoName = True

    ############################################################
    # Set name columns from provided comma-separated list of column numbers
    ############################################################

    def setNameColumns (self,nameColumnString):
        nameColumnString = nameColumnString.strip()
        if nameColumnString == '0':
            self.nameColumns = -1
        else:
            nameCols = nameColumnString.split(",")
            nameColumnCount = len(nameCols)
            colIndex = -1
            self.nameColumns = array.array('I', [0] * nameColumnCount)
            for nameCol in nameCols:
                nameCol = nameCol.strip()
                if nameCol.isdigit():
                    colIndex += 1
                    self.nameColumns[colIndex] = int(nameCol)
                else:
                    self.statusWriter.outputErrorAndExit( "Non-integer name column index specified ... {}".format(nameCol))

    ############################################################
    # Set seq columns from provided comma-separated list of column numbers
    ############################################################

    def setSeqColumns (self,seqColumnString):
        seqColumnString = seqColumnString.strip()
        seqCols = seqColumnString.split(",")
        seqColumnCount = len(seqCols)
        colIndex = -1
        self.seqColumns = array.array('I', [0] * seqColumnCount)
        for seqCol in seqCols:
            seqCol = seqCol.strip()
            if seqCol.isdigit():
                colIndex += 1
                self.seqColumns[colIndex] = int(seqCol)
            else:
                self.statusWriter.outputErrorAndExit( "Non-integer seq column index specified ... {}".format(seqCol))

    ############################################################
    # Set qual columns from provided comma-separated list of column numbers
    ############################################################

    def setQualColumns (self,qualColumnString):
        qualColumnString = qualColumnString.strip()
        qualCols = qualColumnString.split(",")
        qualColumnCount = len(qualCols)
        colIndex = -1
        self.qualColumns = array.array('I', [0] * qualColumnCount)
        for qualCol in qualCols:
            qualCol = qualCol.strip()
            if qualCol.isdigit():
                colIndex += 1
                self.qualColumns[colIndex] = int(qualCol)
            else:
                self.statusWriter.outputErrorAndExit( "Non-integer qual column index specified ... {}".format(qualCol))

    ############################################################
    # Set group columns from provided comma-separated list of column numbers
    ############################################################

    def setGroupColumns (self,groupColumnString):
        groupColumnString = groupColumnString.strip()
        groupCols = groupColumnString.split(",")
        groupColumnCount = len(groupCols)
        colIndex = -1
        self.groupColumns = array.array('I', [0] * groupColumnCount)
        for groupCol in groupCols:
            groupCol = groupCol.strip()
            if groupCol.isdigit():
                colIndex += 1
                self.groupColumns[colIndex] = int(groupCol)
            else:
                self.statusWriter.outputErrorAndExit( "Non-integer group column index specified ... {}".format(groupCol))

    ############################################################
    # Remove clip columns if not being used
    ############################################################

    def removeClipColumns ( self ):
        del self.tbl['SEQUENCE']['CLIP_QUALITY_LEFT']
        del self.tbl['SEQUENCE']['CLIP_QUALITY_RIGHT']

    ############################################################
    # Remove label columns if not being used
    ############################################################

    def removeLabelColumns ( self ):
        del self.tbl['SEQUENCE']['LABEL']
        del self.tbl['SEQUENCE']['LABEL_START']
        del self.tbl['SEQUENCE']['LABEL_LEN']

    ############################################################
    # Remove nanopore columns if not being used
    ############################################################

    def removeNanoporeColumns ( self ):
        del self.tbl['SEQUENCE']['READ_NO']
        del self.tbl['SEQUENCE']['CHANNEL']

    ############################################################
    # Remove absolid columns if not being used
    ############################################################

    def removeColorSpaceColumns ( self ):
        del self.tbl['SEQUENCE']['CSREAD']
        del self.tbl['SEQUENCE']['CS_KEY']

    ############################################################
    # Remove non-absolid columns if not being used
    ############################################################

    def removeNonColorSpaceColumns ( self ):
        del self.tbl['SEQUENCE']['READ']

    ############################################################
    # Determin if sharq can be executed
    ############################################################
    def isSharqCompatible(self):
        if ( self.useSharq and
             self.fileType=='fastq' and
             not self.mixedTypes and
             not self.mixedDeflines and
             not self.ignAndDiscardNames and
             not self.ignoreNames and
             not self.discardBarcodes and
             self.discardBCchars == 0 and
             not self.convertDeflineSpaces and
             not self.appendBCtoName and
             not self.seqColumns and
             not self.ignTailCharsNum and
             not self.ignLeadCharsNum and
             not self.spotGroupProvided and
             not self.labelsProvided and
             not self.lengthsProvided and
             not self.offsetProvided and
             not self.setClips and
             not self.duplicateReads and
             not self.logOdds and
             not self.ignCharsAfterPlus and
             not self.fabricateNames and
             not self.removeLastChar and
             not self.removeLineQuotes and
             not self.fastqWithoutAtSign and
             not self.fastqWithGtrThan and
             not self.ignoreSGforPairing and
             not self.noPairFileCheck and
             not self.seqOnlyFile and
             not self.genbankFile and
             not self.spotNumAtNameStart and
             not self.addFileToName and
             not self.convertDeflineSpaces and
             not self.convertEmptyDeflines and
             not self.badFlatNumQual and
             not self.retainAltNum and
             not self.genericDeflines and
             not self.removeSeqSpaces and
             not self.removeBadChars and
             not self.concatPairFiles and
             not self.bcSpaceOffset and
             not self.useFilenameForSG and
             not self.filesToIgnore and
             not self.filesToTruncate and
             not self.spacesInFilenames and
             not self.dashFileNames):
            return True
        else:
            return False

    ############################################################
    # Set nanopore columns READ_NUMBER and CHANNEL
    ############################################################

    @staticmethod
    def setNanoporeColumns ( fastq, dst ):
        dst['READ_NO']['data'] = array.array( 'I', [ int(fastq.defline.readNo) ] )
        dst['CHANNEL']['data'] = array.array( 'I', [ int(fastq.defline.channel) ] )

    ############################################################
    # Outputs spot count (need to clean thie up!)
    ############################################################

    def outputCount ( self, readCount, fastq, fastq2=None, fastq3=None, fastq4=None, fastq5=None, fastq6=None ):

        fileObj = fileObjects[fastq.filename]
        bytesRead = fileObj.tell()
        bytesRead2 = 0
        percent2 = None

        if fastq2:
            fileObj2 = fileObjects[fastq2.filename]
            bytesRead2 = fileObj2.tell()
        else:
            fileObj2 = None

        # Case where orphan file is larger than total size

        if self.dumpOrphans:
            size = fileSizes[fastq.filename]
            percent = ( int ( bytesRead  / size * 1000 ) ) / 10

        # Total size is sum of 1st files in each pair file set

        else:
            percent = ( int ( ( bytesRead + self.cumulativeBytesRead ) / self.totalSize * 1000 ) ) / 10
            if fileObj2:
                percent2 = ( int ( ( bytesRead2 + self.cumulativeBytesRead2 ) / self.totalSize2 * 1000 ) ) / 10

        if ( percent2 and
             percent2 < percent ):
            percent = percent2

        percentInt = int(percent)

        if ( fastq6 and fastq5 and fastq4 and fastq3 and fastq2 and
             fastq.filename != fastq6.filename ) :
            if ( fastq.spotCount == fastq2.spotCount and
                 fastq.spotCount == fastq3.spotCount and
                 fastq.spotCount == fastq4.spotCount and
                 fastq.spotCount == fastq5.spotCount and
                 fastq.spotCount == fastq6.spotCount ):
                if readCount:
                    self.statusWriter.outputInfo("Reading spot number {} from {}, {}, {}, {}, {}, and {} ({}%)"
                                                 .format(fastq.spotCount,fastq.filename,fastq2.filename,fastq3.filename,fastq4.filename,fastq5.filename,fastq6.filename,percent), percentInt )
                else:
                    self.statusWriter.outputInfo("Wrote approx. spot number {} from {}, {}, {}, {}, {}, and {} ({}%)"
                                                 .format(fastq.spotCount,fastq.filename,fastq2.filename,fastq3.filename,fastq4.filename,fastq5.filename,fastq6.filename,percent), percentInt )
            else:
                if readCount:
                    self.statusWriter.outputInfo("Reading spot number {} from {}, spot number {} from {}, spot number {} from {}, spot number {} from {}, spot number {} from {}, and spot number {} from {} ({}%)"
                                                 .format(fastq.spotCount,fastq.filename,fastq2.spotCount,fastq2.filename,fastq3.spotCount,fastq3.filename,fastq4.spotCount,fastq4.filename,fastq5.spotCount,fastq5.filename,fastq6.spotCount,fastq6.filename,percent), percentInt )
                else:
                    self.statusWriter.outputInfo("Wrote approx. spot number {} from {}, approx. spot number {} from {}, approx. spot number {} from {}, approx. spot number {} from {}, spot number {} from {}, and approx. spot number {} from {} ({}%)"
                                                 .format(fastq.spotCount,fastq.filename,fastq2.spotCount,fastq2.filename,fastq3.spotCount,fastq3.filename,fastq4.spotCount,fastq4.filename,fastq5.spotCount,fastq5.filename,fastq6.spotCount,fastq6.filename,percent), percentInt )

        elif ( fastq5 and fastq4 and fastq3 and fastq2 and
             fastq.filename != fastq5.filename ) :
            if ( fastq.spotCount == fastq2.spotCount and
                 fastq.spotCount == fastq3.spotCount and
                 fastq.spotCount == fastq4.spotCount and
                 fastq.spotCount == fastq5.spotCount ):
                if readCount:
                    self.statusWriter.outputInfo("Reading spot number {} from {}, {}, {}, {}, and {} ({}%)"
                                                 .format(fastq.spotCount,fastq.filename,fastq2.filename,fastq3.filename,fastq4.filename,fastq5.filename,percent), percentInt )
                else:
                    self.statusWriter.outputInfo("Wrote approx. spot number {} from {}, {}, {}, {}, and {} ({}%)"
                                                 .format(fastq.spotCount,fastq.filename,fastq2.filename,fastq3.filename,fastq4.filename,fastq5.filename,percent), percentInt )
            else:
                if readCount:
                    self.statusWriter.outputInfo("Reading spot number {} from {}, spot number {} from {}, spot number {} from {}, spot number {} from {}, and spot number {} from {} ({}%)"
                                                 .format(fastq.spotCount,fastq.filename,fastq2.spotCount,fastq2.filename,fastq3.spotCount,fastq3.filename,fastq4.spotCount,fastq4.filename,fastq5.spotCount,fastq5.filename,percent), percentInt )
                else:
                    self.statusWriter.outputInfo("Wrote approx. spot number {} from {}, approx. spot number {} from {}, approx. spot number {} from {}, approx. spot number {} from {}, and approx. spot number {} from {} ({}%)"
                                                 .format(fastq.spotCount,fastq.filename,fastq2.spotCount,fastq2.filename,fastq3.spotCount,fastq3.filename,fastq4.spotCount,fastq4.filename,fastq5.spotCount,fastq5.filename,percent), percentInt )

        elif ( fastq4 and fastq3 and fastq2 and
             fastq.filename != fastq4.filename ) :
            if ( fastq.spotCount == fastq2.spotCount and
                 fastq.spotCount == fastq3.spotCount and
                 fastq.spotCount == fastq4.spotCount ):
                if readCount:
                    self.statusWriter.outputInfo("Reading spot number {} from {}, {}, {}, and {} ({}%)"
                                                 .format(fastq.spotCount,fastq.filename,fastq2.filename,fastq3.filename,fastq4.filename,percent), percentInt )
                else:
                    self.statusWriter.outputInfo("Wrote approx. spot number {} from {}, {}, {}, and {} ({}%)"
                                                 .format(fastq.spotCount,fastq.filename,fastq2.filename,fastq3.filename,fastq4.filename,percent), percentInt )
            else:
                if readCount:
                    self.statusWriter.outputInfo("Reading spot number {} from {}, spot number {} from {}, spot number {} from {}, and spot number {} from {} ({}%)"
                                                 .format(fastq.spotCount,fastq.filename,fastq2.spotCount,fastq2.filename,fastq3.spotCount,fastq3.filename,fastq4.spotCount,fastq4.filename,percent), percentInt )
                else:
                    self.statusWriter.outputInfo("Wrote approx. spot number {} from {}, approx. spot number {} from {}, approx. spot number {} from {}, and approx. spot number {} from {} ({}%)"
                                                 .format(fastq.spotCount,fastq.filename,fastq2.spotCount,fastq2.filename,fastq3.spotCount,fastq3.filename,fastq4.spotCount,fastq4.filename,percent), percentInt )

        elif ( fastq3 and fastq2 and
               fastq.filename != fastq3.filename ) :
            if ( fastq.spotCount == fastq2.spotCount and
                 fastq.spotCount == fastq3.spotCount ):
                if readCount:
                    self.statusWriter.outputInfo("Reading spot number {} from {}, {}, and {} ({}%)"
                                                 .format(fastq.spotCount,fastq.filename,fastq2.filename,fastq3.filename,percent), percentInt )
                else:
                    self.statusWriter.outputInfo("Wrote approx. spot number {} from {}, {}, and {} ({}%)"
                                                 .format(fastq.spotCount,fastq.filename,fastq2.filename,fastq3.filename,percent), percentInt )
            else:
                if readCount:
                    self.statusWriter.outputInfo("Reading spot number {} from {}, spot number {} from {}, and spot number {} from {} ({}%)"
                                                 .format(fastq.spotCount,fastq.filename,fastq2.spotCount,fastq2.filename,fastq3.spotCount,fastq3.filename,percent), percentInt )
                else:
                    self.statusWriter.outputInfo("Wrote approx. spot number {} from {}, approx. spot number {} from {}, and approx. spot number {} from {} ({}%)"
                                                 .format(fastq.spotCount,fastq.filename,fastq2.spotCount,fastq2.filename,fastq3.spotCount,fastq3.filename,percent), percentInt )

        elif ( fastq2 and
               fastq.filename != fastq2.filename ):
            if fastq.spotCount == fastq2.spotCount:
                if readCount:
                    self.statusWriter.outputInfo("Reading spot number {} from {} and {} ({}%)"
                                                 .format(fastq.spotCount,fastq.filename,fastq2.filename,percent), percentInt )
                else:
                    self.statusWriter.outputInfo("Wrote approx. spot number {} from {} and {} ({}%)"
                                                 .format(fastq.spotCount,fastq.filename,fastq2.filename,percent), percentInt )
            else:
                if readCount:
                    self.statusWriter.outputInfo("Reading spot number {} from {}, spot number {} from {} ({}%)"
                                                 .format(fastq.spotCount,fastq.filename,fastq2.spotCount,fastq2.filename,percent), percentInt )
                else:
                    self.statusWriter.outputInfo("Wrote approx. spot number {} from {}, approx. spot number {} from {} ({}%)"
                                                 .format(fastq.spotCount,fastq.filename,fastq2.spotCount,fastq2.filename,percent), percentInt )
        else:
            if readCount:
                self.statusWriter.outputInfo("Reading spot number {} from {} ({}%)"
                                             .format(fastq.spotCount,fastq.filename,percent), percentInt )
            else:
                self.statusWriter.outputInfo("Wrote approx. spot number {} from {} ({}%)"
                                             .format(fastq.spotCount,fastq.filename,percent), percentInt )

############################################################
# Process command line arguments
############################################################

def processArguments():

    dashFileName = re.compile("^-.*(fastq|fq|txt|qual|fasta|csfasta|fa|fna|out)(.gz|.bz2|)$")

    if len(sys.argv[1:]) == 0 :
        usage( None , 0 )

    append1=False
    append2=False
    append3=False
    append4=False
    append5=False
    append6=False
    global filesToTruncate

    for arg in sys.argv[1:]:
        if append1:
            sw.read1PairFiles += arg
            if sw.read1PairFiles.endswith("\\"):
                sw.read1PairFiles = sw.read1PairFiles[:-1] + " "
            else:
                append1 = False
        elif append2:
            sw.read2PairFiles += arg
            if sw.read2PairFiles.endswith("\\"):
                sw.read2PairFiles = sw.read2PairFiles[:-1] + " "
            else:
                append2 = False
        elif append3:
            sw.read3PairFiles += arg
            if sw.read3PairFiles.endswith("\\"):
                sw.read3PairFiles = sw.read3PairFiles[:-1] + " "
            else:
                append3 = False
        elif append4:
            sw.read4PairFiles += arg
            if sw.read4PairFiles.endswith("\\"):
                sw.read4PairFiles = sw.read4PairFiles[:-1] + " "
            else:
                append4 = False
        elif append5:
            sw.filesToIgnore += arg
            if sw.filesToIgnore.endswith("\\"):
                sw.filesToIgnore = sw.filesToIgnore[:-1] + " "
            else:
                append5 = False
        elif append6:
            sw.filesToTruncate += arg
            if sw.filesToTruncate.endswith("\\"):
                sw.filesToTruncate = sw.filesToTruncate[:-1] + " "
            else:
                append6 = False
        elif ( arg[0] == '-' and
             ( "=" in arg or
               not ( dashFileName.match(arg) ) ) ) :
            if arg[2:9] == 'offset=':
                sw.setQualityOffset ( arg[9:], True )
            elif arg[2:9] == 'output=':
                sw.outdir = arg[9:]
            elif arg[2:10] == 'quality=':
                sw.setQualityOffset ( arg[10:], True )
            elif arg[2:11] == 'readLens=':
                sw.setReadLengths ( arg[11:] )
            elif arg[2:12] == 'readTypes=':
                sw.setReadTypes ( arg[12:] )
            elif arg[2:12] == 'spotGroup=':
                sw.setSpotGroup ( arg[12:] )
            elif arg[2:] == 'orphanReads':
                sw.orphanReads = True
                sw.allowEarlyFileEnd = True
            elif arg[2:] == 'logOdds':
                sw.logOdds = True
            elif arg[2:] == 'ignoreNames':
                sw.setIgnoreNames()
            elif arg[2:] == 'ignAndDiscardNames':
                sw.setIgnAndDiscardNames()
            elif arg[2:] == 'fabricateNames':
                sw.fabricateNames = True
            elif arg[2:17] == 'read1PairFiles=':
                sw.read1PairFiles = arg[17:]
                if sw.read1PairFiles.endswith("\\"):
                    append1=True
                    sw.read1PairFiles = sw.read1PairFiles[:-1] + " "
                    sw.spacesInFilenames = True
            elif arg[2:17] == 'read2PairFiles=':
                sw.read2PairFiles = arg[17:]
                if sw.read2PairFiles.endswith("\\"):
                    append2=True
                    sw.read2PairFiles = sw.read2PairFiles[:-1] + " "
                    sw.spacesInFilenames = True
            elif arg[2:17] == 'read3PairFiles=':
                sw.read3PairFiles = arg[17:]
                if sw.read3PairFiles.endswith("\\"):
                    append3=True
                    sw.read3PairFiles = sw.read3PairFiles[:-1] + " "
                    sw.spacesInFilenames = True
            elif arg[2:17] == 'read4PairFiles=':
                sw.read4PairFiles = arg[17:]
                if sw.read4PairFiles.endswith("\\"):
                    append4=True
                    sw.read4PairFiles = sw.read4PairFiles[:-1] + " "
                    sw.spacesInFilenames = True
            elif arg[2:17] == 'read5PairFiles=':
                sw.read5PairFiles = arg[17:]
            elif arg[2:17] == 'read6PairFiles=':
                sw.read6PairFiles = arg[17:]
            elif arg[2:17] == 'read1QualFiles=':
                sw.read1QualFiles = arg[17:]
            elif arg[2:17] == 'read2QualFiles=':
                sw.read2QualFiles = arg[17:]
            elif arg[2:16] == 'filesToIgnore=':
                sw.filesToIgnore = arg[16:]
                if sw.filesToIgnore.endswith("\\"):
                    append5=True
                    sw.filesToIgnore = sw.filesToIgnore[:-1] + " "
            elif arg[2:18] == 'filesToTruncate=':
                sw.filesToTruncate = arg[18:]
                sw.allowEarlyFileEnd = True
                if sw.filesToTruncate.endswith("\\"):
                    append6=True
                    sw.filesToTruncate = sw.filesToTruncate[:-1] + " "
            elif arg[2:17] == 'fileSpotGroups=':
                sw.fileSpotGroupsProvided = arg[17:]
                sw.spotGroupProvided = True
            elif arg[2:] == 'useFilenameForSG':
                sw.useFilenameForSG = True
                sw.spotGroupProvided = True
            elif arg[2:] == 'profile':
                sw.profile = True
            elif arg[2:11] == 'platform=':
                sw.setPlatform ( arg[11:] )
            elif arg[2:13] == 'readLabels=':
                sw.setReadLabels ( arg[13:] )
            elif arg[2:16] == 'maxErrorCount=':
                sw.maxErrorCount = int ( arg[16:] )
            elif arg[2:17] == 'maxErrorOutput=':
                sw.maxErrorOutput = int ( arg[17:] )
            elif arg[2:17] == 'maxSearchCount=':
                sw.maxSearchCount = int ( arg[17:] )
                if ( sw.maxErrorCount < sw.maxSearchCount and
                     sw.maxErrorCount != 0 ):
                    sw.maxErrorCount = sw.maxSearchCount
            elif arg[2:18] == 'maxSeqLineCount=':
                sw.maxSeqLineCountAllowed = int ( arg[18:] )
            elif arg[2:16] == 'maxDeflineLen=':
                sw.maxDeflineLen = int ( arg[16:] )
            elif arg[2:] == 'mixedDeflines':
                sw.mixedDeflines = True
            elif arg[2:] == 'genericDeflines':
                sw.genericDeflines = True
                sw.mixedDeflines = True
            elif arg[2:] == 'mixedTypes':
                sw.mixedTypes = True
            elif arg[2:] == 'isMultiLine':
                sw.isMultiLine = True
            elif arg[2:9] == 'schema=':
                sw.setSchema ( arg[9:] )
            elif arg[2:10] == 'xml-log=':
                statusWriter.setXmlLog( arg[10:] )
            elif arg[2:15] == 'ignLeadChars=':
                sw.setIgnLeadChars( arg[15:] )
            elif arg[2:15] == 'ignTailChars=':
                sw.setIgnTailChars( arg[15:] )
            elif arg[2:] == 'ignCharsAfterPlus':
                sw.ignCharsAfterPlus = True
            elif arg[2:] == 'discardBarcodes':
                sw.setDiscardBarcodes()
            elif arg[2:17] == 'discardBCchars=':
                sw.setDiscardBCchars( arg[17:] )
            elif arg[2:16] == 'appendBCtoName':
                sw.appendBCtoName = True
            elif arg[2:15] == 'bcSpaceOffset':
                sw.bcSpaceOffset = True
            elif arg[2:17] == 'concatPairFiles':
                sw.concatPairFiles = True
            elif arg[2:16] == 'removeBadChars':
                sw.removeBadChars = True
            elif arg[2:17] == 'removeSeqSpaces':
                sw.removeSeqSpaces = True
            elif arg[2:22] == 'convertDeflineSpaces':
                sw.convertDeflineSpaces = True
            elif arg[2:20] == 'ignoreSGforPairing':
                sw.ignoreSGforPairing = True
            elif arg[2:15] == 'addFileToName':
                sw.addFileToName = True
            elif arg[2:] == 'useAndDiscardNames':
                sw.setUseAndDiscardNames()
            elif arg[2:] == 'badFlatNumQual':
                sw.badFlatNumQual = True
            elif arg[2:] == 'convertEmptyDeflines':
                sw.convertEmptyDeflines = True
                sw.mixedDeflines = True
            elif arg[2:14] == 'nameColumns=':
                sw.setNameColumns ( arg[14:] )
            elif arg[2:13] == 'seqColumns=':
                sw.setSeqColumns ( arg[13:] )
            elif arg[2:14] == 'qualColumns=':
                sw.setQualColumns ( arg[14:] )
            elif arg[2:15] == 'groupColumns=':
                sw.setGroupColumns ( arg[15:] )
            elif arg[2:20] == 'abiLastPrimerBase=':
                sw.abiLastPrimerBase = arg[20:]
            elif arg[2:] == 'spotNumAtNameStart':
                sw.spotNumAtNameStart = True
            elif arg[2:] == 'allowEarlyFileEnd':
                sw.allowEarlyFileEnd = True
            elif arg[2:] == 'seqOnlyFile':
                sw.seqOnlyFile = True
                sw.genericDeflines = True
            elif arg[2:] == 'genbankFile':
                sw.genbankFile = True
                sw.genericDeflines = True
                sw.multiLine = True
                sw.removeSeqSpaces = True
            elif arg[2:] == 'keepEmptySeqs':
                sw.keepEmptySeqs = True
            elif arg[2:] == 'noPairFileCheck':
                sw.noPairFileCheck = True
            elif arg[2:] == 'fastqWithGtrThan':
                sw.fastqWithGtrThan = True
            elif arg[2:] == 'fastqWithoutAtSign':
                sw.fastqWithoutAtSign = True
            elif arg[2:] == 'removeLastChar':
                sw.removeLastChar = True
            elif arg[1:2] == 'z=':
                statusWriter.setXmlLog( arg[2:] )
            elif arg[2:] == 'useSharq':
                sw.useSharq = True
            elif arg[2:] == 'doNotUseSharq':
                sw.useSharq = False
            elif ( arg[2:] == 'help' or
                   arg[1:] == 'h'):
                usage(None,0)
            elif ( arg[2:] == 'version' or
                   arg[1:] == 'V'):
                sys.stderr.write("\nfastq-load.py.{}\n\n".format(version))
                sys.exit(0)
            else:
                usage( "Unrecognized option ... " + arg, 1 )
        else:
            if os.path.isdir(arg):
                for fastq in sh.find(arg, "-iname", "*.fastq"):
                    fastq = fastq.strip()
                    if os.path.isfile(fastq):
                        filePaths[fastq] = fastq
                for fastq in  sh.find(arg, "-iname", "*.fastq.gz"):
                    fastq = fastq.strip()
                    if os.path.isfile(fastq):
                        filePaths[fastq] = fastq
                for fastq in  sh.find(arg, "-iname", "*.fq"):
                    fastq = fastq.strip()
                    if os.path.isfile(fastq):
                        filePaths[fastq] = fastq
                for fastq in  sh.find(arg, "-iname", "*.fq"):
                    fastq = fastq.strip()
                    if os.path.isfile(fastq):
                        filePaths[fastq] = fastq
            elif re.match("^[.]_",arg):
                statusWriter.outputInfo( "Ignoring file beginning with '._': {}".format(arg))
            else:
                if not sw.dashFileNames and dashFileName.match(arg):
                    sw.dashFileNames = True
                filePaths[os.path.basename(arg)] = arg

    # Ensure read lengths and types are consistent

    if ( sw.lengthsProvided and
         ( len(sw.readLengths) != len(sw.readTypes) or
           ( sw.labelsProvided and
             ( len(sw.readLengths) != len(sw.labelLengths) ) ) ) ):
        statusWriter.outputErrorAndExit ( "Read types, read lengths, and read labels (if provided)\nmust have the same number of values" )

    # Ensure not orphan reads and ignore names

    if ( ( sw.ignoreNames or
           sw.ignAndDiscardNames ) and
         sw.orphanReads ):
        statusWriter.outputErrorAndExit ( "Cannot specify ignoreNames|ignAndDiscardNames and orphanReads at the same time" )

    # Process files to ignore

    moreFilesToIgnore = {}
    if sw.filesToIgnore:
        for filename in sw.filesToIgnore.split(","):
            if re.search("^\*",filename):
                suffix = filename[1:]
                for fileNameWithPath in filePaths:
                    if fileNameWithPath.endswith(suffix):
                        moreFilesToIgnore[fileNameWithPath] = 1
            elif re.search("\*$",filename):
                prefix = filename[:-1]
                for fileNameWithPath in filePaths:
                    if fileNameWithPath.startswith(prefix):
                        moreFilesToIgnore[fileNameWithPath] = 1
            else:
                if not filename in filePaths:
                    filename = getActualFileName(filename)
                del filePaths[filename]
    for filename in moreFilesToIgnore:
        del filePaths[filename]

    # Process files to truncate

    if sw.filesToTruncate:
        for filenameAndLine in sw.filesToTruncate.split(","):
            ( filename, lineNumber ) = filenameAndLine.split(":")
            if not filename in filePaths:
                filename = getActualFileName(filename)
            filesToTruncate[filename] = int(lineNumber)

    # Ensure fastq file paths were provided

    if len(filePaths) == 0 :
        statusWriter.outputErrorAndExit ( "No fastq files were specified" )

    # Setup logging and work directory

    logging.basicConfig(level=logging.DEBUG)

############################################################
# Get actual filename if filename in list is incorrect
############################################################
def getActualFileName (filename):
    extractSuffixMatch = re.compile("^(.*)(.gz|.bz2)$")
    if extractSuffixMatch.match(filename):
        m = extractSuffixMatch.match(filename)
        (filenameTrunc, suffix ) = m.groups()
        if filenameTrunc in filePaths:
            return filenameTrunc
        else:
            statusWriter.outputErrorAndExit( "Unrecognized file name in file list(1): {}".format(filename) )
    elif filename+'.gz' in filePaths:
        return filename+'.gz'
    elif filename+'.bz2' in filePaths:
        return filename+'.bz2'
    else:
        statusWriter.outputErrorAndExit( "Unrecognized file name in file list(2): {}".format(filename) )

############################################################
# Process provided file lists
############################################################

def processPairAndQualLists ():
    # Process read1, read2, read3, read4, read5, and read6 pair filename lists if provided

    sw.pairFilesProvided = True

    if ( sw.read6PairFiles and
         sw.read5PairFiles and
         sw.read4PairFiles and
         sw.read3PairFiles and
         sw.read2PairFiles and
         sw.read1PairFiles ):
        for (filename1,filename2,filename3,filename4,filename5,filename6) in zip(sw.read1PairFiles.split(","),sw.read2PairFiles.split(","),sw.read3PairFiles.split(","),sw.read4PairFiles.split(","),sw.read5PairFiles.split(","),sw.read6PairFiles.split(",")):
            if ( not filename1 or filename1 == "-" or
                 not filename2 or filename2 == "-" or
                 not filename3 or filename3 == "-" or
                 not filename4 or filename4 == "-" or
                 not filename5 or filename5 == "-" or
                 not filename6 or filename6 == "-" ):
                statusWriter.outputErrorAndExit( "6 fastq files must be provided when read6 files are indicated" )
            else:
                if not filename1 in fileHandles:
                    filename1 = getActualFileName(filename1)
                if not filename2 in fileHandles:
                    filename2 = getActualFileName(filename2)
                if not filename3 in fileHandles:
                    filename3 = getActualFileName(filename3)
                if not filename4 in fileHandles:
                    filename4 = getActualFileName(filename4)
                if not filename5 in fileHandles:
                    filename5 = getActualFileName(filename5)
                if not filename6 in fileHandles:
                    filename6 = getActualFileName(filename6)
                fileReadPairs[filename1] = filename2
                fileReadPairs2[filename1] = filename3
                fileReadPairs3[filename1] = filename4
                fileReadPairs4[filename1] = filename5
                fileReadPairs5[filename1] = filename6
                sw.keepEmptySeqs = True
                sw.useSharq = False

    elif ( sw.read5PairFiles and
           sw.read4PairFiles and
           sw.read3PairFiles and
           sw.read2PairFiles and
           sw.read1PairFiles ):
        for (filename1,filename2,filename3,filename4,filename5) in zip(sw.read1PairFiles.split(","),sw.read2PairFiles.split(","),sw.read3PairFiles.split(","),sw.read4PairFiles.split(","),sw.read5PairFiles.split(",")):
            if ( not filename1 or filename1 == "-" or
                 not filename2 or filename2 == "-" or
                 not filename3 or filename3 == "-" or
                 not filename4 or filename4 == "-" or
                 not filename5 or filename5 == "-" ):
                statusWriter.outputErrorAndExit( "5 fastq files must be provided when read5 files are indicated" )
            else:
                if not filename1 in fileHandles:
                    filename1 = getActualFileName(filename1)
                if not filename2 in fileHandles:
                    filename2 = getActualFileName(filename2)
                if not filename3 in fileHandles:
                    filename3 = getActualFileName(filename3)
                if not filename4 in fileHandles:
                    filename4 = getActualFileName(filename4)
                if not filename5 in fileHandles:
                    filename5 = getActualFileName(filename5)
                fileReadPairs[filename1] = filename2
                fileReadPairs2[filename1] = filename3
                fileReadPairs3[filename1] = filename4
                fileReadPairs4[filename1] = filename5
                sw.keepEmptySeqs = True
                sw.useSharq = False

    elif ( sw.read4PairFiles and
           sw.read3PairFiles and
           sw.read2PairFiles and
           sw.read1PairFiles ):
        for (filename1,filename2,filename3,filename4) in zip(sw.read1PairFiles.split(","),sw.read2PairFiles.split(","),sw.read3PairFiles.split(","),sw.read4PairFiles.split(",")):
            if ( not filename1 or filename1 == "-" or
                 not filename2 or filename2 == "-" or
                 not filename3 or filename3 == "-" or
                 not filename4 or filename4 == "-" ):
                statusWriter.outputErrorAndExit( "4 fastq files must be provided when read4 files are indicated" )
            else:
                if not filename1 in fileHandles:
                    filename1 = getActualFileName(filename1)
                if not filename2 in fileHandles:
                    filename2 = getActualFileName(filename2)
                if not filename3 in fileHandles:
                    filename3 = getActualFileName(filename3)
                if not filename4 in fileHandles:
                    filename4 = getActualFileName(filename4)
                fileReadPairs[filename1] = filename2
                fileReadPairs2[filename1] = filename3
                fileReadPairs3[filename1] = filename4
                sw.keepEmptySeqs = True

    elif ( sw.read3PairFiles and
           sw.read2PairFiles and
           sw.read1PairFiles ):
        for (filename1,filename2,filename3) in zip(sw.read1PairFiles.split(","),sw.read2PairFiles.split(","),sw.read3PairFiles.split(",")):
            if ( not filename1 or filename1 == "-" or
                 not filename2 or filename2 == "-" or
                 not filename3 or filename3 == "-" ):
                statusWriter.outputErrorAndExit ( "3 fastq files must be provided when read3 files are indicated" )
            else:
                if not filename1 in fileHandles:
                    filename1 = getActualFileName(filename1)
                if not filename2 in fileHandles:
                    filename2 = getActualFileName(filename2)
                if not filename3 in fileHandles:
                    filename3 = getActualFileName(filename3)
                fileReadPairs[filename1] = filename2
                fileReadPairs2[filename1] = filename3
                sw.keepEmptySeqs = True

    elif ( sw.read2PairFiles and
           sw.read1PairFiles ):

        for (filename1,filename2) in zip(sw.read1PairFiles.split(","),sw.read2PairFiles.split(",")):
            if ( not filename1 or
                 filename1 == "-"):
                statusWriter.outputInfo(f"Assuming this file is a fragment file {filename2}... ")
            elif ( not filename2 or
                   filename2 == "-" ):
                statusWriter.outputInfo(f"Assuming this file is a fragment file {filename1}... ")
            else:
                if not filename1 in fileHandles:
                    filename1 = getActualFileName(filename1)
                if not filename2 in fileHandles:
                    filename2 = getActualFileName(filename2)
                fileReadPairs[filename1] = filename2
                if sw.ignoreNames or sw.ignAndDiscardNames:
                    sw.keepEmptySeqs = True

##    Can't remember what this case was about
##
##    elif ( not sw.read1QualFiles and
##           not sw.read2QualFiles and
##           ( sw.read1PairFiles or
##             sw.read2PairFiles or
##             sw.read3PairFiles or
##             sw.read4PairFiles ) ):
##        statusWriter.outputErrorAndExit ( "Providing only read1, read2, read3 or read4 pair file lists is not valid" )

    # Process read 1 seq and qual filename lists if provided

    if ( sw.read1PairFiles and
         sw.read1QualFiles ):

        sw.qualFilesProvided = True
        for (filename1,filename2) in zip(sw.read1PairFiles.split(","),sw.read1QualFiles.split(",")):
            if ( not filename1 or
                 filename1 == "-" ):
                continue
            elif ( not filename2 or
                   filename2 == "-" ):
                statusWriter.outputErrorAndExit ( "Filename does not have a qual file even though read1QualFiles provided ... {}".format(filename1) )
            else:
                if not filename1 in fileHandles:
                    filename1 = getActualFileName(filename1)
                if not filename2 in fileHandles:
                    filename2 = getActualFileName(filename2)
                refineSeqQual (filename1, fileHandles[filename1], filename2, fileHandles[filename2] )

    # Process read 2 seq and qual filename lists if provided

    if ( sw.read2PairFiles and
         sw.read2QualFiles ):

        sw.qualFilesProvided = True
        for (filename1,filename2) in zip(sw.read2PairFiles.split(","),sw.read2QualFiles.split(",")):
            if ( not filename1 or
                 filename1 == "-" ):
                continue
            elif ( not filename2 or
                   filename2 == "-" ):
                statusWriter.outputErrorAndExit ( "Filename does not have a qual file even though read2QualFiles provided ... {}".format(filename1) )
            else:
                if not filename1 in fileHandles:
                    filename1 = getActualFileName(filename1)
                if not filename2 in fileHandles:
                    filename2 = getActualFileName(filename2)
                refineSeqQual (filename1, fileHandles[filename1], filename2, fileHandles[filename2] )

############################################################
# Determine if eight line fastq
############################################################

def isEightLineOrSameNameFastq (filename, handle, fileType=None):
    if fileType == "multiLine":
        fastq1 = getFastqInstance ( filename )
        fastq1.isMultiLine = True
        fastq2 = getFastqInstance ( filename )
        fastq2.isMultiLine = True
    elif fileType in ("seqQual", "fasta"):
        fastq1 = FastaFastqReader (filename, handle)
        fastq2 = FastaFastqReader (filename, handle)
    elif fileType in ("multiLineSeqQual", "multiLineFasta"):
        fastq1 = FastaFastqReader (filename, handle)
        fastq1.isMultiLine = True
        fastq2 = FastaFastqReader (filename, handle)
        fastq2.isMultiLine = True
    else:
        fastq1 = getFastqInstance ( filename )
        fastq2 = getFastqInstance ( filename )

    fastq1.isMultiLine = True # Just in case
    fastq2.isMultiLine = True # Just in case
    spotNames = {}
    spotCount = 0
    isEightLine = False
    fastq1.setStatus(sw)
    fastq2.setStatus(sw)
    foundLargeAltNum = False
    while True:
        spotCount += 1
        fastq1.read()
        fastq2.lineCount = fastq1.lineCount
        fastq2.lineCountQual = fastq1.lineCountQual
        if fastq1.defline.name:
            if fastq1.defline.name not in spotNames:
                spotNames[fastq1.defline.name] = 0
            spotNames[fastq1.defline.name] += 1
            if fastq1.isMultiLine:
                fastq2.savedDeflineStringSeq = fastq1.savedDeflineStringSeq
            fastq2.read()
            fastq1.lineCount = fastq2.lineCount
            fastq1.lineCountQual = fastq2.lineCountQual
            if fastq2.defline.name:
                if fastq2.defline.name not in spotNames:
                    spotNames[fastq2.defline.name] = 0
                spotNames[fastq2.defline.name] += 1
                if fastq2.isMultiLine:
                    fastq1.savedDeflineStringSeq = fastq2.savedDeflineStringSeq
        if ( not isEightLine and
             Defline.isPairedDeflines ( fastq1.defline, fastq2.defline, False, sw.ignoreSGforPairing, sw.discardBCchars ) ):
            if ( fastq1.defline.readNum and
                 fastq1.defline.readNum  == fastq2.defline.readNum and
                 fastq1.seq == fastq2.seq ):
                sw.duplicateReads = True
                if not fastq1.filename in fileWithDupReads:
                    statusWriter.outputInfo("Duplicate reads exist within {}".format(fastq1.filename) )
                    fileWithDupReads[fastq1.filename] = 1
            else:
                if (fastq1.defline.poreRead and # Look for 2D read, too
                    ( fastq1.defline.poreRead == "2D" or
                      fastq2.defline.poreRead == "2D" ) ):
                    filePore2D[filename] = filename
                    sw.orphanReads = True
                    sw.allowEarlyFileEnd = True
                isEightLine = True
        elif ( isEightLine and
               not filename in filePore2D and
               fastq1.defline.poreRead and  # Look for 2D read, too
               ( fastq1.defline.poreRead == "2D" or
                 fastq2.defline.poreRead == "2D" ) ):
            filePore2D[filename] = filename
            sw.orphanReads = True
            sw.allowEarlyFileEnd = True

        if ( not foundLargeAltNum and
             fastq1.defline.retainAltNum ):
            isEightLine = False
            foundLargeAltNum = True

        if (spotCount == 10000 or
            fastq1.eof):
            break

    for name in spotNames:
        if spotNames[name] > 6:
            isEightLine = False
            if fastq1.seqParser.isColorSpace:
                statusWriter.outputWarning("Ignoring spot names because file {} has repeated names (e.g. {}). Paired files must be specified."
                                           .format(filename,name) )
                sw.fabricateNames = True
                sw.isColorSpace = True
                sw.ignoreNames = True
            else:
                statusWriter.outputErrorAndExit("Please specify 'ignoreNames' or 'ignAndDiscardNames' or 'convertDeflineSpaces'. File {} has repeated names (e.g. {}). Paired files must be specified if ignoring names."
                                         .format(filename,name) )
            break

    return isEightLine

############################################################
# Check seqQual files for multiline, eightline
############################################################

def refineSeqQual ( filename, handle, filename2, handle2 ):
    fileSkip[filename2] = 1
    fileQualPairs[filename] = filename2
    fastq = SeqQualFastqReader (filename,handle,filename2,handle2)
    fastq.setStatus(sw)
    fastq.isMultiLine = True
    if fastq.isMultiLineFastq():
        sw.isMultiLine = True
        fileTypes[filename] = "multiLineSeqQual"
    else:
        fileTypes[filename] = "seqQual"

    # Check for '8 line' seq qual

    if ( sw.isEightLine or
         ( not ( sw.ignoreNames or sw.ignAndDiscardNames ) and
           isEightLineOrSameNameFastq ( filename, handle, fileTypes[filename] ) ) ):
        fileReadPairs[filename] = filename
        if fileTypes[filename] == "multiLineSeqQual":
            fileTypes[filename] = "multiLineEightLineSeqQual"
        else:
            fileTypes[filename] = "eightLineSeqQual"

    statusWriter.outputInfo("File {} is identified as {}".format(filename,fileTypes[filename] ) )
    statusWriter.outputInfo("File {} is identified as {} (seq file is {})".format(filename2,fileTypes[filename],filename) )

############################################################
# Check fasta files for multiline, eightline
############################################################

def refineFasta ( filename, handle ):
    if sw.genbankFile:
        fileTypes[filename] = "multiLineFasta"
        statusWriter.outputInfo("File {} is identified as a GenBank format file".format(filename) )
    else:
        fastq = FastaFastqReader (filename,handle)
        fastq.setStatus(sw)
        fastq.isMultiLine = True
        if fastq.isMultiLineFasta():
            sw.isMultiLine = True
            fileTypes[filename] = "multiLineFasta"
        else:
            fileTypes[filename] = "fasta"

        # Check for '8 line' fasta

        if ( not sw.seqOnlyFile and
             ( sw.isEightLine or
               ( not ( sw.ignoreNames or sw.ignAndDiscardNames ) and
                 isEightLineOrSameNameFastq ( filename, handle, fileTypes[filename] ) ) ) ):
            fileReadPairs[filename] = filename
            if fileTypes[filename] == "multiLineFasta":
                fileTypes[filename] = "multiLineEightLineFasta"
            else:
                fileTypes[filename] = "eightLineFasta"

        statusWriter.outputInfo("File {} is identified as {}".format(filename,fileTypes[filename]) )

############################################################
# Check fasta files for multiline, eightline
############################################################

def refineFastq ( filename, handle ):
    fastq = getFastqInstance ( filename )
    fastq.isMultiLine = True
    if fastq.isMultiLineFastq():
        sw.isMultiLine = True
        if ( not sw.concatPairFiles and
             ( sw.isEightLine or
               ( not ( sw.ignoreNames or sw.ignAndDiscardNames) and
                 isEightLineOrSameNameFastq ( filename, handle, "multiLine" ) ) ) ):
            fileReadPairs[filename] = filename
            fileTypes[filename] = "multiLineEightLine"
            sw.useSharq = False
        else:
            fileTypes[filename] = "multiLine"
    elif ( not sw.concatPairFiles and
             ( sw.isEightLine or
               ( not ( sw.ignoreNames or sw.ignAndDiscardNames) and
                 isEightLineOrSameNameFastq ( filename, handle ) ) ) ):
        fileReadPairs[filename] = filename
        fileTypes[filename] = "eightLine"
    else:
        fileTypes[filename] = "normal"

    statusWriter.outputInfo("File {} is identified as {}".format(filename, fileTypes[filename]) )

############################################################
# Determine file types
# Check for separate seq/qual files and single/eight line fastq
# Also check for fasta only file (i.e. no corresponding qual file)
# Check for multi-line seq/qual
############################################################

def setFileTypes ():
    startFileNum = 0
    skipLines = 0
    for filename in sorted(fileHandles):
        startFileNum += 1
        if ( filename in fileSkip or
             filename in fileQualPairs ):
            continue
        defline = Defline('')
        defline.setStatus ( sw )
        handle = fileHandles[filename]
        lineCount = 0
        FastqReader.processHeader(handle,defline,sw.seqOnlyFile,sw.removeLastChar)
        fileDeflineString = None
        if sw.seqOnlyFile:
            deflineChar = '>'
        else:
            if re.search ( "[._-]Index[._]" , filename ):
                sw.foundIndex = True
            elif re.search ( "[._-]I1[._]" , filename ):
                sw.foundI1 = True
            elif re.search ( "[._-]R1[._]" , filename ):
                sw.foundR1 = True
            elif re.search ( "[._-]R2[._]" , filename ):
                sw.foundR2 = True

            while True:
                fileDeflineString = handle.readline().strip()
                try:
                    fileDeflineString = fileDeflineString.decode("utf-8") # In case reading from bz2 file
                except:
                    pass
                lineCount += 1

                # Check for quoted lines (won't work for single line cases)

                if ( fileDeflineString.startswith('"') and
                     fileDeflineString.endswith('"') and
                     ( fileDeflineString[1:2] in "@>" ) ):
                    deflineChar = fileDeflineString[1:2]
                    sw.removeLineQuotes = True
                elif ( sw.fastqWithoutAtSign and
                       defline.parseDeflineString( fileDeflineString ) ):
                    deflineChar = "@"
                elif fileDeflineString:
                    deflineChar = fileDeflineString[0]
                else:
                    deflineChar = ''

                if ( sw.seqColumns or
                     sw.genbankFile or
                     ( fileDeflineString and re.match("LOCUS\s+(.+)",fileDeflineString) ) ):
                    break
                elif deflineChar not in '@>':
                    fastq = SingleLineFastqReader (filename,handle)
                    fastq.setStatus(sw)
                    fastq.outputStatus = False
                    if fastq.isSingleLineFastq(skipLines):
                        break
                    else:
                        skipLines = lineCount
                        if ( not sw.maxErrorOutput or
                             lineCount < sw.maxErrorOutput ):
                            statusWriter.outputInfo("Skipping unrecognized line {} at start of {}".format(lineCount,filename) )
                        fastq.restart()
                        for i in range(0,skipLines):
                            handle.readline()
                else:
                    break

                if ( ( sw.maxSearchCount != 0 and lineCount > sw.maxSearchCount ) or
                     lineCount > 100000 or
                     fastq.isEof() ):
                    statusWriter.outputErrorAndExit( "First {} lines (or all lines if EOF encountered) in this file were not recognized as a defline ... {}"
                                                     .format(sw.maxSearchCount,filename) )


        # Check for genbank format (treated like fasta)

        if ( sw.genbankFile or
             ( fileDeflineString and re.match("^LOCUS\s+(.+)",fileDeflineString) ) ):
            sw.genbankFile = True
            sw.genericDeflines = True
            sw.multiLine = True
            sw.removeSeqSpaces = True
            refineFasta ( filename, handle )
            sw.fileType = "fasta"

        # Check for 'single line' fastq

        elif not deflineChar or not (deflineChar in '@>') or sw.seqColumns:
            if sw.seqColumns:
                fileTypes[filename] = "singleLineColumns"
                if not sw.nameColumns:
                    sw.nameColumns = -1
                sw.fileType = "singleLine"
            else:
                fastq = SingleLineFastqReader (filename,handle)
                fastq.setStatus(sw)
                if fastq.isSingleLineFastq(skipLines):
                    if fastq.delim == "\t":
                        fileTypes[filename] = "singleLineQseq"
                        statusWriter.outputInfo("File {} is identified as single line fastq Illumina qseq or export".format(filename) )
                    else:
                        fileTypes[filename] = "singleLine"
                        statusWriter.outputInfo("File {} is identified as single line fastq".format(filename) )
                        sw.fileType = "singleLine"
                else:
                    statusWriter.outputErrorAndExit( "Defline without > or @ is not recognizable as single line fastq ... {}"
                                                     .format(fileDeflineString) )

        # Check for 'interleaved' or '8-line' fastq (8-line with ignored/discarded names must be user specified)
        # and/or 'multi-line' fastq. Check for 10x defline variation.

        elif ( deflineChar == "@" or
               ( sw.fastqWithGtrThan and deflineChar == '>') ):
            if ( fileDeflineString.find("BX:") != -1 and
                 fileDeflineString.find("BY:") != -1 ):
                sw.set10x()
            refineFastq (filename, handle)
            sw.fileType = "fastq"

        # If seqOnlyFile, then the file will be treated as containing a single fasta sequence

        elif sw.seqOnlyFile:
            refineFasta ( filename, handle )
            sw.fileType = "fasta"

        # Check for separate seq/qual files (deflineChar == ">")
        # Also check for fasta only file

        else:
            fileNum = 0
            for filename2 in sorted(fileHandles):
                fileNum += 1

                if ( fileNum <= startFileNum or
                     filename2 in fileQualPairs or
                     filename2 in fileSkip ):
                    continue
                else:

                    # Look for matching seq/qual file in remaining files

                    defline2 = Defline('')
                    handle2 = fileHandles[filename2]
                    FastqReader.processHeader(handle2,defline2,sw.seqOnlyFile,sw.removeLastChar)
                    file2DeflineString = handle2.readline().strip()
                    try:
                        file2DeflineString = file2DeflineString.decode("utf-8")
                    except:
                        pass

                    if ( file2DeflineString.startswith('"') and
                         file2DeflineString.endswith('"') and
                         ( file2DeflineString[1:2] in "@>" ) ):
                        deflineChar2 = file2DeflineString[1:2]
                    else:
                        deflineChar2 = file2DeflineString[0]

                    if deflineChar2 != ">":
                        handle2.seek(0)
                        continue
                    else:
                        fastq = SeqQualFastqReader (filename,handle,filename2,handle2)
                        fastq.setStatus(sw)
                        isSeqQualFastq = fastq.isSeqQualFastq()
                        handle2.seek(0)
                        if isSeqQualFastq == 1:
                            refineSeqQual(filename,handle,filename2,handle2)
                            sw.fileType = "split"
                            break
                        elif isSeqQualFastq == 2:
                            refineSeqQual(filename2,handle2,filename,handle)
                            sw.fileType = "split"
                            break

            if ( filename not in fileTypes and
                 filename not in fileSkip and
                 filename not in fileQualPairs ):
                fastq = FastaFastqReader (filename,handle)
                fastq.setStatus(sw)
                if fastq.isFastaFastq():
                    refineFasta ( filename, handle )
                    sw.fileType = "fasta"

                else:
                    fastq.restart()
                    fastq.read()
                    if not fastq.seq:
                        qual = fastq.prevLine
                    else:
                        qual = fastq.seq
                    if sw.badFlatNumQual:
                        qual = Qual.fixBadFlatNumQual(qual)
                    parsedQual = Qual(qual, 0)
                    if parsedQual.isValid :
                        statusWriter.outputErrorAndExit( "Unable to associate seq file with this qual file ... {}".format(filename) )
                    else:
                        statusWriter.outputErrorAndExit( "File with '>' defline character is unrecognized as containing sequence or quality ... {}"
                                                         .format(filename) )

        handle.seek(0)

############################################################
# Based on fileType, return fastq instance for provided filename
# Types are normal, singleLine, eightLine, seqQual, fasta.
############################################################

def getFastqInstance (filename):
    reader = None
    if filename in fileTypes:
        file_type = fileTypes[filename]
        if file_type in ("normal", "eightLine"):
            reader = FastqReader (filename,fileHandles[filename])
        elif file_type in ("multiLine", "multiLineEightLine"):
            reader = FastqReader (filename,fileHandles[filename])
            reader.isMultiLine = True
        elif file_type in ("seqQual", "eightLineSeqQual"):
            filenameQual = fileQualPairs[filename]
            reader = SeqQualFastqReader (filename,fileHandles[filename],filenameQual,fileHandles[filenameQual])
        elif file_type in ("multiLineSeqQual", "multiLineEightLineSeqQual"):
            filenameQual = fileQualPairs[filename]
            reader = SeqQualFastqReader (filename,fileHandles[filename],filenameQual,fileHandles[filenameQual])
            reader.isMultiLine = True
        elif file_type in ("fasta", "eightLineFasta"):
            reader = FastaFastqReader (filename,fileHandles[filename])
        elif file_type in ("multiLineFasta", "multiLineEightLineFasta"):
            reader = FastaFastqReader (filename,fileHandles[filename])
            reader.isMultiLine = True
        elif file_type == "singleLineQseq":
            reader = SingleLineFastqReader (filename,fileHandles[filename])
            reader.delim = "\t"
        elif file_type == "singleLineColumns":
            reader = SingleLineFastqReader (filename,fileHandles[filename])
            reader.delim = None
            reader.nameColumns = sw.nameColumns
            reader.seqColumns = sw.seqColumns
            reader.qualColumns = sw.qualColumns
            reader.groupColumns = sw.groupColumns
        elif file_type == "singleLine":
            reader = SingleLineFastqReader (filename,fileHandles[filename])
    else:
        reader = FastqReader (filename,fileHandles[filename])

    reader.setStatus(sw)

    if filename in filesToTruncate:
        reader.truncateLine = filesToTruncate[filename]

    return reader

############################################################
# Pair up files containing read pairs
############################################################

def setFilePairs ():

    if sw.ignoreNames or sw.ignAndDiscardNames or sw.noPairFileCheck:
        return

    if sw.typesProvided and not sw.read1PairFiles:
        sw.useFileOrder = 1
    elif ( ( sw.foundIndex or sw.foundI1 )
           and sw.foundR1 and sw.foundR2
           and (len(sw.readTypes) < 3) ):
        sw.useFileOrder = 1

    if sw.foundIndex and sw.foundR1 and sw.foundR2:
        sw.ignoreSGforPairing = True

    startFileNum = 0
    for filename1 in sorted(filePaths):
        if sw.concatPairFiles:
            fileReadPairs[filename1] = filename1
            fastq1 = getFastqInstance(filename1)
            fileConcatStart[filename1] = fastq1.findConcatSecReadStart()
            if fileConcatStart[filename1] == 0:
                statusWriter.outputErrorAndExit( "Unable to find concatenated file second read start ... {}".format(filename1) )
            else:
                statusWriter.outputInfo( "Concatenated file {} second read start line is {} at position {}"
                                         .format(filename1,fastq1.prevLineCount,fileConcatStart[filename1] ) )
            continue

        startFileNum += 1
        if ( filename1 in fileSkip or
            ( filename1 in fileReadPairs and not ( filename1 == fileReadPairs[filename1]))): # Excluding 8-line fastq to allow pairing with regular file
            continue
        fastq1 = getFastqInstance(filename1)
        fastq1.read()
        fileNum = 0

        # Search for corresponding pair file in the rest of the files

        for filename2 in sorted(fileHandles):
            fileNum += 1
            if ( fileNum <= startFileNum or
                (filename2 in fileReadPairs and not ( filename2 == fileReadPairs[filename2]) ) or # Excluding 8-line fastq to allow pairing with regular file
                 filename2 in fileSkip ):
                continue
            else:
                fastq2 = getFastqInstance(filename2)
                fastq2.read()

                isPairedDeflines = Defline.isPairedDeflines ( fastq1.defline, fastq2.defline, False, sw.ignoreSGforPairing, sw.discardBCchars, sw.useFileOrder )

                # Ensure filename with first read (if specified) is captured in fileReadPairs

                if isPairedDeflines:
                    processPairedDeflines ( filename1, filename2, fastq1.defline, fastq2.defline, isPairedDeflines )
                    if isPairedDeflines == 2: # In case more than 2 files are paired together
                        fastq1 = fastq2
                        filename1 = filename2

    # Check for pairs with orphan reads (if any unpaired files still exist)
    # Also checks for eight line fastq instances with orphans
    # Case count is dealing with the altNum situation (alternate spot or read numbering approaches)

    startFileNum = 0
    for filename1 in sorted(filePaths):
        startFileNum += 1
        if ( filename1 in fileSkip or
             filename1 in fileReadPairs or
             filename1 in filePore2D ):
            continue
        fastq1 = getFastqInstance(filename1)
        file1Reads = {}
        file1ReadNums = {}
        file1ReadSeqs = {}
        read1Count = 0
        caseCount = 0
        while True:
            fastq1.read()
            if fastq1.eof:
                break
            read1Count+=1
            if fastq1.defline.name in file1Reads:
                caseCount += 1
                if caseCount < 5:
                    continue
                elif ( fastq1.defline.readNum and
                       fastq1.defline.readNum == file1ReadNums[fastq1.defline.name] and
                       fastq1.seq == file1ReadSeqs[fastq1.defline.name] ):
                    sw.duplicateReads = True
                    if not fastq1.filename in fileWithDupReads:
                        statusWriter.outputInfo("Duplicate reads exist within {}".format(fastq1.filename) )
                        fileWithDupReads[fastq1.filename] = 1
                else:
                    sw.orphanReads = True
                    sw.allowEarlyFileEnd = True


                    # Possible to have concluded eightline but still checking for nanopore 2d

                    if not filename1 in fileReadPairs:
                        statusWriter.outputInfo("File {} is eight line fastq with orphan reads".format(filename1) )

                    fileReadPairs[filename1] = filename1

                    file_type = fileTypes[filename1]
                    if file_type == "normal":
                        fileTypes[filename1] = "eightLine"
                    elif file_type == "seqQual":
                        fileTypes[filename1] = "eightLineSeqQual"
                    elif file_type == "fasta":
                        fileTypes[filename1] = "eightLineFasta"
                    elif file_type == "multiLine":
                        fileTypes[filename1] = "multiLineEightLine"
                    elif file_type == "multiLineSeqQual":
                        fileTypes[filename1] = "multiLineEightLineSeqQual"
                    elif file_type == "multiLineFasta":
                        fileTypes[filename1] = "multiLineEightLineFasta"
                    elif file_type in ("eightLine", "eightLineSeqQual",
                                       "eightLineFasta", "multiLineEightLine",
                                       "multiLineEightLineSeqQual",
                                       "multiLineEightLineFasta"):
                        pass
                    else:
                        statusWriter.outputErrorAndExit ("Unsupported eight line fastq file with orphan reads ... {}".format(filename1) )

                    if fastq1.defline.poreRead:
                        defline1 = Defline(file1Reads[fastq1.defline.name])
                        if ( defline1.poreRead == "2D" or
                             fastq1.defline.poreRead == "2D" ):
                            filePore2D[filename1] = filename1

                    if ( filename1 in fileReadPairs and
                         ( not fastq1.defline.poreRead or
                           filename1 in filePore2D ) ):
                        break
            else:
                file1Reads[fastq1.defline.name] = fastq1.deflineStringSeq
                file1ReadNums[fastq1.defline.name] = fastq1.defline.readNum
                file1ReadSeqs[fastq1.defline.name] = fastq1.seq

            if read1Count == sw.maxPairCheckCount:
                break

        if fileTypes[filename1] in ("eightLine", "eightLineSeqQual",
                                    "eightLineFasta", "multiLineEightLine",
                                    "multiLineEightLineSeqQual",
                                    "multiLineEightLineFasta"):
            continue

        else:
            fileNum = 0

            # Search for corresponding pair file in the rest of the files (allowing for orphans)

            for filename2 in sorted(fileHandles):
                fileNum += 1
                if ( fileNum <= startFileNum or
                     filename2 in fileSkip ):
                    continue
                else:
                    read2Count = 0
                    fastq2 = getFastqInstance(filename2)
                    while True:
                        fastq2.read()
                        if fastq2.eof:
                            break
                        read2Count += 1
                        if fastq2.defline.name in file1Reads:

                            # Reset fastq1 defline object to state for retained fastq1.deflineStringSeq

                            fastq1.defline.parseDeflineString(file1Reads[fastq2.defline.name])

                            # Determine which file has lesser read (e.g. read1 vs read2)

                            isPairedDeflines = Defline.isPairedDeflines ( fastq1.defline, fastq2.defline, False, sw.ignoreSGforPairing, sw.discardBCchars )

                            # Use 'isPairedDeflines' values to retain lesser read

                            if isPairedDeflines:
                                sw.orphanReads = True
                                sw.allowEarlyFileEnd = True
                                statusWriter.outputInfo("Files {} and {} are paired with orphan reads".format(filename1,filename2) )
                                processPairedDeflines ( filename1, filename2, fastq1.defline, fastq2.defline, isPairedDeflines )
                                if ( ( filename1 in fileReadPairs and
                                       ( not fastq1.defline.poreRead or
                                         filename1 in filePore2D ) ) or
                                     filename2 in fileReadPairs or
                                     filename2 in filePore2D ):
                                    break
                        if read2Count == sw.maxPairCheckCount:
                            break
                if ( ( filename1 in fileReadPairs and
                       ( not fastq1.defline.poreRead or
                         filename1 in filePore2D ) ) or
                     filename2 in fileReadPairs or
                     filename2 in filePore2D ):
                    break

    # Check for nanopore 2D files still not paired (if any unpaired 2D files still exist)
    # Addressing case where template/complement reads in one file and 2D reads in another

    startFileNum = 0
    for filename1 in sorted(filePaths):
        startFileNum += 1
        if ( filename1 in fileSkip or
             filename1 in fileReadPairs or
             filename1 in fileReadPairs2 or
             filename1 in fileReadPairs3 or
             filename1 in filePore2D ):
            continue
        fastq1 = getFastqInstance(filename1)
        file1Reads = {}
        read1Count = 0
        poreRead2Dfound = False
        while True:
            fastq1.read()

            # Filename1 must be a non-paired nanopore file with poreRead value of '2D'

            if ( fastq1.eof or
                 not fastq1.defline.poreRead or
                 ( fastq1.defline.poreRead and
                   not fastq1.defline.poreRead == "2D" ) ):
                break
            read1Count+=1
            file1Reads[fastq1.defline.name] = fastq1.deflineStringSeq
            poreRead2Dfound = True
            if read1Count == sw.maxPairCheckCount:
                break

        # If filename1 is not a non-paired nanopore file with poreRead value of '2D' then move on

        if not poreRead2Dfound:
            continue

        # Filename1 is a non-paired nanopore file with a poreRead values of '2D', so try to find pair file

        else:
            for filename2 in sorted(fileHandles):
                if ( filename2 == filename1 or
                     filename2 in fileSkip or
                     filename2 in filePore2D ):
                    continue
                else:
                    read2Count = 0
                    fastq2 = getFastqInstance(filename2)
                    while True:
                        fastq2.read()
                        if fastq2.eof:
                            break
                        read2Count += 1
                        if fastq2.defline.name in file1Reads:
                            fastq1.defline.parseDeflineString ( file1Reads[fastq2.defline.name] )
                            isPairedDeflines = Defline.isPairedDeflines ( fastq1.defline, fastq2.defline, False, sw.ignoreSGforPairing, sw.discardBCchars )
                            if isPairedDeflines:
                                sw.orphanReads = True
                                sw.allowEarlyFileEnd = True
                                processPairedDeflines ( filename1, filename2, fastq1.defline, fastq2.defline, isPairedDeflines )
                                if filename2 in filePore2D:
                                    statusWriter.outputInfo("Files {} and {} are paired with orphan reads".format(filename2,filename1) )
                                    break
                        if read2Count == sw.maxPairCheckCount:
                            break
                if filename2 in filePore2D:
                    break

############################################################
def processPairedDeflines ( filename1, filename2, defline1, defline2, isPairedDeflines ):

    # Check for nanopore 2D

    if defline1.poreRead:

        if isPairedDeflines == 1:
            fileSkip[filename2] = 1
            if defline2.poreRead == "complement":
                fileReadPairs[filename1] = filename2
                sw.foundComplement = True
            elif defline2.poreRead == "2D":
                filePore2D[filename1] = filename2
                sw.found2D = True
            else:
                statusWriter.outputErrorAndExit("Unrecognized poreRead type2 ... {}\t{}\t{}\t{}\t{}"
                                                .format(defline2.poreRead,defline1.name,defline2.name,filename1,filename2) )
        else:
            fileSkip[filename1] = 1
            if defline1.poreRead == "complement":
                fileReadPairs[filename2] = filename1
                sw.foundComplement = True
            elif defline1.poreRead == "2D":
                filePore2D[filename2] = filename1
                sw.found2D = True
            else:
                statusWriter.outputErrorAndExit("Unrecognized poreRead type1 ... {}\t{}\t{}\t{}\t{}"
                                                .format(defline1.poreRead,defline2.name,defline1.name,filename1,filename2) )

    # Process paired files other than nanopore

    else:

        if isPairedDeflines == 1:
            if not filename1 in fileReadPairs:
                fileReadPairs[filename1] = filename2
            elif not filename1 in fileReadPairs2:
                fileReadPairs2[filename1] = filename2
            elif not filename1 in fileReadPairs3:
                fileReadPairs3[filename1] = filename2
            fileSkip[filename2] = 1
        else:
            if not filename2 in fileReadPairs:
                fileReadPairs[filename2] = filename1
            elif not filename2 in fileReadPairs2:
                fileReadPairs2[filename2] = filename1
            elif not filename2 in fileReadPairs3:
                fileReadPairs3[filename2] = filename1
            fileSkip[filename1] = 1

        # Check for consistent file types (may be an issue with multiLine type
        # paired with non-multiLine type)

        type1 = fileTypes[filename1]
        type2 = fileTypes[filename2]
        if type1 != type2:
            if ( not sw.mixedTypes and
                type1 != type2 ):
                statusWriter.outputErrorAndExit( "Paired files have different file types ... {} is {}, {} is {}"
                                                .format(filename1,type1,filename2,type2) )
            elif (  isPairedDeflines == 1 and
                    type1 == 'normal' and
                    type2 == 'eightLine'):
                del fileReadPairs[filename2]
                if not filename1 in fileReadPairs2:
                    fileReadPairs2[filename1] = filename2
                else:
                    fileReadPairs3[filename1] = filename2

############################################################
# Check for nanopore files that contain all three read types that were
# not captured by the eight-line fastq check or 2D only nanopore files
# Also check for color space files to handle them appropriately
############################################################

def checkForColorSpaceAndPlatform():

    abLabelSet = False

    for filename in sorted(filePaths):
        if filename in fileSkip:
            continue
        fastq = getFastqInstance(filename)
        fastq.read()
        parsedSeq = Seq(fastq.seqOrig) # fastq.seq has cs key removed so need seqOrig to detect color space
        parsedQual = Qual(fastq.qual,parsedSeq.length)

        # Check for presence of illumina reads

        if (fastq.defline.platform == "ILLUMINA" and
            sw.platformString == ''):
            sw.setPlatform("ILLUMINA")

        # Check for presence of pacbio reads

        elif (fastq.defline.platform == "PACBIO" and
            sw.platformString == ''):
            sw.setPlatform("PACBIO")

        # Check for presence of absolid reads

        elif ( fastq.defline.panel or
             parsedSeq.isColorSpace ):
            sw.setPlatform("ABSOLID")
            if parsedSeq.isColorSpace:
                sw.isColorSpace = True
                statusWriter.outputInfo( "Sequence values for file {} are color space".format(filename) )
                if parsedQual.length == parsedSeq.length + 1:
                    sw.extraQualForCSkey = True

            if fastq.defline.tagType:
                fileLabels[filename] = fastq.defline.tagType
                if fileTypes[filename] in ("eightLine", "eightLineSeqQual",
                                           "eightLineFasta", "multiLineEightLine",
                                           "multiLineEightLineSeqQual",
                                           "multiLineEightLineFasta"):
                    fastq.read()
                    fileLabels[filename] += "," + fastq.defline.tagType
                elif filename in fileReadPairs:
                    pairFilename = fileReadPairs[filename]
                    fastq2 = getFastqInstance(pairFilename)
                    fastq2.read()
                    fileLabels[filename] += "," + fastq2.defline.tagType
                    if fastq2.defline.tagType == "BC":
                        sw.readTypes[1] = sw.READ_TYPE_TECHNICAL
                if not abLabelSet:
                    sw.setReadLabels(fileLabels[filename])
                    abLabelSet = True

        # Check for presence of nanopore reads

        elif fastq.defline.platform == "NANOPORE":
            if sw.platformString:
                if sw.platformString != "NANOPORE":
                    statusWriter.outputErrorAndExit("Unable to mix nanopore reads with other platforms")
            else:
                sw.setPlatform("NANOPORE")

            if filename in filePore2D :
                pass
            else:
                fastq.restart()
                foundTemplate = False
                foundComplement = False
                found2D = False
                foundNoLabel = False
                readCount = 0
                while True:
                    fastq.read()
                    readCount += 1
                    if fastq.eof:
                        break
                    elif fastq.defline.poreRead == "template":
                        foundTemplate = True
                    elif fastq.defline.poreRead == "complement":
                        foundComplement = True
                    elif fastq.defline.poreRead == "2D":
                        found2D = True
                    else:
                        foundNoLabel = True

                    if (readCount == sw.maxPairCheckCount or
                        fastq.eof or
                        ( foundTemplate and foundComplement and found2D ) ):
                        break

                # Template and complement in the same file

                if ( foundTemplate and
                     foundComplement):
                    fileReadPairs[filename] = filename
                    sw.foundTemplate = True
                    sw.foundComplement = True

                # Template and 2D in the same file

                if ( foundTemplate and
                     found2D):
                    filePore2D[filename] = filename
                    sw.foundTemplate = True
                    sw.found2D = True

                # Template only in this file

                if ( foundTemplate and
                     not foundComplement and
                     not found2D and
                     not sw.orphanReads ):
                    sw.nanoporeTemplateOnly = True
                    sw.foundTemplate = True

                # Complement only in this file

                if ( foundComplement and
                     not foundTemplate and
                     not found2D ):
                    sw.foundComplement = True

                # Template and complement only in this file

                if ( foundComplement and
                     foundTemplate and
                     not found2D ):
                    sw.nanoporeTemplateComplementOnly = True
                    sw.foundTemplate = True
                    sw.foundComplement = True

                # Template and complement not found in this file

                if ( not ( foundTemplate or
                           foundComplement ) ):

                    # 2D found in this file

                    if ( found2D and
                         not sw.orphanReads ) :
                        sw.nanopore2Donly = True
                        sw.found2D = True

                # Check for presence of unlabeled reads (forces load of SEQUENCE table only)

                if foundNoLabel:
                    statusWriter.outputInfo("Nanopore reads without labels are present in {}".format(filename))
                    if filename in fileReadPairs:
                        del fileReadPairs[filename]
                    if filename in filePore2D:
                        del filePore2D[filename]
                    sw.nanoporeNoLabels = True

    # In case template, complement, and 2d are divided into separate files

    if sw.platformString == "NANOPORE":
        if ( ( sw.foundTemplate or
               sw.foundComplement ) and
             sw.nanopore2Donly ):
            sw.nanopore2Donly = False

        if ( ( sw.foundComplement or
               sw.found2D ) and
             sw.nanoporeTemplateOnly ):
            sw.nanoporeTemplateOnly = False

        if ( sw.found2D and
             sw.nanoporeTemplateComplementOnly ):
            sw.nanoporeTemplateComplementOnly = False

        if ( ( sw.found2D or sw.foundTemplate ) and
             sw.nanoporeNoLabels ):
            statusWriter.outputInfo("A mix of labeled and unlabeled nanopore reads are loaded at the same quality level")
            sw.appendPoreReadToName = True
        elif sw.nanopore2Donly:
            statusWriter.outputInfo("Only Nanopore 2D reads are present")
        elif sw.nanoporeTemplateOnly:
            statusWriter.outputInfo("Only Nanopore template reads are present")
        elif sw.nanoporeTemplateComplementOnly:
            statusWriter.outputInfo("Only Nanopore template and complement reads are present")

############################################################
# Update qual range for determination of min and max qual
# Not setting clip status yet
############################################################

def updateQualRangeAndClipStatus ( fastq1, fastq2, fastq3, fastq4, fastq5, fastq6 ):

    qual = Qual('',0)
    if ( fastq1.eof or
         ( fastq2 and
           fastq2.eof ) or
         ( fastq3 and
           fastq3.eof ) or
         ( fastq4 and
           fastq4.eof ) or
         ( fastq5 and
           fastq5.eof ) or
         ( fastq6 and
           fastq6.eof ) ):
        return

    if fastq6:
        if ( fastq1.qual and
             fastq2.qual and
             fastq3.qual and
             fastq4.qual and
             fastq5.qual and
             fastq6.qual and
             ( " " in fastq1.qual or
               ( fastq1.length == 1 and qual.isInt(fastq1.qual) ) ) and
             ( " " in fastq2.qual or
               ( fastq2.length == 1 and qual.isInt(fastq2.qual) ) ) and
             ( " " in fastq3.qual or
               ( fastq3.length == 1 and qual.isInt(fastq3.qual) ) ) and
             ( " " in fastq4.qual or
               ( fastq4.length == 1 and qual.isInt(fastq4.qual) ) ) and
             ( " " in fastq5.qual or
               ( fastq5.length == 1 and qual.isInt(fastq5.qual) ) ) and
             ( " " in fastq6.qual or
               ( fastq6.length == 1 and qual.isInt(fastq6.qual) ) ) and
             ( len(fastq1.qual) > 1 or
               len(fastq2.qual) > 1 or
               len(fastq3.qual) > 1 or
               len(fastq4.qual) > 1 or
               len(fastq5.qual) > 1 or
               len(fastq6.qual) > 1 ) ):
            combinedQual = fastq1.qual + " " + fastq2.qual + " " + fastq3.qual + " " + fastq4.qual + " " + fastq5.qual + " " + fastq6.qual
        else:
            combinedQual = fastq1.qual + fastq2.qual + fastq3.qual + fastq4.qual + fastq5.qual + fastq6.qual
        qual.parseQual( combinedQual, ( fastq1.length + fastq2.length + fastq3.length + fastq4.length + fastq5.length + fastq6.length ) ) # Providing seq length sum which may be different than length of combined qual

    elif fastq5:
        if ( fastq1.qual and
             fastq2.qual and
             fastq3.qual and
             fastq4.qual and
             fastq5.qual and
             ( " " in fastq1.qual or
               ( fastq1.length == 1 and qual.isInt(fastq1.qual) ) ) and
             ( " " in fastq2.qual or
               ( fastq2.length == 1 and qual.isInt(fastq2.qual) ) ) and
             ( " " in fastq3.qual or
               ( fastq3.length == 1 and qual.isInt(fastq3.qual) ) ) and
             ( " " in fastq4.qual or
               ( fastq4.length == 1 and qual.isInt(fastq4.qual) ) ) and
             ( " " in fastq5.qual or
               ( fastq5.length == 1 and qual.isInt(fastq5.qual) ) ) and
             ( len(fastq1.qual) > 1 or
               len(fastq2.qual) > 1 or
               len(fastq3.qual) > 1 or
               len(fastq4.qual) > 1 or
               len(fastq5.qual) > 1 ) ):
            combinedQual = fastq1.qual + " " + fastq2.qual + " " + fastq3.qual + " " + fastq4.qual + " " + fastq5.qual
        else:
            combinedQual = fastq1.qual + fastq2.qual + fastq3.qual + fastq4.qual + fastq5.qual
        qual.parseQual( combinedQual, ( fastq1.length + fastq2.length + fastq3.length + fastq4.length + fastq5.length ) ) # Providing seq length sum which may be different than length of combined qual

    elif fastq4:
        if ( fastq1.qual and
             fastq2.qual and
             fastq3.qual and
             fastq4.qual and
             ( " " in fastq1.qual or
               ( fastq1.length == 1 and qual.isInt(fastq1.qual) ) ) and
             ( " " in fastq2.qual or
               ( fastq2.length == 1 and qual.isInt(fastq2.qual) ) ) and
             ( " " in fastq3.qual or
               ( fastq3.length == 1 and qual.isInt(fastq3.qual) ) ) and
             ( " " in fastq4.qual or
               ( fastq4.length == 1 and qual.isInt(fastq4.qual) ) ) and
             ( len(fastq1.qual) > 1 or
               len(fastq2.qual) > 1 or
               len(fastq3.qual) > 1 or
               len(fastq4.qual) > 1 ) ):
            combinedQual = fastq1.qual + " " + fastq2.qual + " " + fastq3.qual + " " + fastq4.qual
        else:
            combinedQual = fastq1.qual + fastq2.qual + fastq3.qual + fastq4.qual
        qual.parseQual( combinedQual, ( fastq1.length + fastq2.length + fastq3.length + fastq4.length ) ) # Providing seq length sum which may be different than length of combined qual

    elif fastq3 and not fastq2:
        if ( fastq1.qual and
             fastq3.qual and
             ( " " in fastq1.qual or
               ( fastq1.length == 1 and qual.isInt(fastq1.qual) ) ) and
             ( " " in fastq3.qual or
               ( fastq3.length == 1 and qual.isInt(fastq3.qual) ) ) and
             ( len(fastq1.qual) > 1 or
               len(fastq3.qual) > 1 ) ):
            combinedQual = fastq1.qual + " " + fastq3.qual
        else:
            combinedQual = fastq1.qual + fastq3.qual
        qual.parseQual( combinedQual, ( fastq1.length + fastq3.length ) ) # Providing seq length sum which may be different than length of combined qual

    elif fastq3:
        if ( fastq1.qual and
             fastq2 and fastq2.qual and
             fastq3 and fastq3.qual and
             ( " " in fastq1.qual or
               ( fastq1.length == 1 and qual.isInt(fastq1.qual) ) ) and
             ( " " in fastq2.qual or
               ( fastq2.length == 1 and qual.isInt(fastq2.qual) ) ) and
             ( " " in fastq3.qual or
               ( fastq3.length == 1 and qual.isInt(fastq3.qual) ) ) and
             ( len(fastq1.qual) > 1 or
               len(fastq2.qual) > 1 or
               len(fastq3.qual) > 1 ) ):
            combinedQual = fastq1.qual + " " + fastq2.qual + " " + fastq3.qual
        else:
            combinedQual = fastq1.qual + fastq2.qual + fastq3.qual
        qual.parseQual( combinedQual, ( fastq1.length + fastq2.length + fastq3.length ) ) # Providing seq length sum which may be different than length of combined qual

    elif fastq2:
        if ( fastq1.qual and
             fastq2.qual and
             ( " " in fastq1.qual or
               ( fastq1.length == 1 and qual.isInt(fastq1.qual) ) ) and
             ( " " in fastq2.qual or
               ( fastq2.length == 1 and qual.isInt(fastq2.qual) ) ) and
             ( len(fastq1.qual) > 1 or
               len(fastq2.qual) > 1 ) ):
            combinedQual = fastq1.qual + " " + fastq2.qual
        else:
            combinedQual = fastq1.qual + fastq2.qual
        qual.parseQual( combinedQual, ( fastq1.length + fastq2.length ) ) # Providing seq length sum which may be different than length of combined qual
    else:
        qual.parseQual( fastq1.qual, fastq1.length )

    if qual.minQual < sw.minQual:
        sw.minQual = qual.minQual
    if qual.maxQual > sw.maxQual:
        sw.maxQual = qual.maxQual

    if ( sw.isNumQual and
         not qual.isNumQual and
         qual.length > 1 ):
        if fastq6:
            statusWriter.outputErrorAndExit(f"Possible mixed numerical and non-numerical quality  ... {fastq1.defline.name}, {fastq1.filename}, {fastq2.defline.name}, {fastq2.filename}, {fastq3.defline.name}, {fastq3.filename}, {fastq4.defline.name}, {fastq4.filename}, {fastq5.defline.name}, {fastq5.filename}, {fastq6.defline.name}, {fastq6.filename}")
        elif fastq5:
            statusWriter.outputErrorAndExit(f"Possible mixed numerical and non-numerical quality  ... {fastq1.defline.name}, {fastq1.filename}, {fastq2.defline.name}, {fastq2.filename}, {fastq3.defline.name}, {fastq3.filename}, {fastq4.defline.name}, {fastq4.filename}, {fastq5.defline.name}, {fastq5.filename}")
        elif fastq4:
            statusWriter.outputErrorAndExit(f"Possible mixed numerical and non-numerical quality  ... {fastq1.defline.name}, {fastq1.filename}, {fastq2.defline.name}, {fastq2.filename}, {fastq3.defline.name}, {fastq3.filename}, {fastq4.defline.name}, {fastq4.filename}")
        elif fastq3 and not fastq2:
            statusWriter.outputErrorAndExit(f"Possible mixed numerical and non-numerical quality  ... {fastq1.defline.name}, {fastq1.filename}, {fastq3.defline.name}, {fastq3.filename}")
        elif fastq3:
            statusWriter.outputErrorAndExit(f"Possible mixed numerical and non-numerical quality  ... {fastq1.defline.name}, {fastq1.filename}, {fastq2.defline.name}, {fastq2.filename}, {fastq3.defline.name}, {fastq3.filename}")
        elif fastq2:
            statusWriter.outputErrorAndExit(f"Possible mixed numerical and non-numerical quality  ... {fastq1.defline.name}, {fastq1.filename}, {fastq2.defline.name}, {fastq2.filename}")
        else:
            statusWriter.outputErrorAndExit(f"Possible mixed numerical and non-numerical quality  ... {fastq1.defline.name}, {fastq1.filename}")
    elif qual.isNumQual:
        sw.isNumQual = qual.isNumQual

    if ( ( not sw.setClips ) and
         sw.lengthsProvided and
         ( fastq1.seq != fastq1.seqOrig.strip() or
           ( fastq2 and
             fastq2.seq != fastq2.seqOrig.strip() ) or
           ( fastq3 and
             fastq3.seq != fastq3.seqOrig.strip() ) or
           ( fastq4 and
             fastq4.seq != fastq4.seqOrig.strip() ) or
           ( fastq5 and
             fastq5.seq != fastq5.seqOrig.strip() ) or
           ( fastq6 and
             fastq6.seq != fastq6.seqOrig.strip() ) ) ):
        if fastq6:
            statusWriter.outputInfo("Setting left and right clips based on {},{} or {},{} or {},{} or {},{} or {},{} or {},{}"
                                    .format(fastq1.defline.name, fastq1.filename, fastq2.defline.name, fastq2.filename, fastq3.defline.name, fastq3.filename, fastq4.defline.name, fastq4.filename, fastq5.defline.name, fastq5.filename, fastq6.defline.name, fastq6.filename) )
        elif fastq5:
            statusWriter.outputInfo("Setting left and right clips based on {},{} or {},{} or {},{} or {},{} or {},{}"
                                    .format(fastq1.defline.name, fastq1.filename, fastq2.defline.name, fastq2.filename, fastq3.defline.name, fastq3.filename, fastq4.defline.name, fastq4.filename, fastq5.defline.name, fastq5.filename) )
        elif fastq4:
            statusWriter.outputInfo("Setting left and right clips based on {},{} or {},{} or {},{} or {},{}"
                                    .format(fastq1.defline.name, fastq1.filename, fastq2.defline.name, fastq2.filename, fastq3.defline.name, fastq3.filename, fastq4.defline.name, fastq4.filename) )
        elif fastq3 and not fastq2:
            statusWriter.outputInfo("Setting left and right clips based on {},{} or {},{}"
                                    .format(fastq1.defline.name, fastq1.filename, fastq3.defline.name, fastq3.filename) )
        elif fastq3:
            statusWriter.outputInfo("Setting left and right clips based on {},{} or {},{} or {},{}"
                                    .format(fastq1.defline.name, fastq1.filename, fastq2.defline.name, fastq2.filename, fastq3.defline.name, fastq3.filename) )
        elif fastq2:
            statusWriter.outputInfo("Setting left and right clips based on {},{} or {},{}"
                                    .format(fastq1.defline.name, fastq1.filename, fastq2.defline.name, fastq2.filename) )
        else:
            statusWriter.outputInfo("Setting left and right clips based on {},{}"
                                    .format(fastq1.defline.name, fastq1.filename) )
        sw.setClips = True

    if ( fastq1 and fastq2 and
         fastq1.defline.readNum and
         fastq1.defline.readNum == fastq2.defline.readNum and
         fastq1.seq == fastq2.seq ):
        sw.duplicateReads = True
        if not fastq1.filename in fileWithDupReads:
            if fastq1.filename == fastq2.filename:
                statusWriter.outputInfo("Duplicate reads exist within {}".format(fastq1.filename) )
                fileWithDupReads[fastq1.filename] = 1
            else:
                statusWriter.outputInfo("Duplicate reads exist between {} and {}".format(fastq1.filename, fastq2.filename) )
                fileWithDupReads[fastq1.filename] = 1
                fileWithDupReads[fastq2.filename] = 1

############################################################
# Ensure lengths of seq and qual are the same
############################################################

def checkFastqLengths ( fastq ):
    if fastq.eof:
        return False
    else:
        if fastq.length != fastq.lengthQual:

            if fastq.lengthQual >= 5:
                fastq.addToBadLines(fastq.qual,fastq.lineCount)
            elif fastq.length >= 5:
                fastq.addToBadLines(fastq.seq,fastq.lineCount)
            elif len(fastq.prevLine) >= 5:
                fastq.addToBadLines(fastq.prevLine,fastq.lineCount)

            # Looks for absence of quality

            if ( not fastq.qual or
                 ( fastq.qual[0:2] == fastq.defline.deflineString[0:2] and
                   fastq.isDeflineString(fastq.qual) ) ):

                if fastq.qual:
                    fastq.savedDeflineString = fastq.qual

                if sw.isNumQual:
                    fastq.qual = "30"
                    fastq.qual += " 30" * (fastq.length - 1)
                else:
                    fastq.qual = fastq.defaultQual * fastq.length

                if ( not sw.maxErrorOutput or
                     sw.errorCount < sw.maxErrorOutput ):
                    statusWriter.outputWarning("Qual not present for spot {} in {} at line {} - created qual containing all values of 30"
                                               .format(fastq.defline.name,fastq.filename,fastq.lineCount) )

            # Addresses seq length > qual length

            elif fastq.length > fastq.lengthQual:
                delta = fastq.length - fastq.lengthQual
                if sw.isNumQual:
                    fastq.qual += " 30" * delta
                else:
                    fastq.qual += fastq.defaultQual * delta

                if ( not sw.maxErrorOutput or
                     sw.errorCount < sw.maxErrorOutput ):
                    statusWriter.outputWarning("Qual length less than seq length for spot {} in {} at line {} - extended qual to seq length with values of 30"
                                               .format(fastq.defline.name,fastq.filename,fastq.lineCount) )

            # Addresses qual length > seq length

            elif fastq.length < fastq.lengthQual:
                if sw.isNumQual:
                    qualChunks = fastq.qual.split()
                    delim = " "
                    fastq.qual = delim.join(qualChunks[0:fastq.length])
                else:
                    fastq.qual = fastq.qual[0:fastq.length]

                if ( not sw.maxErrorOutput or
                     sw.errorCount < sw.maxErrorOutput ):
                    statusWriter.outputWarning("Qual length greater than seq length for spot {} in {} at line {} - truncated qual to seq length"
                                               .format(fastq.defline.name,fastq.filename,fastq.lineCount) )

            fastq.lengthQual = fastq.length
            return True
        else:
            return False

############################################################
# Read a fastq file or a fastq file pair
############################################################

def processFastqSpots ( determineOffsetAndClip, fastq1, fastq2=None, fastq3=None, fastq4=None, fastq5=None, fastq6=None ):

    sw.errorCount = 0
    isEightLine = False
    isEightLine2 = False
    isMultiLineEightLine = False
    if ( fastq2 and
         fastq1.filename == fastq2.filename and
         not sw.concatPairFiles ):
        isEightLine = True
        if fastq1.isMultiLine:
            isMultiLineEightLine = True
    if ( fastq2 and
         fastq3 and
         fastq2.filename == fastq3.filename):
        isEightLine2 = True

    # Read from fastq files

    while True:

        if not fastq1.eof:
            fastq1.read()
            if ( fastq1.qualHandle and
                 fastq1.seqQualDeflineMismatch ):
                statusWriter.outputErrorAndExit( "Mismatch between seq and qual deflines ... filename,lineCount,deflineStringSeq,lineCountQual,deflineStringQual\t{}\t{}\t{}\t{}\t{}"
                                                 .format(fastq1.filename,fastq1.lineCount,fastq1.deflineStringSeq,fastq1.lineCountQual,fastq1.deflineStringQual))

        if ( fastq2 and
             not fastq2.eof ):

            if isEightLine:
                fastq2.lineCount = fastq1.lineCount
            if isMultiLineEightLine:
                fastq2.savedDeflineStringSeq = fastq1.savedDeflineStringSeq

            fastq2.read()
            if ( fastq2.qualHandle and
                 fastq2.seqQualDeflineMismatch ):
                statusWriter.outputErrorAndExit( "Mismatch between seq and qual deflines ... filename,lineCount,deflineStringSeq,lineCountQual,deflineStringQual\t{}\t{}\t{}\t{}\t{}"
                                                 .format(fastq2.filename,fastq2.lineCount,fastq2.deflineStringSeq,fastq1.lineCountQual,fastq2.deflineStringQual))

            if isEightLine:
                fastq1.lineCount = fastq2.lineCount
            elif isEightLine2:
                fastq3.lineCount = fastq2.lineCount
            if isMultiLineEightLine:
                fastq1.savedDeflineStringSeq = fastq2.savedDeflineStringSeq

        # Nanopore 2D fastq - Assuming this will end sooner than fastq1 or fastq2
        # Could all be the same file, too.
        # Alternatively, now allowing for third fastq file that is not nanopore 2D
        # Not allowing orphans in fastq3 and fastq4 and fastq5 and fastq6
        # Not allowing split seq/qual for fastq3/4

        if ( fastq3 and
             not fastq3.eof ):
            fastq3.read()
            if isEightLine2:
                fastq2.lineCount = fastq3.lineCount

        if ( fastq4 and
             not fastq4.eof ):
            fastq4.read()

        if ( fastq5 and
             not fastq5.eof ):
            fastq5.read()

        if ( fastq6 and
             not fastq6.eof ):
            fastq6.read()

        # Exit if at end of fastq1, fastq2 (if it exists), fastq3 (if it exists), and fastq4 (if it exists)
        # Possible to have empty quality so checking spot name, too.

        if ( fastq1.eof and
             ( fastq2 is None or
               fastq2.eof ) and
             ( fastq3 is None or
               fastq3.eof ) and
             ( fastq4 is None or
               fastq4.eof ) and
             ( fastq5 is None or
               fastq5.eof ) and
             ( fastq6 is None or
               fastq6.eof ) ):
            break

        # Accumulate extent of values for offset determination

        elif determineOffsetAndClip:
            updateQualRangeAndClipStatus ( fastq1, fastq2, fastq3, fastq4, fastq5, fastq6 )
            if ( fastq1.spotCount == 100000 or
                 fastq1.eof):
                break

        # Write spot for archive generation
        # Process extras in one file as fragments

        else:

            # Check for emergence of orphan reads (due to corruption for example - N/A if fastq3 or more)

            if ( not ( sw.ignoreNames or sw.ignAndDiscardNames or fastq3 ) and
                 not sw.orphanReads and
                 fastq1 and
                 fastq2 and
                 ( not fastq1.eof ) and
                 ( not fastq2.eof ) and
                 fastq1.defline.name != fastq2.defline.name ):
                sw.orphanReads = True
                sw.allowEarlyFileEnd = True
                if fastq1.filename == fastq2.filename:
                    statusWriter.outputInfo("File {} contains paired reads with orphan reads starting at line {}".format(fastq1.filename,fastq1.lineCount) )
                else:
                    statusWriter.outputInfo("File {} and {} are paired with orphan reads starting at line {} and line {} respectively".format(fastq1.filename,fastq2.filename,fastq1.lineCount,fastq2.lineCount) )

            # If fastq2 and at end of fastq2, set fastq2 to None
            # Set new unchanging values

            if ( fastq2 and
                 fastq2.eof ):
                if ( fastq3 and
                     ( not fastq1.eof ) and
                        ( not fastq1.defline.platform == 'NANOPORE') ):
                    if sw.allowEarlyFileEnd:
                        break
                    else:
                        statusWriter.outputErrorAndExit("Expected matching deflines but {} ended early at line {}. Use '--allowEarlyFileEnd' to allow load to finish."
                                                        .format(fastq2.filename,fastq2.lineCount) )
                else:
                    if sw.allowEarlyFileEnd or fastq1.defline.platform == 'NANOPORE':
                        fastq2 = None
                        sw.setUnchangingSpotValues()
                    else:
                        statusWriter.outputErrorAndExit("{} ended early at line {}. Use '--allowEarlyFileEnd' to allow load to proceed with orphan reads."
                                                        .format(fastq2.filename,fastq2.lineCount) )

            # If fastq2 and at end of fastq1, and set fastq1 to fastq2
            # and fastq2 to None. Set new unchanging values

            elif ( fastq2 and
                   fastq1.eof ):
                if ( fastq3 and
                     ( not fastq2.eof ) and
                     ( not fastq1.defline.platform == 'NANOPORE') ):
                    if sw.allowEarlyFileEnd:
                        break
                    else:
                        statusWriter.outputErrorAndExit("Expected matching deflines but {} ended early at line {}. Use '--allowEarlyFileEnd' to allow load to finish."
                                                        .format(fastq1.filename,fastq1.lineCount) )
                else:
                    if sw.allowEarlyFileEnd or fastq1.defline.platform == 'NANOPORE':
                        fastq1 = fastq2
                        fastq2 = None
                        sw.setUnchangingSpotValues()
                    else:
                        statusWriter.outputErrorAndExit("{} ended early at line {}. Use '--allowEarlyFileEnd' to allow load to proceed with orphan reads."
                                                        .format(fastq1.filename,fastq1.lineCount) )

            elif ( fastq3 and
                   fastq3.eof and
                   ( not fastq1.eof ) and
                   ( not fastq1.defline.platform == 'NANOPORE') ):
                if sw.allowEarlyFileEnd:
                    break
                else:
                    statusWriter.outputErrorAndExit("Expected matching deflines but {} ended early at line {}. Use '--allowEarlyFileEnd' to allow load to finish."
                                                    .format(fastq3.filename,fastq3.lineCount) )
            elif ( fastq4 and
                   fastq4.eof and
                   ( not fastq1.eof ) ):
                if sw.allowEarlyFileEnd:
                    break
                else:
                    statusWriter.outputErrorAndExit("Expected matching deflines but {} ended early at line {}. Use '--allowEarlyFileEnd' to allow load to finish."
                                                    .format(fastq4.filename,fastq4.lineCount) )

            elif ( fastq5 and
                   fastq5.eof and
                   ( not fastq1.eof ) ):
                if sw.allowEarlyFileEnd:
                    break
                else:
                    statusWriter.outputErrorAndExit("Expected matching deflines but {} ended early at line {}. Use '--allowEarlyFileEnd' to allow load to finish."
                                                    .format(fastq5.filename,fastq5.lineCount) )

            elif ( fastq6 and
                   fastq6.eof and
                   ( not fastq1.eof ) ):
                if sw.allowEarlyFileEnd:
                    break
                else:
                    statusWriter.outputErrorAndExit("Expected matching deflines but {} ended early at line {}. Use '--allowEarlyFileEnd' to allow load to finish."
                                                    .format(fastq6.filename,fastq6.lineCount) )

            # Check for  mismatch between seq and qual lengths

            mismatchLineCount = None
            mismatchFileName = None
            mismatchFastq = None
            if checkFastqLengths ( fastq1 ):
                sw.errorCount += 1
                mismatchLineCount = fastq1.lineCount
                mismatchFileName = fastq1.filename
                mismatchFastq = fastq1

            if ( fastq2 and
                 checkFastqLengths ( fastq2 ) ):
                sw.errorCount += 1
                mismatchLineCount = fastq2.lineCount
                mismatchFileName = fastq2.filename
                mismatchFastq = fastq2

            if ( fastq3 and
                 checkFastqLengths ( fastq3 ) ):
                sw.errorCount += 1
                mismatchLineCount = fastq3.lineCount
                mismatchFileName = fastq3.filename
                mismatchFastq = fastq3

            if ( fastq4 and
                 checkFastqLengths ( fastq4 ) ):
                sw.errorCount += 1
                mismatchLineCount = fastq4.lineCount
                mismatchFileName = fastq4.filename
                mismatchFastq = fastq4

            if ( fastq5 and
                 checkFastqLengths ( fastq5 ) ):
                sw.errorCount += 1
                mismatchLineCount = fastq5.lineCount
                mismatchFileName = fastq5.filename
                mismatchFastq = fastq5

            if ( fastq6 and
                 checkFastqLengths ( fastq6 ) ):
                sw.errorCount += 1
                mismatchLineCount = fastq6.lineCount
                mismatchFileName = fastq6.filename
                mismatchFastq = fastq6

            if ( sw.maxErrorCount and
                 sw.errorCount > sw.maxErrorCount ):
                mismatchFastq.processBadLines()
                if mismatchFastq.corruptionStart:
                    statusWriter.outputErrorAndExit( "Excessive errors occurred due to repetitive corruption starting in {} at line {}"
                                                     .format(mismatchFileName, mismatchFastq.corruptionStart ) )
                else:
                    statusWriter.outputErrorAndExit( "Excessive errors: {} cases occurred where sequence length does not match quality length (last case in {} at line {})"
                                                     .format(sw.errorCount,mismatchFileName,mismatchLineCount) )

            # Write spot to archive

            sw.writeSpot ( fastq1, fastq2, fastq3, fastq4, fastq5, fastq6 )

############################################################
# Process files (determine offset if specified)
############################################################

def processFiles (determineOffsetAndClip):

    # Open spot writer for actual archive generation

    if not determineOffsetAndClip:
        sw.openGeneralWriter()
        sw.writingToArchive = True
    else:
        sw.totalSize = 0
        sw.totalSize2= 0
        sw.totalSize3= 0
        sw.totalSize4= 0

    # Process each file or file pair/triplet/quad

    fastq1 = None
    fastq2 = None
    fastq3 = None
    fastq4 = None
    fastq5 = None
    fastq6 = None

    for filename1 in sorted(filePaths):

        if filename1 in fileSkip:
            continue

        else:

            # Set spotGroup if necessary

            if sw.fileSpotGroupsProvided:
                sw.spotGroup = fileSpotGroups[filename1]
            elif sw.useFilenameForSG:
                sw.spotGroup = FastqReader.getFilenameTrunc(filename1)

            # Retain truncated filename if adding file name to spot name

            if sw.addFileToName:
                sw.filenameTrunc = FastqReader.getFilenameTrunc(filename1)

            # Get fastq instances

            fastq1 = getFastqInstance(filename1)
            fastq1.read()
            if determineOffsetAndClip:
                if not fastq1.defline.sharqCompatible:
                    sw.useSharq = False
                #  Determine total file size for first files in pair sets while also determining quality offset/clip
                sw.totalSize += fileSizes[fastq1.filename]
            else:
                statusWriter.outputInfo("File {} platform based on defline is {}".format(filename1,fastq1.defline.platform) )
                fastq1.outputStatus = True
            fastq1.restart()

            if filename1 in fileConcatStart:
                fastq1.concatOne = True
                fastq1.concatPos = fileConcatStart[filename1]

            fastq2 = None
            fastq3 = None
            fastq4 = None
            fastq5 = None
            fastq6 = None
            filename2 = ""
            filename3 = ""
            filename4 = ""
            filename5 = ""
            filename6 = ""

            if filename1 in fileReadPairs:
                filename2 = fileReadPairs[filename1]
                fastq2 = getFastqInstance(filename2)
                if not determineOffsetAndClip:
                    statusWriter.outputInfo("File {} platform based on defline is {}".format(filename2,fastq1.defline.platform) )
                    fastq2.outputStatus = True
                else:
                    sw.totalSize2 += fileSizes[fastq2.filename]

                if filename1 in fileConcatStart:
                    fastq2.handle = getFileHandle ( filename1 )
                    fastq2.handle.seek(fileConcatStart[filename1])
                    fastq2.concatTwo = True

            if filename1 in filePore2D:
                filename3 = filePore2D[filename1]
                fastq3 = getFastqInstance(filename3)
                if not determineOffsetAndClip:
                    statusWriter.outputInfo("File {} platform based on defline is {}".format(filename3,fastq1.defline.platform) )
                    fastq3.outputStatus = True

            elif filename1 in fileReadPairs2:
                filename3 = fileReadPairs2[filename1]
                fastq3 = getFastqInstance(filename3)
                if not determineOffsetAndClip:
                    statusWriter.outputInfo("File {} platform based on defline is {}".format(filename3,fastq1.defline.platform) )
                    fastq3.outputStatus = True
                else:
                    sw.totalSize3 += fileSizes[fastq3.filename]

            if filename1 in fileReadPairs3:
                filename4 = fileReadPairs3[filename1]
                fastq4 = getFastqInstance(filename4)
                if not determineOffsetAndClip:
                    statusWriter.outputInfo("File {} platform based on defline is {}".format(filename4,fastq1.defline.platform) )
                    fastq4.outputStatus = True
                else:
                    sw.totalSize4 += fileSizes[fastq4.filename]

            if filename1 in fileReadPairs4:
                filename5 = fileReadPairs4[filename1]
                fastq5 = getFastqInstance(filename5)
                if not determineOffsetAndClip:
                    statusWriter.outputInfo("File {} platform based on defline is {}".format(filename5,fastq1.defline.platform) )
                    fastq5.outputStatus = True

            if filename1 in fileReadPairs5:
                filename6 = fileReadPairs5[filename1]
                fastq6 = getFastqInstance(filename6)
                if not determineOffsetAndClip:
                    statusWriter.outputInfo("File {} platform based on defline is {}".format(filename6,fastq1.defline.platform) )
                    fastq6.outputStatus = True

            # Reset fastq1 and set unchanging values for some optimization for archive generation

            if not determineOffsetAndClip:
                sw.setUnchangingSpotValues(fastq2,fastq3)

            # Read fastq files

            processFastqSpots ( determineOffsetAndClip, fastq1, fastq2, fastq3, fastq4, fastq5, fastq6 )
            if not determineOffsetAndClip:
                if fastq1.discardCount > 0:
                    statusWriter.outputInfo( "{} lines were discarded from {}".format(fastq1.discardCount,fastq1.filename) )
                if fastq2 and not (fastq2.filename == fastq1.filename ) and fastq2.discardCount > 0:
                    statusWriter.outputInfo( "{} lines were discarded from {}".format(fastq2.discardCount,fastq2.filename) )
                if fastq3 and fastq3.discardCount > 0:
                    statusWriter.outputInfo( "{} lines were discarded from {}".format(fastq3.discardCount,fastq3.filename) )
                if fastq4 and fastq4.discardCount > 0:
                    statusWriter.outputInfo( "{} lines were discarded from {}".format(fastq4.discardCount,fastq4.filename) )
                if fastq5 and fastq5.discardCount > 0:
                    statusWriter.outputInfo( "{} lines were discarded from {}".format(fastq5.discardCount,fastq5.filename) )
                if fastq6 and fastq6.discardCount > 0:
                    statusWriter.outputInfo( "{} lines were discarded from {}".format(fastq6.discardCount,fastq6.filename) )

            # Output orphan spots (for now, traversing the files again
            # to be consistent with the order of latf-load). I could
            # just dump the orphan reads collected by the sw object.
            # Not handling orphan reads for non-2D fastq3 or fastq4

            if ( sw.orphanReads and
                 not determineOffsetAndClip ):

                statusWriter.outputInfo( "Reprocessing fastq files to output discovered orphan reads ... ")
                fastq1.outputStatus = False
                if fastq2:
                    fastq2.outputStatus = False
                if fastq3:
                    fastq3.outputStatus = False
                if fastq4:
                    fastq4.outputStatus = False
                if fastq5:
                    fastq5.outputStatus = False
                if fastq6:
                    fastq6.outputStatus = False

                # Read fastq files again (keeping order consistent with latf-load)
                # Possibility of orphan 2D reads exists

                sw.dumpOrphans = True
                fastq1.restart()
                processFastqSpots ( False, fastq1 )

                if ( fastq2 and
                     ( sw.concatPairFiles or
                       not fastq1.filename == fastq2.filename ) ) :
                    fastq2.restart()
                    if sw.concatPairFiles:
                        fastq2.handle.seek(fileConcatStart[filename1])
                    processFastqSpots ( False, fastq2 )

                # Reset dumpOrphans and read hashes for case where multiple
                # mixed files/pairs are processed (may not be the correct
                # approach)

                sw.dumpOrphans = False
                sw.pairedRead1 = {}
                sw.pairedRead2 = {}

                if sw.platformString == "NANOPORE":

                    sw.dumpOrphans2D = True
                    fastq1.restart()
                    processFastqSpots ( False, fastq1 )

                    if ( fastq3  and
                         not fastq1.filename == fastq3.filename ) :
                        fastq3.restart()
                        processFastqSpots ( False, fastq3 )

                    sw.dumpOrphans2D = False
                    sw.nanopore2Dread = {}

            # Indicate quality range for file

            if determineOffsetAndClip:

                if ( sw.minQual == 1000 and
                     sw.maxQual == 0 ):
                    statusWriter.outputInfo( "Qual range not determined for file {}"
                                             .format(filename1) )
                elif ( not filename2 or
                     filename1 == filename2 ):
                    statusWriter.outputInfo( "Qual range using 0 or 33 offset for file {} is {},{}"
                                             .format(filename1,sw.minQual,sw.maxQual) )
                elif ( filename6 and
                       filename1 in fileReadPairs5 ):
                    statusWriter.outputInfo( "Qual range using 0 or 33 offset for files {}, {}, {}, {}, {}, and {} is {},{}"
                                             .format(filename1,filename2,filename3,filename4,filename5,filename6,sw.minQual,sw.maxQual) )

                elif ( filename5 and
                       filename1 in fileReadPairs4 ):
                    statusWriter.outputInfo( "Qual range using 0 or 33 offset for files {}, {}, {}, {}, and {} is {},{}"
                                             .format(filename1,filename2,filename3,filename4,filename5,sw.minQual,sw.maxQual) )

                elif ( filename4 and
                       filename1 in fileReadPairs3 ):
                    statusWriter.outputInfo( "Qual range using 0 or 33 offset for files {}, {}, {}, and {} is {},{}"
                                             .format(filename1,filename2,filename3,filename4,sw.minQual,sw.maxQual) )

                elif ( filename3 and
                       ( filename1 in fileReadPairs2 or
                         filename1 in filePore2D ) ):
                    statusWriter.outputInfo( "Qual range using 0 or 33 offset for files {}, {}, and {} is {},{}"
                                             .format(filename1,filename2,filename3,sw.minQual,sw.maxQual) )
                else :
                    statusWriter.outputInfo( "Qual range using 0 or 33 offset for files {} and {} is {},{}"
                                             .format(filename1,filename2,sw.minQual,sw.maxQual) )

            # Indicate end of files for archive generation

            elif ( not filename2 or
                   filename1 == filename2 ):
                statusWriter.outputInfo( "End of file {}".format(filename1) )
                sw.cumulativeBytesRead += fileSizes[fastq1.filename]

            elif ( filename6 and
                   filename1 in fileReadPairs5 ):
                statusWriter.outputInfo( "End of files {}, {}, {}, {}, {}, and {}".format(filename1,filename2,filename3,filename4,filename5,filename6) )
                sw.cumulativeBytesRead += fileSizes[fastq1.filename]
                sw.cumulativeBytesRead2 += fileSizes[fastq2.filename]

            elif ( filename5 and
                   filename1 in fileReadPairs4 ):
                statusWriter.outputInfo( "End of files {}, {}, {}, {}, and {}".format(filename1,filename2,filename3,filename4,filename5) )
                sw.cumulativeBytesRead += fileSizes[fastq1.filename]
                sw.cumulativeBytesRead2 += fileSizes[fastq2.filename]

            elif ( filename4 and
                   filename1 in fileReadPairs3 ):
                statusWriter.outputInfo( "End of files {}, {}, {}, and {}".format(filename1,filename2,filename3,filename4) )
                sw.cumulativeBytesRead += fileSizes[fastq1.filename]
                sw.cumulativeBytesRead2 += fileSizes[fastq2.filename]

            elif ( filename3 and
                   ( filename1 in fileReadPairs2 or
                     filename1 in filePore2D ) ):
                statusWriter.outputInfo( "End of files {}, {}, and {}".format(filename1,filename2,filename3) )
                sw.cumulativeBytesRead += fileSizes[fastq1.filename]
                sw.cumulativeBytesRead2 += fileSizes[fastq2.filename]

            else:
                statusWriter.outputInfo( "End of files {} and {}".format(filename1,filename2) )
                sw.cumulativeBytesRead += fileSizes[fastq1.filename]
                sw.cumulativeBytesRead2 += fileSizes[fastq2.filename]

    if determineOffsetAndClip:
        if sw.isNumQual:

            # Take into account color space

            if ( sw.isColorSpace and
                 sw.minQual == -1 ):
                sw.minQual = 0

            # Set log odds first if necessary (setQualityOffset dependency)
            # '-1' only is likely used for dot or N qualities

            if sw.minQual < -1:
                sw.logOdds = True
                statusWriter.outputInfo( "Logodds quality type" )
            elif sw.minQual == -1:
                sw.changeNegOneQual = True

            # Set offset to 0

            sw.setQualityOffset ( 0, False )
            statusWriter.outputInfo( "Treating quality as numeric" )

        elif ( ( sw.minQual > 25 and
                 sw.maxQual > 45 ) ) :

            # Set log odds first if necessary (setQualityOffset dependency)
            # '-1' only is likely used for dot or N qualities

            if ( sw.minQual + 33 - 64 ) < -1 :
                sw.logOdds = True
                statusWriter.outputInfo( "Logodds quality type2" )
            elif ( sw.minQual + 33 - 64 ) == -1:
                sw.changeNegOneQual = True

            # Set offset to 64

            sw.setQualityOffset ( 64, False )
            statusWriter.outputInfo( "Treating quality as ascii with an offset of 64" )

        else:
            sw.setQualityOffset ( 33, False )
            statusWriter.outputInfo( "Treating quality as ascii with an offset of 33" )

        if ( fastq1.defline.retainAltNum or
             ( fastq2 and fastq2.defline.retainAltNum ) or
             ( fastq3 and fastq3.defline.retainAltNum ) or
             ( fastq4 and fastq4.defline.retainAltNum ) or
             ( fastq5 and fastq5.defline.retainAltNum ) or
             ( fastq6 and fastq6.defline.retainAltNum ) or
             ( fastq1.defline.deflineType == fastq1.defline.ALTNUM and
               not fastq2 and
               not fastq1.filename in fileReadPairs ) ):
            sw.retainAltNum = True

    else:
        sw.gw = None # close stream and flush

############################################################
# Extract fastq from AB1 for processing
############################################################
def extractFastqFromAb1 (filename):
    filePathAb1 = filePaths[filename]
    filePathFastq = filePathAb1.replace('.ab1', '.fastq')
    progDir = os.path.dirname(os.path.abspath(__file__))
    extractFastqPath = "{}/extract_fastq".format(progDir)
    command = [extractFastqPath]
    command.append("-output")
    command.append("{}".format(filePathFastq))
    command.append(filePathAb1)
    statusWriter.outputInfo("     {}".format(command))
    try:
        subprocess.check_call(command)
        statusWriter.outputInfo("ab1 extract_fastq conversion succeeded")
    except subprocess.CalledProcessError as e:
        statusWriter.outputErrorAndExit("ab1 extract_fastq failed with return code {}".format(e.returncode))
    return filePathFastq

############################################################
# Open file handle for filename
############################################################

def getFileHandle (filename):
    handle = None

    filePath = filePaths[filename]
    if " " in filePath:
        filePath = r'%s' % filePath

    if filename.endswith('.ab1'):
        filePath = extractFastqFromAb1(filename)
        filePaths[filename] = filePath

    fileSizes[filename] = os.path.getsize(filePath)

    result = subprocess.check_output( ["file", "-bikL", "--", filePath ] )

    if b"gzip" in result:
        fileObj = open(filePaths[filename], 'rb')
        fileObjects[filename] = fileObj
        cmd = "zcat '" + filePath +  "' | file -bikL -- -"
        result = subprocess.check_output( ["bash", "-c", cmd])
        if b"charset=utf-16le" in result:
            handle = gzip.open ( fileObj, 'rt', encoding='utf-16le' )
            sw.removeBadChars = True # Ensures Byte Order Mark is removed
            sw.useSharq = False
        else:
            handle = gzip.open ( fileObj, 'rt' )
        sw.isCompressed = 1
    elif b"bzip" in result:
        fileObj = open(filePaths[filename], 'rb')
        fileObjects[filename] = fileObj
        cmd = "bzcat '" + filePath +  "' | file -bikL -- -"
        result = subprocess.check_output( ["bash", "-c", cmd])
        if b"charset=utf-16le" in result:
            handle = bz2.open ( fileObj, 'rt', encoding='utf-16le' )
            sw.removeBadChars = True # Ensures Byte Order Mark is removed
            sw.useSharq = False
        else:
            handle = bz2.open( fileObj, 'rt' )
        sw.isCompressed = 1
    elif b"application/x-xz" in result:
        fileObj = open(filePaths[filename], 'rb')
        fileObjects[filename] = fileObj
        handle = lzma.open( fileObj, 'rt' )
        sw.isCompressed = 1
        sw.useSharq = False
    elif b"text/rtf" in result:
        fileObj = open(filePaths[filename], 'rb')
        fileObjects[filename] = fileObj
        doc = Rtf15Reader.read(fileObj)
        handle = PlaintextWriter.write(doc)
        sw.useSharq = False
    elif b"application/zip;" in result:
        statusWriter.outputErrorAndExit( "Unable to handle files in a zip archive ... {} is {}".format(filePath,result) )
    elif b"application/x-compress;" in result:
        statusWriter.outputErrorAndExit( "Unable to handle x-compress files ... {} is {}".format(filePath,result) )
    elif b"charset=utf-16le" in result:
        handle = open ( filePaths[filename], encoding='utf-16le' )
        fileObjects[filename] = handle
        sw.removeBadChars = True # Ensures Byte Order Mark is removed
        sw.useSharq = False
    else: # If neither gz or bz2 or xz,then just try text
        handle = open ( filePaths[filename], encoding='utf-8', errors='ignore' )
        fileObjects[filename] = handle

    return handle

############################################################
# Set read types for 10x submission based on file names and
# read 1 lengths (i.e. technical vs biological)
############################################################
def set10xReadTypes():
    for filename1 in sorted(filePaths):

        type1 = None
        type2 = None
        type3 = None
        type4 = ""
        filenameR1 = None
        filenameR2 = None

        if filename1 in fileSkip:
            continue

        else:
            if ( re.search ( "[._-]I1[._]" , filename1 ) or
                 re.search ( "[._-]Index[._]" , filename1 ) ):
                type1 = "T"
            elif re.search ( "[._-]I2[._]" , filename1 ):
                type1 = "T"
            elif re.search ( "[._-]R1[._]" , filename1 ):
                type1 = "R1"
                filenameR1 = filename1
            elif re.search ( "[._-]R2[._]" , filename1 ):
                type1 = "R2"
                filenameR2 = filename1

            if filename1 in fileReadPairs:
                filename2 = fileReadPairs[filename1]
                if ( re.search ( "[._-]I1[._]" , filename2 ) or
                     re.search ( "[._-]Index[._]" , filename2 ) ):
                    type2 = "T"
                elif re.search ( "[._-]I2[._]" , filename2 ):
                    type2 = "T"
                elif re.search ( "[._-]R1[._]" , filename2 ):
                    type2 = "R1"
                    if filenameR1 is not None:
                        statusWriter.outputErrorAndExit( "Two R1 files are paired together ... {} and {}"
                                                         .format(filenameR1,filename2) )
                    else:
                        filenameR1 = filename2
                elif re.search ( "[._-]R2[._]" , filename2 ):
                    type2 = "R2"
                    if filenameR2 is not None:
                        statusWriter.outputErrorAndExit( "Two R2 files are paired together ... {} and {}"
                                                         .format(filenameR2,filename2) )
                    else:
                        filenameR2 = filename2

            if filename1 in fileReadPairs2:
                filename3 = fileReadPairs2[filename1]
                if ( re.search ( "[._-]I1[._]" , filename3 ) or
                     re.search ( "[._-]Index[._]" , filename3 ) ):
                    type3 = "T"
                elif re.search ( "[._-]I2[._]" , filename3 ):
                    type3 = "T"
                elif re.search ( "[._-]R1[._]" , filename3 ):
                    type3 = "R1"
                    if filenameR1 is not None:
                        statusWriter.outputErrorAndExit( "Two R1 files are paired together ... {} and {}"
                                                         .format(filenameR1,filename3) )
                    else:
                        filenameR1 = filename3
                elif re.search ( "[._-]R2[._]" , filename3 ):
                    type3 = "R2"
                    if filenameR2 is not None:
                        statusWriter.outputErrorAndExit( "Two R2 files are paired together ... {} and {}"
                                                         .format(filenameR2,filename3) )
                    else:
                        filenameR2 = filename3

            if filename1 in fileReadPairs3: # Assuming biological for 4th file in set, if it exists
                type4 = "B"

            if type1 and type2 and type3:
                fastqR1 = getFastqInstance(filenameR1)
                fastqR1.read()

                if fastqR1.length > 40:
                    typeR1 = "B"
                else:
                    typeR1 = "T"

                if type1 == "R1":
                    type1 = typeR1
                elif type2 == "R1":
                    type2 = typeR1
                elif type3 == "R1":
                    type3 = typeR1

                # If R2 is not file 4

                if filenameR2 is not None:
                    fastqR2 = getFastqInstance(filenameR2)
                    fastqR2.read()

                    if fastqR2.length > 40:
                        typeR2 = "B"
                    else:
                        typeR2 = "T"

                    if type1 == "R2":
                        type1 = typeR2
                    elif type2 == "R2":
                        type2 = typeR2
                    elif type3 == "R2":
                        type3 = typeR2

                readTypes = type1 + type2 + type3 + type4
                sw.setReadTypes(readTypes)
                break

############################################################
# Check provided pair file sets
############################################################

def checkProvidedPairs():
    for filename1 in fileReadPairs:
        filename2 = fileReadPairs[filename1]
        if filename1 != filename2:
            fileSkip[filename2] = 1
            type1 = fileTypes[filename1]
            type2 = fileTypes[filename2]
            if (not sw.mixedTypes and
                    type1 != type2):
                statusWriter.outputErrorAndExit("Paired files have different file types ... {} is {}, {} is {}"
                                                .format(filename1, type1, filename2, type2))

    for filename1 in fileReadPairs2:
        filename3 = fileReadPairs2[filename1]
        if filename1 != filename3:
            fileSkip[filename3] = 1
            type1 = fileTypes[filename1]
            type3 = fileTypes[filename3]
            if (not sw.mixedTypes and
                    type1 != type3):
                statusWriter.outputErrorAndExit("Triplet files have different file types ... {} is {}, {} is {}"
                                                .format(filename1, type1, filename3, type3))

    for filename1 in fileReadPairs3:
        filename4 = fileReadPairs3[filename1]
        if filename1 != filename4:
            fileSkip[filename4] = 1
            type1 = fileTypes[filename1]
            type4 = fileTypes[filename4]
            if (not sw.mixedTypes and
                    type1 != type4):
                statusWriter.outputErrorAndExit("Quad files have different file types ... {} is {}, {} is {}"
                                                .format(filename1, type1, filename4, type4))

    for filename1 in fileReadPairs4:
        filename5 = fileReadPairs4[filename1]
        if filename1 != filename5:
            fileSkip[filename5] = 1
            type1 = fileTypes[filename1]
            type5 = fileTypes[filename5]
            if (not sw.mixedTypes and
                    type1 != type5):
                statusWriter.outputErrorAndExit("Five file set has different file types ... {} is {}, {} is {}"
                                                .format(filename1, type1, filename5, type5))

    for filename1 in fileReadPairs5:
        filename6 = fileReadPairs5[filename1]
        if filename1 != filename6:
            fileSkip[filename6] = 1
            type1 = fileTypes[filename1]
            type6 = fileTypes[filename6]
            if (not sw.mixedTypes and
                    type1 != type6):
                statusWriter.outputErrorAndExit("Six file set has different file types ... {} is {}, {} is {}"
                                                .format(filename1, type1, filename6, type6))

############################################################
# Execute load with sharq
############################################################

def processFilesWithSharq ():
    statusWriter.outputInfo("Loading with sharq command line:")
    progDir = os.path.dirname(os.path.abspath(__file__))
    sharqPath = "{}/sharq".format(progDir)
    command = [sharqPath]
    if sw.typesProvided:
        command.append("--readTypes={}".format(sw.readTypeString))
    if sw.read1PairFiles:
        command.append("--read1PairFiles={}".format(sw.read1PairFiles))
    if sw.read2PairFiles:
        command.append("--read2PairFiles={}".format(sw.read2PairFiles))
    if sw.read3PairFiles:
        command.append("--read3PairFiles={}".format(sw.read3PairFiles))
    if sw.read4PairFiles:
        command.append("--read4PairFiles={}".format(sw.read4PairFiles))
    if sw.allowEarlyFileEnd:
        command.append("--allowEarlyFileEnd")
    if sw.useAndDiscardNames:
        command.append("--useAndDiscardNames")
    command.append("--log-level=info")
    command.append("--output={}".format(sw.outdir))
    for filename in filePaths:
        command.append("{}".format(filePaths[filename]))
    statusWriter.outputInfo("     {}".format(command))
    try:
        subprocess.check_call(command)
        statusWriter.outputInfo("sharq load succeeded")
    except subprocess.CalledProcessError as e:
        pass

############################################################
# Generate archive from provided fastq files
############################################################

def generateArchive():

    # Open files to be processed

    for filename in filePaths:
        fileHandles[filename] = getFileHandle ( filename )

    # Process file lists if provided

    if sw.read1PairFiles:
        processPairAndQualLists()

    # Set file types
    # (normal, singleLine, eightLine, multiLine, multiLineEightLine,
    #  seqQual, multiLineSeqQual, fasta, multiLineFasta, solid)

    setFileTypes()

    # Pair up files containing read pairs (unless pairs were provided)

    if not sw.pairFilesProvided:
        setFilePairs()

    # Check provided pair files

    else:
        checkProvidedPairs()

    # Check for properly characterized absolid or nanopore files

    checkForColorSpaceAndPlatform()

    # Assign barcodes to files if necessary

    if sw.fileSpotGroupsProvided:
        spotGroups = sw.fileSpotGroupsProvided.split(",")
        for filename in sorted(filePaths):
            if filename in fileSkip:
                continue
            else:
                fileSpotGroups[filename] = spotGroups.pop(0)

    # Check for unsupported combination of file sets

    if sw.mixedTypes:
        pass
    elif ( len ( fileReadPairs5 ) > 0 and
         ( ( len ( filePaths ) % 6 ) != 0 or
           len ( fileReadPairs5 ) != len ( fileReadPairs4 ) or
           len ( fileReadPairs5 ) != len ( fileReadPairs3 ) or
           len ( fileReadPairs5 ) != len ( fileReadPairs2 ) or
           len ( fileReadPairs5 ) != len ( fileReadPairs ) ) ):
        statusWriter.outputErrorAndExit("Mix of one/two/three/four/five and six file sets cannot be processed")
    elif ( not ( len ( fileReadPairs5 ) > 0 ) and
           len ( fileReadPairs4 ) > 0 and
         ( ( len ( filePaths ) % 5 ) != 0 or
           len ( fileReadPairs4 ) != len ( fileReadPairs3 ) or
           len ( fileReadPairs4 ) != len ( fileReadPairs2 ) or
           len ( fileReadPairs4 ) != len ( fileReadPairs ) ) ):
        statusWriter.outputErrorAndExit("Mix of one/two/three/four and five file sets cannot be processed")
    elif ( not ( len ( fileReadPairs4 ) > 0 ) and
           len ( fileReadPairs3 ) > 0 and
         ( ( len ( filePaths ) % 4 ) != 0 or
           len ( fileReadPairs3 ) != len ( fileReadPairs2 ) or
           len ( fileReadPairs3 ) != len ( fileReadPairs ) ) ):
        statusWriter.outputErrorAndExit("Mix of one/two/three and four file sets cannot be processed")
    elif ( not ( len ( fileReadPairs3 ) > 0 ) and
           len ( fileReadPairs2 ) > 0 and
         ( ( len ( filePaths ) % 3 ) != 0 or
           len ( fileReadPairs2 ) != len ( fileReadPairs ) ) ):
        statusWriter.outputErrorAndExit("Mix of one/two and three file sets cannot be processed")

    # Try to set readTypes for 10x submission

    if ( sw.platformString != "NANOPORE" and
         ( len ( sw.readTypes ) < 3 ) and
         ( len ( fileReadPairs2 ) > 0 ) ):
        set10xReadTypes()

    # Prescan to determine quality offset/logodds/etc. and sw.totalSize (of 1st files in each pair set)

    processFiles ( True )

    # Check for large file size (except for specific platforms)

    if not ( sw.platformString.lower() in ("illumina",
                                           "nanopore",
                                           "pacbio") or
             sw.ignAndDiscardNames or
             sw.useAndDiscardNames ) :

        oneGig = 1024 * 1024 * 1024
        if ( sw.isCompressed and
             ((sw.totalSize / oneGig) > 25 or
              (sw.totalSize2 / oneGig) > 25 or
              (sw.totalSize3 / oneGig) > 25 or
              (sw.totalSize4 / oneGig) > 25)):
            sw.setUseAndDiscardNames()
            statusWriter.outputInfo("Using and discarding names due to compressed file sizes.")
        elif ((sw.totalSize / oneGig) > 90 or
              (sw.totalSize2 / oneGig) > 90 or
              (sw.totalSize3 / oneGig) > 90 or
              (sw.totalSize4 / oneGig) > 90):
            sw.setUseAndDiscardNames()
            statusWriter.outputInfo("Using and discarding names due to file sizes.")

    # Read/Parse/Write fastq data

    if sw.isSharqCompatible():
        processFilesWithSharq()
    else:
        processFiles ( False )

############################################################
# Get general loader process id
############################################################

def getGeneralLoaderPid():
    pid = str ( os.getpid() )
    username = getpass.getuser()
    glPid = None
    with open(os.devnull, 'w') as devnull:
        for line in subprocess.check_output( ["/usr/sbin/lsof", "-b", "-u", username], stderr=devnull ).split(b'\n'):
            line = line.decode("utf-8")
            lsofChunks = line.split()
            if ( len ( lsofChunks) > 8 and
                 lsofChunks[0] in ('fastq-loa', 'python') and
                 lsofChunks[1] == pid and
                 lsofChunks[8] == 'pipe'):
                pipeId = lsofChunks[7]
                if pipeId is not None:
                    for line2 in subprocess.check_output( ["/usr/sbin/lsof", "-b", "-u", username], stderr=devnull ).split(b'\n'):
                        line2 = line2.decode("utf-8")
                        lsofChunks = line2.split()
                        if ( len ( lsofChunks) > 8 and
                             lsofChunks[0] == 'general-l' and
                             lsofChunks[7] == pipeId and
                             lsofChunks[8] == 'pipe'):
                            glPid = lsofChunks[1]
                            break
                break

    return glPid

############################################################
# Wait for general loader to finish
############################################################

def waitForGeneralLoader(glPid):
    if glPid is not None:
        glPidPath = "/proc/" + str(glPid) + "/cmdline"
        sys.stdout.flush()
        while True:
            time.sleep(0.1)
            if not os.path.exists(glPidPath):
                break
    while True:
        time.sleep(0.1)
        if ( not os.path.exists( sw.outdir ) or
             time.time() - os.path.getmtime(sw.outdir) > 20 ):
            break

############################################################
# Execute generic fastq loader
############################################################

def main():
    global sw,statusWriter
    sw = FastqSpotWriter()
    statusWriter = StatusWriter(version)
    sw.statusWriter = statusWriter

    sw.useSharq = True
    processArguments()

    glPid = getGeneralLoaderPid()

    if sw.profile:
        pr = cProfile.Profile()
        pr.enable()
        generateArchive()
        pr.disable()
        ps = pstats.Stats(pr, stream=sys.stderr).sort_stats('cumulative')
        ps.print_stats()
    else:
        generateArchive()

# Wait for general-loader to finish

    waitForGeneralLoader(glPid)

# Output load complete indication

    if sw.errorCount > 0:
        statusWriter.outputInfo( "{} cases occurred where quality values were truncated/augmented/created due to a mismatch with sequence length".format(sw.errorCount) )
    statusWriter.outputInfo("Load complete")
    if statusWriter.xmlLogHandle:
        statusWriter.closeXmlLog()

if __name__ == '__main__':
    main()
