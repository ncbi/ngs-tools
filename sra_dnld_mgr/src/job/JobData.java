/* ===========================================================================
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
=========================================================================== */
package job;

/* ---------------------------------------------------------------------------
    JobData ... derived from IniFile

    to_blocks( long )  ....... bined value into block-nr
    get_full_source()  ....... full path of source ( ACC or file )
    get_short_source() ....... just accession of filename without extension
    get_exportpath() ......... path where output has to be written to
    get_output_extension() ... format-specific extension for output
    get_output_filename() .... export_path + short_source + output_extension
    get_state() .............. the state the job is in ( enum JobState )
    is_runnable() ............ READY or PAUSED
    get_format() ............. output-format ( enum JobFormat )
    get_progress() ........... how far we have come ( row-id for BIO, byte for DNLD )
    get_max() ................ max. row-id for BIO or file-size for DNLD
    is_complete() ............ did we reach max?
    get_int_max() ............ max. value in blocks
    get_completion_str() ..... string showing how many percent done

--------------------------------------------------------------------------- */

import Bio.BioSpec;
import data.*;
import Bio.BioAccessionType;
import Bio.BioReadType;
import java.io.File;
import java.nio.file.*;
import java.text.*;
import java.util.concurrent.locks.*;

public final class JobData extends IniFile
{
    private static final String SOURCE = "SOURCE";
    private static final String EXPORTPATH = "EXPORTPATH";
    private static final String DOWNLOADPATH = "DOWNLOADPATH";
    private static final String STATE = "STATE";
    private static final String FORMAT = "FORMAT";
    private static final String SUBFORMAT = "SUBFORMAT";    
    private static final String PROGRESS = "PROGRESS";
    private static final String MAX = "MAX";
    private static final String MD5 = "MD5";
    private static final String BIO_TYPE = "BIO_TYPE";
    private static final String LINE_ENDING = "LINE_ENDING";
    private static final String LINE_WRAP = "LINE_WRAP";
    private static final String USE_LINE_WRAP = "USE_LINE_WRAP";
    private static final String FIXED_QUAL = "FIXED_QUAL";
    private static final String USE_FIXED_QUAL = "USE_FIXED_QUAL";
    private static final String FASTQ_DUMP_NAME = "FASTQ_DUMP_NAME";
    private static final String RUNTIME = "RUNTIME";
    private static final String REJECTED = "REJECTED";
    private static final String F_MIN_READ_LEN = "MIN_READ_LEN";
    private static final String F_USE_MIN_READ_LEN = "USE_MIN_READ_LEN";
    private static final String F_MAX_N = "MAX_N";
    private static final String F_USE_MAX_N = "USE_MAX_N";
    private static final String F_END_N = "END_N";
    private static final String F_USE_END_N = "USE_END_N";
    private static final String F_SPOTGROUP = "SPOTGROUP";
    private static final String F_USE_SPOTGROUP = "USE_SPOTGROUP";
    private static final String F_READ_TYPE = "READ_TYPE";
    private static final String F_USE_READ_TYPE = "USE_READ_TYPE";
    
    private static final LineEndings DFLT_LINE_ENDING = LineEndings.POSIX;
    private static final int DFLT_LINE_WRAP = 75;
    private static final boolean DFLT_USE_LINE_WRAP = true;
    private static final int DFLT_FIXED_QUAL = 30;
    private static final boolean DFLT_USE_FIXED_QUAL = false;
    private static final boolean DFLT_FASTQ_DUMP_NAME = false;
    private static final long DFLT_RUNTIME = 0;
    private static final long DFLT_REJECTED = 0;
    private static final int DFLT_F_MIN_READ_LEN = 50;
    private static final boolean DFLT_F_USE_MIN_READ_LEN = false;
    private static final int DFLT_F_MAX_N = 10;
    private static final boolean DFLT_F_USE_MAX_N = false;
    private static final int DFLT_F_END_N = 10;
    private static final boolean DFLT_F_USE_END_N = false;
    private static final String DFLT_F_SPOTGROUP = "";
    private static final boolean DFLT_F_USE_SPOTGROUP = false;
    private static final BioReadType DFLT_F_READ_TYPE = BioReadType.READ_TYPE_ALL;
    private static final boolean DFLT_F_USE_READ_TYPE = false;
    
    private static final int BLOCKSIZE = 1024;

    private final Lock lock;
    private final Lock state_lock;
    private final String sep;
    private long current_row;
    private long max_row;
    
    public static int to_blocks( long value )
    {
        long res = value;
        res /= BLOCKSIZE;
        return ( int ) res;
    }

    private String get_str_value( String key )
    {
        String res = "";
        if ( lock.tryLock() )
        {
            res = get_str( key, "" );
            lock.unlock();
        }
        return res;
    }

    public final String get_full_source()
    {
        return get_str_value( SOURCE );
    }
    
    public final String get_short_source()
    {
        String s = get_full_source();
        String fsep = FileSystems.getDefault().getSeparator();
        if ( s.contains( fsep ) )
        {
            Path p = FileSystems.getDefault().getPath( s );
            int n = p.getNameCount();
            s = p.getName( n - 1 ).toString();
            if ( s.contains( "." ) )
            {
                String x[] = s.split( "\\." );
                if ( x.length > 0 ) s = x[ 0 ];
            }
        }
        return s;
    }
   
    public final String get_exportpath() { return get_str_value( EXPORTPATH ); }
    public final String get_downloadpath() { return get_str_value( DOWNLOADPATH ); }
    
    public final String get_output_extension()
    {
        String res = "txt";
        JobFormat jf = get_format();
        switch( jf )
        {
            case DOWNLOAD            : res = "sra"; break;
            case FASTA               : res = "fasta.txt"; break;
            case FASTQ               : res = "fastq.txt"; break;
        }
        return res;
    }
    
    public final String get_output_filename()
    {
        JobFormat jf = get_format();
        BioAccessionType bt = get_bio_type();
        boolean is_dnld = ( jf == JobFormat.DOWNLOAD );
        boolean is_ref  = ( bt.equals( BioAccessionType.REF_SEQUENCE ) );
        if ( is_ref )
        {
            return String.format( "%s%s%s",
                    is_dnld ? get_downloadpath() : get_exportpath(),
                    sep,
                    get_short_source() );
        }
        return String.format( "%s%s%s.%s",
                is_dnld ? get_downloadpath() : get_exportpath(),
                sep,
                get_short_source(),
                get_output_extension() );
    }
    
    public final String get_output_filename( final int suffix_nr )
    {
        return String.format( "%s%s%s_%d.%s",
                get_exportpath(),
                sep,
                get_short_source(),
                suffix_nr, get_output_extension() );
    }

    public final JobState get_state()
    {
        boolean done = false;
        JobState res = JobState.INVALID;
        while ( !done )
        {
            if ( state_lock.tryLock() )
            {
                res = JobState.from_ordinal( get_int( STATE, JobState.INVALID.to_ordinal() ) );
                state_lock.unlock();
                done = true;
            }
            else
            {
                try { Thread.sleep( 10 ); }
                catch ( InterruptedException ex ) { }
            }
        }
        return res;
    }

    public final boolean is_runnable()
    {
        JobState js = get_state();
        return ( js == JobState.READY || js == JobState.PAUSED );
    }
    
    public final JobFormat get_format()
    {
        JobFormat res = JobFormat.INVALID;
        if ( lock.tryLock() )
        {
            res = JobFormat.from_ordinal( get_int( FORMAT, JobFormat.INVALID.to_ordinal() ) );
            lock.unlock();
        }
        return res;
    }

    public final JobSubFormat get_subformat()
    {
        JobSubFormat res = JobSubFormat.INVALID;
        if ( lock.tryLock() )
        {
            res = JobSubFormat.from_ordinal( get_int( SUBFORMAT, JobSubFormat.INVALID.to_ordinal() ) );
            lock.unlock();
        }
        return res;
   }
    
    public final long get_progress()
    {
        long res = 0;
        if ( lock.tryLock() )
        {
            res = current_row;
            if ( res < 0 ) res = get_long( PROGRESS, 0 );
            lock.unlock();
        }
        return res;
    }

    public final long get_max()
    {
        long res = 0;
        if ( lock.tryLock() )
        {
            res = max_row;
            if ( res < 0 ) res = get_long( MAX, 0 );
            lock.unlock();
        }
        return res;
    }

    public final LineEndings get_line_ending()
    { return LineEndings.from_ordinal( get_int( LINE_ENDING, DFLT_LINE_ENDING.to_ordinal() ) ); }
    public final String get_line_ending_str()
    { return LineEndings.from_ordinal( get_int( LINE_ENDING, DFLT_LINE_ENDING.to_ordinal() ) ).to_line_ending(); }

    public final int get_line_wrap() { return get_int( LINE_WRAP, DFLT_LINE_WRAP ); }
    public final Boolean get_use_line_wrap() { return get_bool( USE_LINE_WRAP, DFLT_USE_LINE_WRAP ); }
    public final int get_fixed_qual() { return get_int( FIXED_QUAL, DFLT_FIXED_QUAL ); }
    public final Boolean get_use_fixed_qual() { return get_bool( USE_FIXED_QUAL, DFLT_USE_FIXED_QUAL ); }
    public final Boolean get_fastq_dump_name() { return get_bool( FASTQ_DUMP_NAME, DFLT_FASTQ_DUMP_NAME ); }
    public final long get_runtime() { return get_long( RUNTIME, DFLT_RUNTIME ); }
    public final long get_rejected() { return get_long( REJECTED, DFLT_REJECTED ); }
    public final String get_md5() { return get_str_value( MD5 );}
    
    public final BioAccessionType get_bio_type() 
    { return BioAccessionType.from_ordinal( get_int( BIO_TYPE, BioAccessionType.INVALID.to_ordinal() ) ); }
    
    /*filter getters: */
    public final int get_min_read_len() { return get_int( F_MIN_READ_LEN, DFLT_F_MIN_READ_LEN ); }
    public final Boolean get_use_min_read_len() { return get_bool( F_USE_MIN_READ_LEN, DFLT_F_USE_MIN_READ_LEN ); }
    public final int get_max_N() { return get_int( F_MAX_N, DFLT_F_MAX_N ); }
    public final Boolean get_use_max_N() { return get_bool( F_USE_MAX_N, DFLT_F_USE_MAX_N ); }
    public final int get_end_N() { return get_int( F_END_N, DFLT_F_END_N ); }
    public final Boolean get_use_end_N() { return get_bool( F_USE_END_N, DFLT_F_USE_END_N ); }
    public final String get_spotgroup() { return get_str( F_SPOTGROUP, DFLT_F_SPOTGROUP ); }
    public final Boolean get_use_spotgroup() { return get_bool( F_USE_SPOTGROUP, DFLT_F_USE_SPOTGROUP ); }
    public final BioReadType get_bio_read_type()
    { return BioReadType.from_ordinal( get_int( F_READ_TYPE, BioReadType.INVALID.to_ordinal() ) ); }
    public final Boolean get_use_bio_read_type() { return get_bool( F_USE_READ_TYPE, DFLT_F_USE_READ_TYPE ); }
    
    public final boolean is_complete()
    {
        long v_max = get_max();
        return ( v_max > 0 && get_progress() >= v_max );
    }

    private static double get_completion( final long progress, final long max_value )
    {
        double res = 0.0;
        if ( progress > 0 && max_value > 0 )
        {
            res = progress;
            res /= max_value;
            res *= 100;
        }
        return res;
    }
    
    public final String get_completion_str()
    {
        long progress = get_progress();
        long max_value = get_max();
        DecimalFormat df = new DecimalFormat( "0" );
        df.setMaximumFractionDigits( 4 );
        String percent = df.format( get_completion( progress, max_value ) );
        return String.format( "%s %% ( %s of %s )",
                percent,
                NumberFormat.getInstance().format( progress ),
                NumberFormat.getInstance().format( max_value ) );
    }
    
    public final boolean is_dnld()
    {
        JobFormat jf = get_format();
        return jf.equals( JobFormat.DOWNLOAD );
    }

    public final boolean is_csra()
    {
        BioAccessionType bio = get_bio_type();
        return bio.equals( BioAccessionType.READ_COLLECTION_ALIGNED );
    }

/* -------------------------------------------------------------------------- */    
    
    private void set_str_value( final String key, final String value )
    {
        if ( lock.tryLock() )
        {
            set_str( key, value );
            lock.unlock();
        }
    }
    
    public final void set_source( final String value ) { set_str_value( SOURCE, value ); }
    public final void set_exportpath( final String value ) { set_str_value( EXPORTPATH, value ); }
    public final void set_downloadpath( final String value ) { set_str_value( DOWNLOADPATH, value ); }
    public final void set_line_ending( final LineEndings value ) {  set_int( LINE_ENDING, value.ordinal() ); }
    public final void set_line_wrap( final int value ) { set_int( LINE_WRAP, value ); }
    public final void set_use_line_wrap( final Boolean value ) { set_bool( USE_LINE_WRAP, value ); }
    public final void set_fixed_qual( final int value ) { set_int( FIXED_QUAL, value ); }
    public final void set_use_fixed_qual( final Boolean value ) { set_bool( USE_FIXED_QUAL, value ); }
    public final void set_fastq_dump_name( final Boolean value ) { set_bool( FASTQ_DUMP_NAME, value ); }
    public final void set_runtime( final long value ) { set_long( RUNTIME, value ); }
    public final void set_rejected( final long value ) { set_long( REJECTED, value ); }
    public final void set_bio_type( final BioAccessionType value ) { set_int( BIO_TYPE, value.ordinal() ); }
    
    /*filter setters: */
    public final void set_min_read_len( final int value ) { set_int( F_MIN_READ_LEN, value ); }
    public final void set_use_min_read_len( final Boolean value ) { set_bool( F_USE_MIN_READ_LEN, value ); }
    public final void set_max_N( final int value ) { set_int( F_MAX_N, value ); }
    public final void set_use_max_N( final Boolean value ) { set_bool( F_USE_MAX_N, value ); }
    public final void set_end_N( final int value ) { set_int( F_END_N, value ); }
    public final void set_use_end_N( final Boolean value ) { set_bool( F_USE_END_N, value ); }
    public final void set_spotgroup( final String value ) { set_str( F_SPOTGROUP, value ); }
    public final void set_use_spotgroup( final Boolean value ) { set_bool( F_USE_SPOTGROUP, value ); }
    public final void set_bio_read_type( final BioReadType value ) { set_int( F_READ_TYPE, value.ordinal() ); }
    public final void set_use_bio_read_type( final Boolean value ) { set_bool( F_USE_READ_TYPE, value ); }
    
    /* set set the values from a spec-class ( name, type, rows ) */
    public final void set_spec( final BioSpec spec )
    {
        set_source( spec.get_accession() );
        set_bio_type( spec.get_type() );
        set_max( spec.get_count(), false );
    }
    
    public final void set_state( final JobState value, final Boolean save )
    {
        if ( state_lock.tryLock() )
        {
            set_int( STATE, value.to_ordinal() );
            state_lock.unlock();
        }
        if ( save ) store();
    }

    public final void set_format( final JobFormat value )
    {
        if ( lock.tryLock() )
        {
            set_int( FORMAT, value.to_ordinal() );
            lock.unlock();
        }
    }

    public final void set_subformat( final JobSubFormat value )
    {
        if ( lock.tryLock() )
        {
            set_int( SUBFORMAT, value.to_ordinal() );
            lock.unlock();
        }
    }
    
    public final void set_progress( final long value, final Boolean save )
    {
        if ( lock.tryLock() )
        {
            current_row = value;
            set_long( PROGRESS, value );
            lock.unlock();
        }
        if ( save ) store();
    }

    public final void set_max( final long value, final Boolean save )
    {
        if ( lock.tryLock() )
        {
            max_row = value;
            set_long( MAX, value );
            lock.unlock();
        }
        if ( save ) store();
    }

    public final void set_md5( final String value, final Boolean save )
    {
        if ( lock.tryLock() )
        {
            set_str_value( MD5, value );
            lock.unlock();
        }
        if ( save ) store();
    }

    private void clear_job()
    {
        set_source( "" );
        set_exportpath( "" );
        set_downloadpath( "" );
        set_state( JobState.INVALID, false );
        set_format( JobFormat.INVALID );
        set_line_ending( LineEndings.AUTOMATIC );
        set_line_wrap( DFLT_LINE_WRAP );
        set_progress( 0, false );
        set_max( 0, false );
        set_runtime( 0 );
        current_row = -1;
        max_row = -1;
    }
    
    public boolean save( final String job_path )
    {
        String filename = String.format( "%s//%s.JOB", job_path, get_short_source() );
        return store_to( filename );
    }

    public boolean change_job_state( final JobState new_state )
    {
        JobState prev_state = get_state();
        boolean res = ( prev_state != new_state );
        if ( res )
        {
            set_state( new_state, true );
            CLogger.logfmt( "job >%s< %s --> %s" ,
                     get_short_source(), prev_state.toString(), new_state.toString() );
        }
        return res;
    }
    
    final void check_paths()
    {
        boolean store_needed = false;
        
        if ( !is_directory( get_exportpath() ) )
        {
            set_exportpath( get_current_dir() );
            store_needed = true;
        }
        
        String p = get_downloadpath();
        if ( p.isEmpty() || !is_directory( get_downloadpath() ) )
        {
            BioAccessionType bt = get_bio_type();
            p = String.format( "%s%sncbi%spublic%s%s",
                    System.getProperty( "user.home" ),
                    sep,
                    sep,
                    sep,
                    bt.equals( BioAccessionType.REF_SEQUENCE ) ? "refseq" : "sra",
                    sep
                    );
            set_downloadpath( p );
            store_needed = true;
        }
        
        if ( store_needed ) store();
    }

    final public boolean does_downloadpath_exist()
    {
        File f = new File( get_downloadpath() );
        return ( f.exists() && f.isDirectory() );
    }

    final public boolean create_downloadpath()
    {
        File f = new File( get_downloadpath() );
        return f.mkdirs();
    }
    
    final public boolean does_exportpath_exist()
    {
        File f = new File( get_exportpath() );
        return ( f.exists() && f.isDirectory() );
    }

    final public boolean create_exportpath()
    {
        File f = new File( get_exportpath() );
        return f.mkdirs();
    }

    final public void set_defaults()
    {
        Settings settings = Settings.getInstance();
        set_exportpath( settings.get_exportpath() );
        set_format( JobFormat.FASTQ );
        set_subformat( JobSubFormat.SPOT );
        set_line_ending( settings.get_line_ending() );
        set_line_wrap( settings.get_line_wrap() );
        set_use_line_wrap( false );
        set_fixed_qual( settings.get_fixed_qual() );
        set_use_fixed_qual( settings.get_use_fixed_qual() );        
        set_max( 0, false ); // 100 k rows... ( because we do not know )
    }
    
    final public void set_defaults( final BioSpec spec )
    {
        set_defaults();
        set_spec( spec );
    }
    
    final void common_init()
    {
        clear();
        clear_job();
        check_paths();
    }
    
    public JobData( final BioSpec spec )
    {
        super( "one job for SRA download" );
        lock = new ReentrantLock();
        state_lock = new ReentrantLock();
        sep = System.getProperty( "file.separator" );
        common_init();        
        set_defaults( spec );
    }
    
    public JobData( final String filename )
    {
        super( "one job for SRA download" );
        lock = new ReentrantLock();
        state_lock = new ReentrantLock();
        sep = System.getProperty( "file.separator" );
        clear();
        clear_job();
        set_valid( load_from( filename ) );
        check_paths();
    }

}
