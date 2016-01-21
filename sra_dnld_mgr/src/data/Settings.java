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
package data;

import java.awt.Rectangle;
import java.io.*;

public class Settings extends IniFile
{
    private static final String JOBPATH = "JOBPATH";
    private static final String EXPORTPATH = "EXPORTPATH";
    private static final String MAXDOWNLOADS = "MAXDOWNLOADS";
    private static final String XPOS = "XPOS";
    private static final String YPOS = "YPOS";
    private static final String WIDTH = "WIDTH";
    private static final String HEIGHT = "HEIGHT";
    private static final String AUTOSTART = "AUTOSTART";
    private static final String RESOLVER = "RESOLVER";
    private static final String DNLD_BLOCK_SIZE = "DNLD_BLOCK_SIZE";
    private static final String CONN_TIMEOUT = "CONN_TIMEOUT";
    private static final String READ_TIMEOUT = "READ_TIMEOUT";
    private static final String USER_AGENT = "USER_AGENT";
    private static final String LINE_ENDING = "LINE_ENDING";
    private static final String LINE_WRAP = "LINE_WRAP";
    private static final String USE_LINE_WRAP = "USE_LINE_WRAP";    
    private static final String FIXED_QUAL = "FIXED_QUAL";
    private static final String USE_FIXED_QUAL = "USE_FIXED_QUAL";
    private static final String PREVIEW_ROWS = "PREVIEW_ROWS";
    private static final String LOG_TO_FILE = "LOG_TO_FILE";
    private static final String LOG_TO_CONS = "LOG_TO_CONS";
    
    private static final int DFLT_MAXDOWNLOADS = 2;
    private static final int DFLT_XPOS = 100;
    private static final int DFLT_YPOS = 100;
    private static final int DFLT_WIDTH = 500;
    private static final int DFLT_HEIGHT = 400;
    private static final String DFLT_RESOLVER = "http://www.ncbi.nlm.nih.gov/Traces/names/names.cgi";
    private static final int DFLT_DNLD_BLOCK_SIZE = 1024 * 128;
    private static final int DFLT_CONN_TIMEOUT = 5000;
    private static final int DFLT_READ_TIMEOUT = 5000;
    private static final String DFLT_USER_AGENT = "sra-dnld-mgr";
    private static final LineEndings DFLT_LINE_ENDING = LineEndings.POSIX;
    private static final int DFLT_LINE_WRAP = 75;
    private static final int DFLT_FIXED_QUAL = 30;
    private static final boolean DFLT_USE_LINE_WRAP = true;
    private static final boolean DFLT_USE_FIXED_QUAL = false;
    private static final int DFLT_PREVIEW_ROWS = 10;
    private static final boolean DFLT_LOG_TO_FILE = false;
    private static final boolean DFLT_LOG_TO_CONS = false;
    
    private static final Settings INSTANCE = new Settings();
    public static Settings getInstance() { return INSTANCE; }
    
    private String get_path_of_file( String fn )
    {
        String dir = "";
        File f = new File( fn );
        try
        {
            File fa = new File( f.getCanonicalPath() );
            dir = fa.getParent();
        }
        catch ( Exception e ) { }
        return dir;
    }
    
    private void set_defaults()
    {
        String dir = get_path_of_file( get_filename() );
        set_jobpath( dir );
        set_exportpath( dir );
        set_maxdownloads( DFLT_MAXDOWNLOADS );
        set_resolver( DFLT_RESOLVER );
        set_dnld_block_size( DFLT_DNLD_BLOCK_SIZE );
        set_conn_timeout( DFLT_CONN_TIMEOUT );
        set_read_timeout( DFLT_READ_TIMEOUT );
        set_user_agent( DFLT_USER_AGENT );
        set_line_ending( DFLT_LINE_ENDING );
        set_line_wrap( DFLT_LINE_WRAP );
        set_use_line_wrap( DFLT_USE_LINE_WRAP );
        set_use_fixed_qual( DFLT_USE_FIXED_QUAL );
        set_preview_rows( DFLT_PREVIEW_ROWS );
        set_log_to_file( DFLT_LOG_TO_FILE );
        set_log_to_cons( DFLT_LOG_TO_CONS );
    }
    
    public final String get_jobpath() { return get_str( JOBPATH, "" ); }
    public final String get_exportpath() { return get_str( EXPORTPATH, "" ); }
    public final int get_maxdownloads() { return get_int( MAXDOWNLOADS, DFLT_MAXDOWNLOADS ); }
    public final boolean get_autostart() { return get_bool( AUTOSTART, false ); }
    
    public final int get_xpos()
    {
        int res = get_int( XPOS, DFLT_XPOS );
        if ( res < 0 ) res = 0;
        return res;
    }
    public final int get_ypos()
    {
        int res =  get_int( YPOS, DFLT_YPOS );
        if ( res < 0 ) res = 0;
        return res;
    }
    
    public final int get_width() { return get_int( WIDTH, DFLT_WIDTH ); }
    public final int get_height() { return get_int( HEIGHT, DFLT_HEIGHT ); }
    public final Rectangle get_position()
    {
        return new Rectangle( get_xpos(), get_ypos(), get_width(), get_height() );
    }
    public final String get_resolver() { return get_str( RESOLVER, DFLT_RESOLVER ); }
    public final int get_dnld_block_size() { return get_int( DNLD_BLOCK_SIZE, DFLT_DNLD_BLOCK_SIZE ); }
    public final int get_conn_timeout() { return get_int( CONN_TIMEOUT, DFLT_CONN_TIMEOUT ); }
    public final int get_read_timeout() { return get_int( READ_TIMEOUT, DFLT_READ_TIMEOUT ); }
    public final String get_user_agent() { return get_str( USER_AGENT, DFLT_USER_AGENT ); }
    public final LineEndings get_line_ending()
    {
        return LineEndings.from_ordinal( get_int( LINE_ENDING, DFLT_LINE_ENDING.ordinal() ) );
    }
    public final String get_line_ending_str()
    {
        return LineEndings.from_ordinal( get_int( LINE_ENDING, DFLT_LINE_ENDING.ordinal() ) ).to_line_ending();
    }
    public final int get_line_wrap() { return get_int( LINE_WRAP, DFLT_LINE_WRAP ); }
    public final int get_fixed_qual() { return get_int( FIXED_QUAL, DFLT_FIXED_QUAL ); }
    public final boolean get_use_line_wrap() { return get_bool( USE_LINE_WRAP, DFLT_USE_LINE_WRAP ); }
    public final boolean get_use_fixed_qual() { return get_bool( USE_FIXED_QUAL, DFLT_USE_FIXED_QUAL ); }
    public final int get_preview_rows() { return get_int( PREVIEW_ROWS, DFLT_PREVIEW_ROWS ); }
    public final boolean get_log_to_file() { return get_bool( LOG_TO_FILE, DFLT_LOG_TO_FILE ); }
    public final boolean get_log_to_cons() { return get_bool( LOG_TO_CONS, DFLT_LOG_TO_CONS ); }
    
    public final void set_jobpath( final String value ) { set_str( JOBPATH, value ); }
    public final void set_exportpath( final String value ) { set_str( EXPORTPATH, value ); }
    public final void set_maxdownloads( final int value ) { set_int( MAXDOWNLOADS, value ); }
    public final void set_autostart( final boolean value ) { set_bool( AUTOSTART, value ); }    
    public final void set_xpos( final int value ) { set_int( XPOS, value ); }
    public final void set_ypos( final int value ) { set_int( YPOS, value ); }
    public final void set_width( final int value ) { set_int( WIDTH, value ); }
    public final void set_height( final int value ) { set_int( HEIGHT, value ); }
    public final void set_position( final Rectangle rect, boolean save )
    {
        set_xpos( rect.x );
        set_ypos( rect.y );
        set_width( rect.width );
        set_height( rect.height );
        if ( save ) store();
    }
    public final void set_resolver( final String value ) { set_str( RESOLVER, value ); }
    public final void set_dnld_block_size( final int value ) { set_int( DNLD_BLOCK_SIZE, value ); }
    public final void set_conn_timeout( final int value ) { set_int( CONN_TIMEOUT, value ); }
    public final void set_read_timeout( final int value ) { set_int( READ_TIMEOUT, value ); }
    public final void set_user_agent( final String value ) { set_str( USER_AGENT, value ); }
    public final void set_line_ending( final LineEndings value ) { set_int( LINE_ENDING, value.ordinal() ); }
    public final void set_line_wrap( final int value ) { set_int( LINE_WRAP, value ); }
    public final void set_fixed_qual( final int value ) { set_int( FIXED_QUAL, value ); }
    public final void set_use_line_wrap( final boolean value ) { set_bool( USE_LINE_WRAP, value ); }
    public final void set_use_fixed_qual( final boolean value ) { set_bool( USE_FIXED_QUAL, value ); }
    public final void set_preview_rows( final int value ) { set_int( PREVIEW_ROWS, value ); }
    public final void set_log_to_file( final boolean value ) { set_bool( LOG_TO_FILE, value ); }
    public final void set_log_to_cons( final boolean value ) { set_bool( LOG_TO_CONS, value ); }
    
    final void check_paths()
    {
        boolean store_needed = false;
        
        if ( !is_directory( get_jobpath() ) )
        {
            set_jobpath( get_current_dir() );
            store_needed = true;
        }

        if ( !is_directory( get_exportpath() ) )
        {
            set_exportpath( get_current_dir() );
            store_needed = true;
        }
        if ( store_needed ) store();
    }

    private Settings()
    {
        super( "settings.ini", "global settings for SRA download" );
        if ( !is_valid() )
        {
            set_defaults();
            set_valid( store() );
        }
        check_paths();
    }
    
    public Settings( final String filename )
    {
        super( filename, "global settings for SRA download" );
        if ( !is_valid() )
        {
            set_defaults();
            set_valid( store() );
        }
        check_paths();
    }
}
