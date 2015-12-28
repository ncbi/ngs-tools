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

import java.util.*;
import java.io.*;

public class IniFile
{
    private String filename;
    private final String hint;
    private final Properties props;
    private boolean valid;

    private int to_int( String s, int dflt )
    {
        int res;
        try { res = Integer.parseInt( s ); }
        catch ( Exception e ) { res = dflt; }
        return res;
    }

    private long to_long( String s, long dflt )
    {
        long res;
        try { res = Long.parseLong( s ); }
        catch ( Exception e ) { res = dflt; }
        return res;
    }

    private boolean to_bool( String s, boolean dflt )
    {
        boolean res;
        try { res = Boolean.parseBoolean( s ); }
        catch ( Exception e ) { res = dflt; }
        return res;
    }

    public final String get_str( String key, String dflt )
    {
        String res = props.getProperty( key );
        if ( res == null )
            res = dflt;
        else if ( res.isEmpty() )
            res = dflt;
        return res;
    }
    
    public final int get_int( String key, int dflt ) { return to_int( props.getProperty( key ), dflt ); }
    public final long get_long( String key, long dflt ) { return to_long( props.getProperty( key ), dflt ); }
    public final boolean get_bool( String key, boolean dflt ) { return to_bool( props.getProperty( key ), dflt ); }
    
    public final void set_str( String key, String value ) { props.setProperty( key, value ); }
    public final void set_int( String key, int value ) { props.setProperty( key, Integer.toString( value ) ); }
    public final void set_long( String key, long value ) { props.setProperty( key, Long.toString( value ) ); }
    public final void set_bool( String key, boolean value ) { props.setProperty( key, Boolean.toString( value ) ); }
    
    public void clear()
    {
        filename = "";
        props.clear();
    }
    
    public final boolean store()
    {
        boolean res = false;
        if ( !filename.isEmpty() )
        {
            try
            {
                FileOutputStream s = new FileOutputStream( filename );
                props.store( s, hint );
                s.close();
                res = true;
            }
            catch ( Exception e ) { }
        }
        return res;
    }

    public final boolean store_to( String new_filename )
    {
        boolean res = false;
        if ( !new_filename.isEmpty() )
        {
            filename = new_filename;
            res = store();
        }
        return res;
    }
            
    public final boolean load()
    {
        boolean res = false;
        if ( !filename.isEmpty() )
        {
            try
            {
                FileInputStream s = new FileInputStream( filename );
                props.load( s );
                s.close();
                res = true;
            }
            catch ( Exception e ) { }
        }
        return res;
    }
    
    public final boolean load_from( String new_filename )
    {
        boolean res = false;
        if ( !new_filename.isEmpty() )
        {
            filename = new_filename;
            res = load();
        }
        return res;
    }
    
    public boolean is_valid() { return valid; }
    public void set_valid( boolean value ) { valid = value; }
    public String get_filename() { return filename; }

    public boolean is_directory( final String path )
    {
        File f = new File( path );
        if ( f.exists() && f.isDirectory() ) return true;
        return false;
    }
    
    public String get_current_dir()
    {
        String res;
        File f = new File( "." );
        try
        {
            res = f.getCanonicalPath();
        }
        catch ( IOException ex )
        {
            res = f.getAbsolutePath();
        }
        return res;
    }
    
    public boolean delete_file()
    {
        boolean res = false;
        if ( !filename.isEmpty() )
        {
            File f = new File( filename );
            res = f.delete();
        }
        return res;
    }
    
    public IniFile( String filename, String hint )
    {
        this.filename = filename;
        this.hint = hint;
        props = new Properties();
        valid = load();
    }

    public IniFile( String hint )
    {
        this.filename = "";
        this.hint = hint;
        props = new Properties();
        valid = false;
    }
    
}
