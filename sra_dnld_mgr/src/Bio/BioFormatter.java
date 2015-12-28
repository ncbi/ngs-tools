/*===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
*  Please cite the author in any work or product based on this material.
*
* ===========================================================================
*
*/

package Bio;

import ngs.*;

public class BioFormatter
{
    protected BioStats stats;
    protected final int line_len;
    private int chars_left_on_line;
    protected boolean header_written;
    protected boolean fastq_dump_style;
    public final String NL;

    protected void write( BioRecord rec, final String s )
    {
        rec.append( s );
        if ( stats != null ) stats.total_bytes += s.length();
    }

    protected void write_newline( BioRecord rec )
    {
        write( rec, NL );
    }
    
    protected void write_fmt( BioRecord rec,
                              final String fmt, Object... args )
    {
        write( rec, String.format( fmt, args ) );
    }
    
    protected void write_clipped( BioRecord rec,
                                  final String s )
    {
        int lb = s.length();
        int idx = 0, remaining = lb;
        while( idx < lb )
        {
            if ( remaining >= line_len )
                write( rec, s.substring( idx, idx + line_len ) );
            else
                write( rec, s.substring( idx, idx + remaining ) );
            write_newline( rec );
            
            idx += line_len;
            remaining -= line_len;
        }
    }

    protected void write_globally_clipped( BioRecord rec,
                                           final String bases )
    {
        if ( line_len == 0 )
            write( rec, bases );
        else
        {
            int bases_len = bases.length();
            if ( bases_len < chars_left_on_line )
            {
                write( rec, bases );
                chars_left_on_line -= bases_len;
            }
            else if ( bases_len == chars_left_on_line )
            {
                write( rec, bases );
                write_newline( rec );
                chars_left_on_line = line_len;
            }
            else
            {
                int ofs = 0;
                int remaining = bases_len;
                while( remaining > 0 )
                {
                    if ( remaining >= chars_left_on_line )
                    {
                        write( rec, bases.substring( ofs, ofs + chars_left_on_line ) + '\n' );
                        ofs += chars_left_on_line;
                        remaining -= chars_left_on_line;
                        chars_left_on_line = line_len;                        
                    }
                    else
                    {
                        write( rec, bases.substring( ofs, ofs + remaining ) );
                        chars_left_on_line -= remaining;
                        remaining = 0;
                    }
                }
            }
        }
    }
    
    protected void write_format( BioRecord rec,
                                 final String bases,
                                 final String fmt,
                                 Object... args )
    {
        write_fmt( rec, fmt, args );
        write_newline( rec );
        
        if ( line_len > 0 )
            write_clipped( rec, bases );
        else
            write( rec, bases );
    }
    
    void format ( BioRecord rec, Read obj, final String bases ) {};
    void format ( BioRecord rec, Fragment obj, final String bases ) {};
    void format ( BioRecord rec, Fragment obj, final String readname, final String bases ) {};
    void format ( BioRecord rec, Alignment obj, final String bases ) {};
    void format ( BioRecord rec, ReferenceSequence obj, final String bases, final long position ) {};
    
    void set_stats( BioStats stats )
    {
        this.stats = stats;
    }
    
    public String reformat_id( final String id )
    {
        String[] parts = id.split( "\\." );
        if ( parts.length > 2 )
            return String.format( "%s.%s", parts[ 0 ], parts[ parts.length - 1 ] );
        else
            return id;
    }
            
    public void set_fastq_dump_style( final Boolean value )
    {
        fastq_dump_style = value;
    }
    
    public BioFormatter( int line_len, final String line_ending )
    {
        this.stats = null;
        this.line_len = line_len;
        chars_left_on_line = line_len;
        header_written = false;
        NL = line_ending;
        fastq_dump_style = false;
    }
}
