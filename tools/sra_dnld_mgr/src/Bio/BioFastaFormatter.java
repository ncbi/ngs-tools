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

public class BioFastaFormatter extends BioFormatter
{
    private static final String HDR_WITH_NAME    = ">%s %s length=%d";
    private static final String HDR_WITHOUT_NAME = ">%s length=%d";
    
    // for READs
    @Override public void format( BioRecord rec, Read obj, final String bases )
    {
        try
        {
            write_format( rec,
                          bases,
                          HDR_WITH_NAME,
                          fastq_dump_style ? reformat_id( obj.getReadId() ) : obj.getReadId(),
                          obj.getReadName(),
                          bases.length() );
            write_newline( rec );
        }
        catch ( ErrorMsg e )
        {
            System.err.println( e.toString() );
        }
    }

    // for FRAGMENTs
    @Override public void format( BioRecord rec, Fragment obj, final String bases )
    {
        try
        {
            write_format( rec,
                          bases,
                          HDR_WITHOUT_NAME,
                          fastq_dump_style ? reformat_id( obj.getFragmentId() ) : obj.getFragmentId(),
                          bases.length() );
            write_newline( rec );
        }
        catch ( ErrorMsg e )
        {
            System.err.println( e.toString() );
        }
    }

    @Override public void format( BioRecord rec, Fragment obj, final String readname, final String bases )
    {
        try
        {
            write_format( rec,
                          bases,
                          HDR_WITH_NAME,
                          fastq_dump_style ? reformat_id( obj.getFragmentId() ) : obj.getFragmentId(),
                          readname,
                          bases.length() );
            write_newline( rec );            
        }
        catch ( ErrorMsg e )
        {
            System.err.println( e.toString() );
        }
    }

    // for ALIGMENTs
    @Override public void format( BioRecord rec, Alignment obj, final String bases )
    {
        try
        {
            write_format( rec,
                          bases,
                          HDR_WITHOUT_NAME,
                          obj.getAlignmentId(),
                          bases.length() );
            write_newline( rec );
        }
        catch ( ErrorMsg e )
        {
            System.err.println( e.toString() );
        }
    }

    // for REFSEQs
    @Override public void format( BioRecord rec, ReferenceSequence obj, final String bases, final long position )
    {
        if ( !header_written )
        {
            try
            {
                write_fmt( rec,
                           HDR_WITHOUT_NAME,
                           obj.getCanonicalName(),
                           obj.getLength() - position );
                write_newline( rec );
                header_written = true;
            }
            catch( ErrorMsg e )
            {
                System.err.println( e.toString() );
            }
        }
        write_globally_clipped( rec, bases );
    }
    
    public BioFastaFormatter( final int line_len, final String line_ending )
    {
        super( line_len, line_ending );
    }
}
