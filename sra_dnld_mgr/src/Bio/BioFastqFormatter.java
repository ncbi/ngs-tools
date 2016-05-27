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

public class BioFastqFormatter extends BioFormatter
{
    private static final String READ_HDR_1   = "@%s %s length=%d";
    private static final String READ_HDR_2   = "+%s %s length=%d";
    private static final String FRAG_HDR_1   = "@%s length=%d";
    private static final String FRAG_HDR_2   = "+%s length=%d";    
    private static final String ALIG_HDR_1   = "@%s length=%d";
    private static final String ALIG_HDR_2   = "+%s length=%d";
    private static final String REFSEQ_HDR_1 = "@%s length=%d";
    private static final String REFSEQ_HDR_2 = "+%s length=%d";
    
    private String qual;
    private final StringBuilder qual_builder;
    private final int fixed_qual;
    private final char fixed_qual_char;
    
    private void make_fixed_qual( int len )
    {
        if ( qual_builder.length() >= len )
            qual = qual_builder.substring( 0, len );
        else
        {
            while ( qual_builder.length() < len )
                qual_builder.append( fixed_qual_char);
            qual = qual_builder.toString();
        }
    }
    
    // for READs
    @Override public void format( BioRecord rec, Read obj, final String bases )
    {
        try
        {
            int len = bases.length();
            String id = fastq_dump_style ? reformat_id( obj.getReadId() ) : obj.getReadId();
            write_format( rec,
                          bases,
                          READ_HDR_1,
                          id,
                          obj.getReadName(),
                          len );
            write_newline( rec );
            
            if ( fixed_qual == 0 )
                qual = obj.getReadQualities();
            else
                make_fixed_qual( len );

            write_format( rec,
                          qual,
                          READ_HDR_2,
                          id,
                          obj.getReadName(),
                          len );
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
            int len = bases.length();
            if ( len > 0 )
            {
                String id = fastq_dump_style ? reformat_id( obj.getFragmentId() ) : obj.getFragmentId();
                write_format( rec,
                              bases,
                              FRAG_HDR_1,
                              id,
                              len );
                write_newline( rec );

                if ( fixed_qual == 0 )
                    qual = obj.getFragmentQualities();
                else
                    make_fixed_qual( len );
                
                write_format( rec,
                              qual,
                              FRAG_HDR_2,
                              id,
                              len );
                write_newline( rec );
            }
        }
        catch ( ErrorMsg e )
        {
            System.err.println( e.toString() );
        }
    }

    // for FRAGMENTs
    @Override public void format( BioRecord rec, Fragment obj, final String readname, final String bases )
    {
        try
        {
            int len = bases.length();
            if ( len > 0 )
            {
                String id = fastq_dump_style ? reformat_id( obj.getFragmentId() ) : obj.getFragmentId();
                write_format( rec,
                              bases,
                              READ_HDR_1,
                              id,
                              readname,
                              len );
                write_newline( rec );

                if ( fixed_qual == 0 )
                    qual = obj.getFragmentQualities();
                else
                    make_fixed_qual( len );
                
                write_format( rec,
                              qual,
                              READ_HDR_2,
                              id,
                              readname,
                              len );
                write_newline( rec );
            }
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
            int len = bases.length();
            if ( len > 0 )
            {
                String id = fastq_dump_style ? reformat_id( obj.getAlignmentId() ) : obj.getAlignmentId();
                write_format( rec,
                              bases,
                              ALIG_HDR_1,
                              id,
                              len );
                write_newline( rec );
                
                if ( fixed_qual == 0 )
                    qual = obj.getClippedFragmentQualities();
                else
                    make_fixed_qual( len );
                
                write_format( rec,
                              qual,
                              ALIG_HDR_2,
                              id,
                              len );
                write_newline( rec );
            }
        }
        catch ( ErrorMsg e )
        {
            System.err.println( e.toString() );
        }
    }

    // for REFSEQs
    @Override public void format( BioRecord rec, ReferenceSequence obj,
                                  String bases,
                                  long position )
    {
        if ( !header_written )
        {
            try
            {
                write_fmt( rec, REFSEQ_HDR_1,
                        obj.getCanonicalName(), obj.getLength() - position );
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
    
    public BioFastqFormatter( final int line_len, final String line_ending )
    {
        super( line_len, line_ending );
        qual = "";
        this.fixed_qual = 0;
        fixed_qual_char = '?';
        qual_builder = new StringBuilder();
    }

    public BioFastqFormatter( final int line_len, final int fixed_qual, final String line_ending )
    {
        super( line_len, line_ending );
        qual = "";
        this.fixed_qual = fixed_qual;
        fixed_qual_char = Character.toChars( fixed_qual + 33 )[ 0 ];
        qual_builder = new StringBuilder();
    }

}
