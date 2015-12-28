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

/*
    Instances of this class are produced or filled out by the 
    BioRunner - class, put into a queue and consumed by the JobThreadBio - class
*/
public class BioRecord
{
    private final StringBuilder content;
    private int length;
    
    public String get() { return content.toString(); }

    public void set_length( final int value ) { length = value; }
    public int get_length() { return length; }
    public boolean is_empty() { return ( length == 0 ); }
    
    public void append( final String s )
    {
        if ( s.length() > 0 )
            content.append( s );
    }
    
    public void clear()
    {
        content.setLength( 0 );
        length = 0;
    }

    public BioRecord()
    {
        content = new StringBuilder();
        length = 0;
    }
}