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
package Bio;

import ngs.Read;

public enum BioReadType
{
    INVALID                     ( 0, 0, "INVALID" ),
    READ_TYPE_FULLY_ALIGNED     ( 1, Read.fullyAligned, "fully aligned" ),
    READ_TYPE_PARTIALLY_ALIGNED ( 2, Read.partiallyAligned, "partially aligned" ),
    READ_TYPE_ALIGNED           ( 3, Read.aligned, "aligned" ),
    READ_TYPE_UNALIGNED         ( 4, Read.unaligned, "unaligned" ),
    READ_TYPE_ALL               ( 5, Read.all, "all" );
    
    private final int value;
    private final int read_type;
    private final String txt;
    private static final BioReadType[] allValues = values();

    public static BioReadType from_ordinal( int i )
    {
        return allValues[ i ];
    }
    public int to_ordinal() { return value; }
    public int to_read_type() { return read_type; }
    
    public static String[] to_StringList()
    {
        String[] res = new String[ allValues.length ];
        int i = 0;
        for ( BioReadType rt : allValues ) res[ i++ ] = rt.txt;
        return res;
    }

    public static String[] to_StringList( BioReadType excluded )
    {
        String[] res = new String[ allValues.length - 1 ];
        int i = 0;
        for ( BioReadType rt : allValues )
        {
            if ( !rt.equals( excluded ) )
                res[ i++ ] = rt.txt;
        }
        return res;
    }

    private BioReadType( int value, int read_type, final String txt )
    {
        this.value = value;
        this.read_type = read_type;
        this.txt = txt;
    }
    
}
