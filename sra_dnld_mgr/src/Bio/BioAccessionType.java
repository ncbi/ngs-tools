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

import java.awt.Color;

public enum BioAccessionType
{
    INVALID                     ( 0, Color.red, "INVALID" ),
    READ_COLLECTION_UNALIGNED   ( 1, Color.decode( "0xF5F7BA" ), "READ_COLLECTION_UNALIGNED" ),
    READ_COLLECTION_ALIGNED     ( 2, Color.decode( "0x7DF0E0" ), "READ_COLLECTION_ALIGNED" ),
    REF_SEQUENCE                ( 2, Color.decode( "0xBBE4F2" ), "REF_SEQUENCE" );

    private final int value;
    private final String txt;
    private final Color color;
    private static final BioAccessionType[] allValues = values();

    public static BioAccessionType from_ordinal( int i ) { return allValues[ i ]; }
    public int to_ordinal() { return value; }
    public Color get_color() { return color; }
    
    public static String[] to_StringList()
    {
        String[] res = new String[ allValues.length ];
        int i = 0;
        for ( BioAccessionType bt : allValues ) res[ i++ ] = bt.txt;
        return res;
    }

    public static String[] to_StringList( BioAccessionType excluded )
    {
        String[] res = new String[ allValues.length - 1 ];
        int i = 0;
        for ( BioAccessionType bt : allValues )
        {
            if ( !bt.equals( excluded ) )
                res[ i++ ] = bt.txt;
        }
        return res;
    }
    
    private BioAccessionType( int value, final Color color, final String txt )
    {
        this.value = value;
        this.txt = txt;
        this.color = color;
    }
}
