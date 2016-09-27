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

import job.JobState;
import static job.JobState.DONE;
import static job.JobState.ERROR;
import static job.JobState.INVALID;
import static job.JobState.PAUSED;
import static job.JobState.READY;
import static job.JobState.RUNNING;

public enum LineEndings
{
    INVALID             ( 0,    "INVALID",      "\n" ),
    AUTOMATIC           ( 1,    "AUTOMATIC",    "\n" ),
    POSIX               ( 2,    "POSIX",        "\n" ),
    WINDOWS             ( 3,    "WINDOWS",      "\r\n" );

    private final int value;
    private final String txt;
    private final String ending;
    private static final LineEndings[] allValues = values();

    public static LineEndings from_ordinal( int i ) { return allValues[ i ]; }
    public static LineEndings from_string( String s )
    {
        if ( s.equals( AUTOMATIC.txt ) )
            return AUTOMATIC;
        else if ( s.equals( POSIX.txt ) )
            return POSIX;
        else if ( s.equals( WINDOWS.txt ) )
            return WINDOWS;
        return INVALID;
    }    

    public int to_ordinal() { return value; }
    public String to_string() { return txt; }
    
    public String to_line_ending()
    {
        if ( value == 1 )
            return System.getProperty( "line.separator" );
        else
            return ending;
    }

    public static String[] to_StringList()
    {
        String[] res = new String[ allValues.length ];
        int i = 0;
        for ( LineEndings le : allValues ) res[ i++ ] = le.txt;
        return res;
    }

    public static String[] to_StringList( LineEndings excluded )
    {
        String[] res = new String[ allValues.length - 1 ];
        int i = 0;
        for ( LineEndings le : allValues )
        {
            if ( !le.equals( excluded ) )
                res[ i++ ] = le.txt;
        }
        return res;
    }
    
    private LineEndings( int value, final String txt, final String ending )
    {
        this.value = value;
        this.txt = txt;
        this.ending = ending;
    }
}
