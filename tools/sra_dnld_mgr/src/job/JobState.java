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

public enum JobState
{
    INVALID (   0,  "INVALID"   ),
    READY   (   1,  "READY"     ),
    RUNNING (   2,  "RUNNING"   ),
    PAUSED  (   3,  "PAUSED"    ),
    ERROR   (   4,  "ERROR"     ),
    DONE    (   5,  "DONE"      ),
    FAILED  (   6,  "FAILED"    );
    
    private final int value;
    private final String txt;
    private static final JobState[] allValues = values();
    
    public static JobState from_ordinal( int n ) { return allValues[ n ]; } 
    public static JobState from_string( String s )
    {
        if ( s.equals( READY.txt ) )
            return READY;
        else if ( s.equals( RUNNING.txt ) )
            return RUNNING;
        else if ( s.equals( PAUSED.txt ) )
            return PAUSED;
        else if ( s.equals( ERROR.txt ) )
            return ERROR;
        else if ( s.equals( DONE.txt ) )
            return DONE;
        else if ( s.equals( FAILED.txt ) )
            return FAILED;
        return INVALID;
    }    
    
    public int to_ordinal() { return value; }
    public String to_string() { return txt; }
    
    private JobState( int value, String txt )
    {
        this.value = value;
        this.txt = txt;
    }
}
