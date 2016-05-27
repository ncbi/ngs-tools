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

import java.text.NumberFormat;

public class TimeDiff
{
    private static final long SecondsInMilli = 1000;
    private static final long MinutesInMilli = SecondsInMilli * 60;
    private static final long HoursInMilli = MinutesInMilli * 60;
    private static final long DaysInMilli = HoursInMilli * 24;
    
    static public String to_str( final long milliseconds )
    {
        long days = milliseconds / DaysInMilli;
        long diff = milliseconds % DaysInMilli;
        long hours = diff / HoursInMilli;
        diff = diff % HoursInMilli;
        long minutes = diff / MinutesInMilli;
        diff = diff % MinutesInMilli;
        long seconds = diff / SecondsInMilli;
        if ( days > 0 )
            return String.format( "%d days, %d hours, %d min %d sec", days, hours, minutes, seconds );
        else if ( hours > 0 )
            return String.format( "%d hours, %d min %d sec", hours, minutes, seconds );
        else if ( minutes > 0 )
            return String.format( "%d min %d sec", minutes, seconds );
        else
            return String.format( "%d sec", seconds );
    }
    
    static public long calc_rpm( final long milliseconds, final long reads )
    {
        long res = 0;
        if ( milliseconds > 0  ) res = ( reads * 60000 ) / milliseconds;
        return res;
        //return String.format( "%s rpm", NumberFormat.getInstance().format( rpm ) );
    }
    
    static public String calc_time_left( final long milliseconds, final long reads, final long max )
    {
        long left = 0;
        if ( milliseconds > 0 && reads > 0 )
        {
            left = ( ( milliseconds * max ) / reads ) - milliseconds;
        }
        return to_str( left );
    }
}
