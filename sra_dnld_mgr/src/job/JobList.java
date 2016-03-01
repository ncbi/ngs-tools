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

import java.io.*;

class ExtensionFilter implements FilenameFilter
{
    private final String extension;
    
    public ExtensionFilter( String extension )
    { this.extension = extension; }
    
    @Override public boolean accept( File dir, String name )
    {
        boolean res = false;
        int last_dot = name.lastIndexOf( '.' );
        if ( last_dot > 0 )
        {
            String ext = name.substring( last_dot );
            res = ( ext.equals( extension ) );
        }
        return res;
    }
}

public class JobList
{
    private final String path;
    private final File dir;
    private final ExtensionFilter job_filter;
    
    public boolean contains_job( String jobname )
    {
        boolean res = false;
        if ( dir.isDirectory() )
        {
            String[] entries = dir.list( job_filter );
            for ( String s : entries )
                res = res | ( s.equals( jobname ) );
        }
        return res;
    }

    public String[] make_list()
    {
        if ( dir.isDirectory() )
            return dir.list( job_filter );
        else
            return null;
    }
    
    public JobList( String path )
    {
        this.path = path;
        dir = new File( path );
        job_filter = new  ExtensionFilter( ".JOB" );
    }
}
