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

import data.CLogger;
import job.JobData;
import ngs.*;

public class BioFragmentDumper extends BioRunDumper
{
    @Override public boolean produce( BioRead read )
    {
        boolean res = false;
        if ( stats != null ) stats.total_objects++;
        FragmentIterator frag_iter = iter;
        try
        {
            read.parse_read_id( iter.getReadId() );
            while ( frag_iter.nextFragment() )
            {
                BioRecord rec = get_empty_bio_record();
                read.put_record( rec );
                String frag_bases = frag_iter.getFragmentBases();
                if ( pass( iter, frag_iter, frag_bases ) )
                {
                    fmt.format( rec, frag_iter, iter.getReadName(), frag_bases );
                    rec.set_length( frag_bases.length() );
                }
                else
                {
                    if ( stats != null ) stats.number_rejected++;
                }
            }
            res = true;
        }
        catch ( ErrorMsg msg )
        {
            CLogger.logfmt( "job >%s< error calling getFragmentBases(): %s",
                    name, msg.toString() );
        }
        return res;
    }

    public BioFragmentDumper( final ReadCollection run,
                              final BioFormatter fmt,
                              final JobData job,
                              long start ) throws ErrorMsg
    {
        super( run, fmt, job, start );
    }

}
