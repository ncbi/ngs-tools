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
import ngs.ErrorMsg;
import ngs.ReferenceSequence;

public class BioRefDumper extends BioDumper
{
    protected final ReferenceSequence ref;
    protected long ref_len;
    private long ref_position;
    private final int ref_chunk = 5000;

    @Override public boolean produce( BioRead read )
    {
        boolean res = false;
        if ( stats != null ) stats.total_objects++;
        
        read.read_id = ref_position;
        
        BioRecord rec = get_empty_bio_record();
        read.put_record( rec );
        try
        {
            String bases = ref.getReferenceBases( ref_position, ref_chunk );
            fmt.format( rec, ref, bases, ref_position );
            rec.set_length( bases.length() );
            res = true;
        }
        catch ( ErrorMsg msg )
        {
            CLogger.logfmt( "job >%s< error calling getReferenceBases(): %s",
                    name, msg.toString() );
        }
        return res;
    }
    
    @Override public boolean next()
    {
        boolean res = false;
        if ( ref != null  )
        {
            ref_position += ref_chunk;
            res = ( ref_position < ref_len );
        }
        return res;
    }

    public BioRefDumper( final ReferenceSequence ref,
                         final BioFormatter fmt,
                         final JobData job,
                         final long start )
    {
        super( fmt, job.get_short_source() );
        this.ref = ref;
        try
        {
            long iter_start = start;
            long iter_count = ref.getLength();

            if ( job.get_use_row_filter() )
            {
                long filter_start = job.get_start_row();
                long filter_count = job.get_row_count();
                long filter_end = filter_start + filter_count - 1;

                if ( iter_start < filter_start )
                    iter_start = filter_start;
                else if ( iter_start > filter_end )
                {
                    CLogger.logfmt( "job >%s< already reached end of filter-row-range", name );
                    iter_count = 0;
                }
                if ( iter_count > 0 )
                    iter_count = filter_end - iter_start + 1;
            }
            
            this.ref_position = iter_start;
            this.ref_len = iter_count;
        }
        catch ( ErrorMsg mgs )
        {
            this.ref_len = 0;
            CLogger.logfmt( "job >%s< error calling getLength(): %s",
                    name, mgs.toString() );
        }
        
    }
    
}
