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

import data.CLogger;
import job.JobData;
import job.JobFormat;
import ngs.ErrorMsg;
import ngs.ReadCollection;
import ngs.ReferenceSequence;

public class BioDumperFactory
{
    private static String NewLine( final JobData job )
    {
        return job.get_line_ending_str();
    }
    
    private static int LineWrap( final JobData job )
    {
        return job.get_use_line_wrap() ? job.get_line_wrap() : 0;
    }
    
    private static BioFormatter make_formatter( final JobData job )
    {
        BioFormatter res = null;
        switch ( job.get_format() )
        {
            case FASTA : res = new BioFastaFormatter( LineWrap( job ), NewLine( job ) ); break;
            case FASTQ : if ( job.get_use_fixed_qual() )
                            res = new BioFastqFormatter( LineWrap( job ), job.get_fixed_qual(), NewLine( job ) );
                         else
                            res = new BioFastqFormatter( LineWrap( job ), NewLine( job ) );
                         break;
        }
        if ( res != null )
            res.set_fastq_dump_style( job.get_fastq_dump_name() );
        return res;
    }

    private static BioDumper make_read_collection_dumper( final JobData job, final long current )
    {
        BioDumper res = null;
        try
        {
            ReadCollection run = gov.nih.nlm.ncbi.ngs.NGS.openReadCollection( job.get_full_source() );
            
            BioFormatter formatter = make_formatter( job );
            if ( formatter != null )
            {
                switch( job.get_subformat() )
                {
                    case SPOT      : res = new BioReadDumper( run, formatter, job, current ); break;
                    case FRAGMENTS : res = new BioFragmentDumper( run, formatter, job, current ); break;
                    case FRAGMENTS_SPLITED : res = new BioFragmentDumper( run, formatter, job, current ); break;
                }
            }
            if ( res != null )
            {
                if ( job.get_use_min_read_len() )
                    res.add_filter( new BioLenFilter( job.get_min_read_len() ) );

                if ( job.get_use_max_N() )
                    res.add_filter( new BioMaxNFilter( job.get_max_N() ) );

                if ( job.get_use_end_N() )
                    res.add_filter( new BioQualEndFilter( job.get_end_N() ) );

                if ( job.get_use_spotgroup() )
                {
                    String spotgroup = job.get_spotgroup();
                    if ( spotgroup.length() > 0  )
                        res.add_filter( new BioSpotGroupFilter( spotgroup ) );
                }
                res.set_read_count_or_size( run.getReadCount() );
            }
        }
        catch ( ErrorMsg ex )
        {
            CLogger.logfmt( "job >%s< error opening read-collection: %s",
                    job.get_short_source(), ex.toString() );
        }
        return res;
    }
    
    private static BioDumper make_reference_dumper( final JobData job, final long current )
    {
        BioDumper res = null;
        try
        {
            ReferenceSequence ref = gov.nih.nlm.ncbi.ngs.NGS.openReferenceSequence( job.get_full_source() );

            BioFormatter formatter = null;
            if ( job.get_format() == JobFormat.FASTA )
            {
                formatter = new BioFastaFormatter( LineWrap( job ), NewLine( job ) );
            }
            if ( formatter != null )
            {
                res = new BioRefDumper( ref, formatter, job, current );
                res.set_read_count_or_size( ref.getLength() );
            }
        }
        catch ( ErrorMsg ex )
        {
            CLogger.logfmt( "job >%s< error opening ReferenceSequence: %s",
                    job.get_short_source(), ex.toString() );
        }
        return res;
    }

    public static BioDumper make_dumper( final JobData job, final long current )
    {
        BioDumper res = null;
        switch( job.get_bio_type() )
        {
            case READ_COLLECTION_UNALIGNED :
            case READ_COLLECTION_ALIGNED   : res = make_read_collection_dumper( job, current ); break;
            case REF_SEQUENCE              : res = make_reference_dumper( job, current ); break;
        }
        return res;
    }

/*    
    public BioDumperFactory()
    {
    }
*/    
}
