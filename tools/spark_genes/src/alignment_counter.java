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
import org.apache.spark.api.java.function.PairFunction;
import ngs.*;
import scala.Tuple2;

public class alignment_counter implements PairFunction < gtf_feature, String, Integer >
{
	private String accession;
	private String curr_ref_name;
	private ReadCollection run;
	private Reference curr_ref;
	private int min_mapq;
	private count_mode cm;
    private boolean initialized;
    private interval_list_cmp cmp;
    private interval_list feature_intervals;
    private interval_list read_intervals;
	
	public alignment_counter( String accession, int min_mapq, count_mode cm )
	{
        this.initialized = false;
		this.accession = accession;
		this.curr_ref_name = "";
		this.curr_ref = null;
		this.run = null;
		this.min_mapq = min_mapq;
		this.cm = cm;
        this.cmp = null;
        this.read_intervals = null;
        this.feature_intervals = null;
	}

    private boolean open_accession()
    {
        boolean res = false;
        try
        {
            run = gov.nih.nlm.ncbi.ngs.NGS.openReadCollection( accession );
            res = true;
        }
        catch ( ErrorMsg err )
        {
            System.out.println( err.toString() );
        }
        return res;
    }

    private void open_reference( gtf_feature f )
    {
        curr_ref_name = f.get_ref();
        try
        {
            curr_ref = run.getReference( curr_ref_name );
        }
        catch ( ErrorMsg err )
        {
            System.out.println( err.toString() );
        }
    }

	public Tuple2< String, Integer > call( gtf_feature f )
	{
        Integer count = 0;

        if ( !initialized )
        {
            cmp = new interval_list_cmp( false );
            read_intervals = new interval_list( 512 );
            feature_intervals = new interval_list( 512 );
            initialized = open_accession();
        }

		if ( initialized )
		{
			/* here we have to make a ngs-iter for the gtf_feature comming in via element-tuple */
			if ( !f.has_ref( curr_ref_name ) )
            {
                curr_ref = null;
                open_reference( f );
            }

			if ( curr_ref != null )
			{
				try 
				{
                    feature_intervals.clear();
                    int n = f.get_ranges();
                    for ( int i = 0; i < n; ++i )
                        feature_intervals.add( f.get_start_at( i ), f.get_len_at( i ), f.get_amb_at( i ) );

					AlignmentIterator al_iter = curr_ref.getAlignmentSlice ( f.get_start(), f.get_len() );
					while ( al_iter.nextAlignment() )
					{
                        boolean alignment_reverse = al_iter.getIsReversedOrientation();
                        if ( alignment_reverse != f.is_reverse() )
                        {
                            //tc.on_other_strand();
                        }
                        else
                        {
                            if ( al_iter.getMappingQuality() < min_mapq )
                            {
                                //tc.rejected_mapq();
                            }
                            else
                            {
                                read_intervals.clear();
                                read_intervals.from_cigar( al_iter.getAlignmentPosition(), al_iter.getShortCigar( false ) );
                                count_result cr = cmp.walk( feature_intervals, read_intervals, cm );
                                //tc.count( cr );
                                if ( cr == count_result.FEATURE )
                                    count++;
                            }
                        }
					}
				}
				catch( ErrorMsg err )
				{
					System.out.println( err.toString() );
				}
			}
		}
		return new Tuple2<>( f.get_id(), count );	
	}
}
