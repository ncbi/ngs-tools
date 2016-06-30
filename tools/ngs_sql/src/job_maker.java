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

import java.util.*;
import ngs.*;

class job_maker
{
	private NGS_TAB_Opt opt;
	ReadCollection ngs_run;
	
	public job_maker( NGS_TAB_Opt opt )
	{
		this.opt = opt;
		try
		{
			ngs_run = gov.nih.nlm.ncbi.ngs.NGS.openReadCollection( opt.get_src() );
        }
        catch ( ErrorMsg err )
        {
			System.err.println( err.toString() );
			ngs_run = null;
		}
	}

	private long get_count( long requested_count, long run_count )
	{
		long res = run_count;
		if ( requested_count > 0 && requested_count < res ) res = requested_count;
		return res;
	}
	
	private boolean make_spot_alig_jobs_from_requested_range( job src, List< job > jobs, long run_count )
	{
		long count = get_count( src.get_count(), run_count );
		if ( count < 1 ) return false;

		long per_section = opt.get_per_section();
		if ( per_section < 1 ) return false;

		long start = src.get_start();
		while ( count > 0 )
		{
			if ( count < per_section ) per_section = count;
			jobs.add( new job( start, per_section ) );
			start += per_section;
			count -= per_section;
		}
		return true;
	}
	
	private boolean make_spot_alig_jobs( List< job > jobs, long run_count )
	{
		boolean res = true;
		for ( job j : opt.get_jobs() )
		{
			if ( res )
				res = res && make_spot_alig_jobs_from_requested_range( j, jobs, run_count );
		}
		return res;
	}

	private boolean make_pileup_jobs_from_ref_iter( Reference ref, job src, List< job > jobs )
	{
		boolean res = false;
		try
		{
			String ref_name = ref.getCanonicalName();
			long count = get_count( src.get_count(), ref.getLength() );
			res = ( count > 0 );
			if ( res )
			{
				long per_section = opt.get_per_section();
				res = ( per_section > 0 );
				if ( res )
				{
					long start = src.get_start();
					while ( count > 0 )
					{
						if ( count < per_section ) per_section = count;
						jobs.add( new job( start, per_section, ref_name ) );
						start += per_section;
						count -= per_section;
					}
				}
			}
		}
		catch ( ErrorMsg err )
		{
			System.err.println( err.toString() );
		}
		return res;
	}
	
	private boolean make_pileup_jobs_from_requested_range( job src, List< job > jobs )
	{
		boolean res = false;
		if ( src.has_name() )
		{
			try
			{
				Reference ref = ngs_run.getReference ( src.get_name() );
				res = make_pileup_jobs_from_ref_iter( ref, src, jobs );
			}
			catch ( ErrorMsg err )
			{
				System.err.println( err.toString() );
			}
		}
		else
		{
			try
			{
				ReferenceIterator ref_iter = ngs_run.getReferences();
				res = true;
				while ( ref_iter.nextReference() && res )
				{
					res = res && make_pileup_jobs_from_ref_iter( ref_iter, src, jobs );
				}
			}
			catch ( ErrorMsg err )
			{
				System.err.println( err.toString() );
			}
		}
		return res;
	}
	
	private boolean make_pileup_jobs( List< job > jobs )
	{
		boolean res = true;
		for ( job j : opt.get_jobs() )
		{
			if ( res )
				res = res && make_pileup_jobs_from_requested_range( j, jobs );
		}
		return res;
	}
	
    private long alig_count()
    {
		if ( ngs_run == null ) return 0;
		try { return ngs_run.getAlignmentCount(); }
        catch ( ErrorMsg err ) { System.err.println( err.toString() ); }
		return 0;
    }

    private long spot_count()
    {
		if ( ngs_run == null ) return 0;	
		try	{ return ngs_run.getReadCount(); }
        catch ( ErrorMsg err ) { System.err.println( err.toString() ); }
		return 0;
    }
	
	public boolean make_jobs( List< job > jobs )
	{
		boolean res = false;
		switch( opt.get_tab_type() )
		{
			case SPOT  	: res = make_spot_alig_jobs( jobs, spot_count() ); break;		
			case FRAG  	: res = make_spot_alig_jobs( jobs, spot_count() ); break;					
			case ALIG  	: res = make_spot_alig_jobs( jobs, alig_count() ); break;
			case REF   	: jobs.add( new job( 1, 1 ) ); res = true; break;
			case RDGRP 	: jobs.add( new job( 1, 1 ) ); res = true; break;
			case PILEUP	: res = make_pileup_jobs( jobs ); break;
		}
		return res;
	}
}
