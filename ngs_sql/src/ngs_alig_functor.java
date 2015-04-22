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
import org.apache.spark.api.java.function.*;
import org.apache.spark.sql.api.java.*;

class ngs_alig_functor extends ngs_functor
{
	public ngs_alig_functor( final String src, EnumSet< ngs_field > fields )
	{
		super( src, fields );
	}
	
	private String alig_read_group( AlignmentIterator iter )
	{
		String res = "";
		try { res = iter.getReadGroup(); }
		catch ( ErrorMsg err ) { ; }
		return res;
	}

	private String alig_rna_orientation( AlignmentIterator iter )
	{
		String res = "";
		try { res = String.valueOf( iter.getRNAOrientation() ); }
		catch ( ErrorMsg err ) { ; }
		return res;
	}
	
	private String alig_mate_id( AlignmentIterator iter )
	{
		String res = "";
		try { res = iter.getMateAlignmentId(); }
		catch ( ErrorMsg err ) { ; }
		return res;
	}

	private String alig_mate_refspec( AlignmentIterator iter )
	{
		String res = "";
		try { res = iter.getMateReferenceSpec(); }
		catch ( ErrorMsg err ) { ; }
		return res;
	}
	
	private boolean alig_mate_reveresed( AlignmentIterator iter )
	{
		boolean res = false;
		try { res = iter.getMateIsReversedOrientation(); }
		catch ( ErrorMsg err ) { ; }
		return res;
	}
	
	public void populate( job job )
	{
		try
		{
			long start = job.get_start();
			long count = job.get_count();
			AlignmentIterator iter = ngs_run.getAlignmentRange( start, count, Alignment.all );
			while ( iter.nextAlignment() )
			{
				Object[] obj = new Object[ fields.size() ];
				int id = 0;
				for ( ngs_field f : fields )
				{
					switch( f )
					{
						case ID			: obj[ id++ ] = iter.getAlignmentId(); break;
						case REFSPEC	: obj[ id++ ] = iter.getReferenceSpec(); break;
						case BASES		: obj[ id++ ] = iter.getAlignedFragmentBases(); break;
						case CIGAR		: obj[ id++ ] = iter.getShortCigar( false ); break;
						case REFPOS		: obj[ id++ ] = iter.getAlignmentPosition(); break;
						case MAPQ		: obj[ id++ ] = iter.getMappingQuality(); break;
						case REFBASES	: obj[ id++ ] = iter.getReferenceBases(); break;
						case READGROUP	: obj[ id++ ] = alig_read_group( iter ); break;
						case READID		: obj[ id++ ] = iter.getReadId(); break;
						case REVERSED	: obj[ id++ ] = iter.getIsReversedOrientation(); break;
						case LENGTH		: obj[ id++ ] = iter.getAlignmentLength(); break;
						case QUALITIES	: obj[ id++ ] = iter.getClippedFragmentQualities(); break;						
						case CLIPPEDFRAGMENTBASES 	: obj[ id++ ] = iter.getClippedFragmentBases(); break;
						case CATEGORY			: obj[ id++ ] = iter.getAlignmentCategory(); break;
						case SOFTCLIPLEFT		: obj[ id++ ] = iter.getSoftClip( Alignment.clipLeft ); break;
						case SOFTCLIPRIGHT		: obj[ id++ ] = iter.getSoftClip( Alignment.clipRight ); break;
						case TEMPLATELENGTH		: obj[ id++ ] = iter.getTemplateLength(); break;
						case CIGARCLIPPED 		: obj[ id++ ] = iter.getShortCigar( true ); break;
						case CIGARLONG			: obj[ id++ ] = iter.getLongCigar( false ); break;
						case CIGARLONGCLIPPED	: obj[ id++ ] = iter.getLongCigar( true ); break;
						case RNAORIENTATION		: obj[ id++ ] = alig_rna_orientation( iter ); break;
						case HASMATE			: obj[ id++ ] = iter.hasMate(); break;
						case MATEID				: obj[ id++ ] = alig_mate_id( iter ); break;
						case MATEREFSPEC		: obj[ id++ ] = alig_mate_refspec( iter ); break;
						case MATEREVERSED		: obj[ id++ ] = alig_mate_reveresed( iter ); break;
						
					}
				}
				iter_res.add( Row.create( obj ) );
			}
		}
		catch ( ErrorMsg err )
		{
			System.out.println( err.toString() );
		}
	}
}