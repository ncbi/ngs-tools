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
import org.apache.spark.sql.api.java.*;

class ngs_pileup_functor extends ngs_functor
{
	public ngs_pileup_functor( final String src, EnumSet< ngs_field > fields )
	{
		super( src, fields );
	}

	public void populate( job job )
	{
		try
		{
			long ref_pos = 0;
			char ref_char = ' ';
			int depth = 0;
		
			boolean per_base_values, per_event_values;
			
			
			String ref_name = job.get_name();
			Reference ref = ngs_run.getReference( ref_name );
			PileupIterator iter = ref.getPileupSlice( job.get_start(), job.get_count() );

			per_base_values =	fields.contains( ngs_field.REFPOS ) ||
								fields.contains( ngs_field.REFBASES ) ||
								fields.contains( ngs_field.DEPTH );

			per_event_values =	fields.contains( ngs_field.MAPQ ) ||
								fields.contains( ngs_field.ID ) ||
								fields.contains( ngs_field.EVENT ) ||
								fields.contains( ngs_field.QUALITIES ) ||
								fields.contains( ngs_field.INDEL );
								
			while ( iter.nextPileup() )
			{
				if ( per_base_values )
				{
					ref_pos  = iter.getReferencePosition();
					ref_char = iter.getReferenceBase();
					depth = iter.getPileupDepth();
				}
				
				if ( per_event_values )
				{
					while ( iter.nextPileupEvent() )
					{
						Object[] obj = new Object[ fields.size() ];
						int id = 0;
						for ( ngs_field f : fields )
						{
							switch( f )
							{
								case REFSPEC		: obj[ id++ ] = ref_name; break;
								case REFPOS			: obj[ id++ ] = ref_pos; break;
								case REFBASES		: obj[ id++ ] = Character.toString( ref_char ); break;
								case DEPTH			: obj[ id++ ] = depth; break;
								case MAPQ 			: obj[ id++ ] = iter.getMappingQuality(); break;
								case ID				: obj[ id++ ] = iter.getAlignmentId(); break;
								case EVENT			: obj[ id++ ] = iter.getEventType(); break;
								case QUALITIES		: obj[ id++ ] = iter.getAlignmentQuality(); break;
								case INDEL   		: obj[ id++ ] = iter.getEventIndelType(); break;
							}
						}
						iter_res.add( Row.create( obj ) );
					}
				}
				else
				{
					Object[] obj = new Object[ fields.size() ];
					int id = 0;
					for ( ngs_field f : fields )
					{
						switch( f )
						{
							case REFSPEC		: obj[ id++ ] = ref_name; break;
							case REFPOS			: obj[ id++ ] = ref_pos; break;
							case REFBASES		: obj[ id++ ] = Character.toString( ref_char ); break;
							case DEPTH			: obj[ id++ ] = depth; break;
						}
					}
					iter_res.add( Row.create( obj ) );
				}
			}
		}
		catch ( ErrorMsg err )
		{
			System.out.println( err.toString() );
		}
	}
}