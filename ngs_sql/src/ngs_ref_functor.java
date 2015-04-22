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

class ngs_ref_functor extends ngs_functor
{
	public ngs_ref_functor( final String src, EnumSet< ngs_field > fields )
	{
		super( src, fields );
	}

	public void populate( job job )
	{
		try
		{
			ReferenceIterator iter = ngs_run.getReferences();
			while ( iter.nextReference() )
			{
				Object[] obj = new Object[ fields.size() ];
				int id = 0;
				for ( ngs_field f : fields )
				{
					switch( f )
					{
						case NAME			: obj[ id++ ] = iter.getCommonName(); break;
						case REFSPEC		: obj[ id++ ] = iter.getCanonicalName(); break;
						case LENGTH			: obj[ id++ ] = iter.getLength(); break;
						case CIRCULAR		: obj[ id++ ] = iter.getIsCircular(); break;
						case NALIGNMENTS	: obj[ id++ ] = iter.getAlignmentCount(); break;
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