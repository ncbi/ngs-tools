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
import org.apache.spark.api.java.*;
import org.apache.spark.api.java.function.*;
import org.apache.spark.sql.api.java.*;


class fastq_functor implements FlatMapFunction < Row, String >
{
	private List< String > res;
	private int name_idx, bases_idx, quality_idx;
	private boolean valid;
	
	public fastq_functor( StructField schema_fields[] )
	{
		res = new LinkedList< String >();
		valid = get_name_base_qual_idx( schema_fields );
	}
	
	private boolean get_name_base_qual_idx( StructField schema_fields[] )
	{
		int i = 0;
		int idx_id = -1;
		int idx_name = -1;
		for ( StructField sf : schema_fields )
		{
			ngs_field f = ngs_field.from_string( sf.getName().toUpperCase() );
			switch( f )
			{
				case BASES		: bases_idx = i; break;
				case QUALITIES	: quality_idx = i; break;
				case ID			: idx_id = i; break;
				case NAME		: idx_name = i; break;
			}
			i++;
		}
		if ( idx_id > -1 )
		{
			name_idx = idx_id;
		}
		else
		{
			if ( idx_name > -1 )
				name_idx = idx_name;
		}
		return ( name_idx > -1 &&  bases_idx > -1 && quality_idx > -1 );
	}
	
	public Iterable< String > call( Row row )
	{
		res.clear();
		if ( valid )
		{
			res.add( row.getString( name_idx ) );
			res.add( row.getString( bases_idx ) );
			res.add( "+" );
			res.add( row.getString( quality_idx ) );
		}
		return res;
	}

}
