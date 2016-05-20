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

enum res_type
{
	STRING, LONG, INT, BYTE, DOUBLE, FLOAT, BOOL, NONE;
}

class res_line_functor implements Function < Row, String >
{
	private StringBuffer sb;
	private List< res_type > typelist;
	private String delim;
	
	public res_line_functor( StructField schema_fields[], String delim )
	{
		sb = new StringBuffer();
		typelist = new ArrayList< res_type >();
		this.delim = delim;
		for ( StructField sf : schema_fields )
		{
			Class c = sf.getDataType().getClass();
			res_type rt = res_type.NONE;
			if ( c == StringType.class )			{ rt = res_type.STRING; }
			else if ( c == IntegerType.class )		{ rt = res_type.INT; }
			else if ( c == LongType.class )		{ rt = res_type.LONG; }
			else if ( c == BooleanType.class )		{ rt = res_type.BOOL; }
			else if ( c == ByteType.class )		{ rt = res_type.BYTE; }
			else if ( c == DoubleType.class )		{ rt = res_type.DOUBLE; }
			else if ( c == FloatType.class )		{ rt = res_type.FLOAT; }
			typelist.add( rt );
		}
	}
	
	public String call( Row row )
	{
		sb.setLength( 0 );
		boolean c = false;
		int i, n = row.length();
		for ( i = 0; i < n; ++i )
		{
			if ( c )
				sb.append( delim );
			else
				c = true;

			switch( typelist.get( i ) )
			{
				case STRING		: sb.append("\""); sb.append( row.getString( i ) ); sb.append("\""); break;
				case LONG		: sb.append( row.getLong( i ) ); break;
				case INT		: sb.append( row.getInt( i ) ); break;
				case BOOL		: sb.append( row.getBoolean( i ) ); break;
				case BYTE		: sb.append( row.getByte( i ) ); break;
				case DOUBLE		: sb.append( row.getDouble( i ) ); break;
				case FLOAT		: sb.append( row.getFloat( i ) ); break;
			}
		}
		return sb.toString();
	}
}
