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
import org.apache.spark.sql.api.java.*;

public enum ngs_field
{
	ID						( DataType.StringType,	EnumSet.of( tt.ALIG, tt.SPOT, tt.FRAG, tt.PILEUP )	),
	REFSPEC					( DataType.StringType,	EnumSet.of( tt.ALIG, tt.REF, tt.PILEUP ) ),
	BASES					( DataType.StringType,	EnumSet.of( tt.ALIG, tt.SPOT, tt.FRAG ) ),
	CIGAR					( DataType.StringType,	EnumSet.of( tt.ALIG ) ),
	CIGARCLIPPED			( DataType.StringType,	EnumSet.of( tt.ALIG ) ),
	CIGARLONG				( DataType.StringType,	EnumSet.of( tt.ALIG ) ),
	CIGARLONGCLIPPED		( DataType.StringType,	EnumSet.of( tt.ALIG ) ),
	REFPOS					( DataType.LongType,	EnumSet.of( tt.ALIG, tt.PILEUP ) ),
	MAPQ					( DataType.IntegerType,	EnumSet.of( tt.ALIG, tt.PILEUP ) ),
	REFBASES				( DataType.StringType,	EnumSet.of( tt.ALIG, tt.PILEUP ) ),
	READGROUP				( DataType.StringType,	EnumSet.of( tt.ALIG, tt.SPOT, tt.FRAG )	),
	READID					( DataType.StringType,	EnumSet.of( tt.ALIG ) ),
	REVERSED				( DataType.BooleanType,	EnumSet.of( tt.ALIG ) ),
	CLIPPEDFRAGMENTBASES	( DataType.StringType,	EnumSet.of( tt.ALIG ) ),
	QUALITIES				( DataType.StringType,	EnumSet.of( tt.ALIG, tt.SPOT, tt.FRAG, tt.PILEUP ) ),
	CATEGORY				( DataType.IntegerType,	EnumSet.of( tt.ALIG, tt.SPOT, tt.FRAG ) ),
	LENGTH					( DataType.LongType,	EnumSet.of( tt.ALIG, tt.REF ) ),
	SOFTCLIPLEFT			( DataType.IntegerType,	EnumSet.of( tt.ALIG ) ),
	SOFTCLIPRIGHT			( DataType.IntegerType,	EnumSet.of( tt.ALIG ) ),
	TEMPLATELENGTH			( DataType.LongType,	EnumSet.of( tt.ALIG ) ),
	RNAORIENTATION			( DataType.StringType,	EnumSet.of( tt.ALIG ) ),
	HASMATE					( DataType.BooleanType,	EnumSet.of( tt.ALIG ) ),
	MATEID					( DataType.StringType,	EnumSet.of( tt.ALIG ) ),
	MATEREFSPEC				( DataType.StringType,	EnumSet.of( tt.ALIG ) ),
	MATEREVERSED			( DataType.BooleanType,	EnumSet.of( tt.ALIG ) ),
	NFRAGMENTS				( DataType.IntegerType,	EnumSet.of( tt.SPOT ) ),
	NAME					( DataType.StringType,	EnumSet.of( tt.SPOT, tt.FRAG, tt.REF, tt.RDGRP ) ),
	CIRCULAR				( DataType.BooleanType,	EnumSet.of( tt.REF ) ),
	NALIGNMENTS				( DataType.LongType,	EnumSet.of( tt.REF ) ),
	DEPTH					( DataType.IntegerType,	EnumSet.of( tt.PILEUP ) ),
	EVENT					( DataType.IntegerType,	EnumSet.of( tt.PILEUP ) ),
	INDEL					( DataType.IntegerType,	EnumSet.of( tt.PILEUP ) ),
	NONE					( DataType.StringType,	EnumSet.noneOf( tt.class ) );

	// private members of each enum-value
	private DataType priv_sql_type;
	private EnumSet< tt > priv_valid_tables;
	
	private ngs_field( DataType sql_type, EnumSet< tt > valid_tables )
	{
		priv_sql_type = sql_type;
		priv_valid_tables = valid_tables;
	}

	public DataType sql_type() { return priv_sql_type; }
	private boolean has_tab( tt t ) { return priv_valid_tables.contains( t ); }
	
	public static ngs_field from_string( final String s )
	{
		try
		{
			return Enum.valueOf( ngs_field.class, s.trim().toUpperCase() );
		}
        catch ( IllegalArgumentException e )
        {
        }
		return NONE;
	}


	private static void make_full_set( EnumSet< ngs_field > res, tt t )
	{
		ngs_field[] all_values = ngs_field.values();
		for ( ngs_field f : all_values )
		{
			if ( f.has_tab( t ) ) res.add( f ); 
		}
	}

	public static EnumSet< ngs_field > make_set( tt t, final String s )
	{
		EnumSet< ngs_field > res = EnumSet.noneOf( ngs_field.class );
		if ( s == null || s.isEmpty() || s.equals( "*" ) )
		{
			make_full_set( res, t );
		}
		else
		{
			String[] splitted = s.split( "," );
			boolean invalid = false;
			for ( String fs : splitted )
			{
				ngs_field f = from_string( fs );
				if ( f != NONE && f.has_tab( t ) )
					res.add( f );
				else
					invalid = true;
			}
			if ( invalid ) res.clear();
		}
		return res;
	}
	
	public static String to_string( EnumSet< ngs_field > fields )
	{
		StringBuffer sb = new StringBuffer();
		boolean c = false;
		for ( ngs_field f : fields )
		{
			if ( c )
				sb.append( "," );
			else
				c = true;
			sb.append( f.toString() );
		}
		return sb.toString();
	}
	
	public boolean can_produce_fastq( EnumSet< ngs_field > fields )
	{
		boolean res = fields.contains( ngs_field.NAME ) || fields.contains( ngs_field.ID );
		res =  res && fields.contains( ngs_field.BASES );
		res =  res && fields.contains( ngs_field.QUALITIES );
		return res;
	}
}
