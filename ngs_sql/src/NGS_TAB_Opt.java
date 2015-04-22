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

class NGS_TAB_Opt
{
	private String src, tab_name, s_fields;
	private tt tab_type;
	private EnumSet< ngs_field > fields;
	private List< job > jobs;
	private long per_section;
	private int idx;
	
	public NGS_TAB_Opt( final int idx, final String key, final String value )
	{
		this.idx = idx;
		jobs = new LinkedList< job >();
		per_section = 1500;
		fields = null;
		tab_type = tt.SPOT;
		tab_name = key;
		s_fields = "";
		src = "";
		define_line( value );
		if ( jobs.size() < 1 ) jobs.add( new job() );
	}

	private boolean handle_field( String part )
	{
		boolean res = true;
		if ( part.endsWith( "}" ) )
		{
			s_fields = s_fields + part.substring( 0, part.length() - 1 );
			res = false;
		}
		else
		{
			s_fields = s_fields + part;
		}
		return res;
	}
	
	private void define_line( final String value )
	{
		String[] tokens = value.split( " " );
		boolean in_temp_fields = false;
		for ( String s : tokens )
		{
			s = s.trim();
			if ( in_temp_fields )
			{
				in_temp_fields = handle_field( s );
			}
			else if ( s.startsWith( "#" ) )
			{
				jobs.add( new job( s.substring( 1 ) ) );
			}
			else if ( s.startsWith( "{" ) )
			{
				in_temp_fields = handle_field( s.substring( 1 ) );
			}
			else if ( s.startsWith( "%" ) )
			{
				per_section = Long.parseLong( s.substring( 1 ) );
			}
			else
			{
				if ( src.length() < 1 )
				{
					String[] src_tokens = s.split( "\\." );
					if ( src_tokens.length > 0 )
						src = src_tokens[ 0 ];
					if ( src_tokens.length > 1 )
						tab_type = tt.from_string( src_tokens[ 1 ] );
				}
			}
		}
		fields = ngs_field.make_set( tab_type, s_fields );
	}
	
	public String get_src() { return src; }
	public String get_tab_name() { return tab_name; }
	public EnumSet< ngs_field > get_fields() { return fields; }
	public tt get_tab_type() { return tab_type; }
	public long get_per_section() { return per_section; }

	public List< job > get_jobs() { return jobs; }
	
	public void report()
	{
		System.out.printf( "SRC[%d]    =%s%n", idx, src );
		System.out.printf( "TABNAME[%d]=%s%n", idx, tab_name );
		System.out.printf( "TABTYPE[%d]=%s%n", idx, tab_type.name() );
		System.out.printf( "FIELDS[%d] =%s%n", idx, ngs_field.to_string( fields ) );
		System.out.printf( "PER_SEC[%d]=%d%n", idx, per_section );

		System.out.printf( "RANGES[%d] =", idx );
		for ( job j : jobs ) System.out.printf( "%s, ", j.toString() );
		System.out.printf( "%n" );
	}
	
	public boolean valid()
	{
		return ( src != null &&
				  fields != null && fields.size() > 0 &&
				  tab_type != tt.NONE &&
				  jobs != null && jobs.size() > 0 );
	}
}
