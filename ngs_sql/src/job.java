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

public class job implements java.io.Serializable
{
    private final long start;
    private final long count;
	private final String name;

	public job()
    {
        this.start = 1;
        this.count = 0;
		this.name = "";
    }

    public job( final long start )
    {
        this.start = start;
        this.count = 0;
		this.name = "";
    }
	
    public job( final long start, final long count )
    {
        this.start = start;
        this.count = count;
		this.name = "";
    }

    public job( final long start, final long count, final String name )
    {
        this.start = start;
        this.count = count;
		this.name = name;
    }

	public job( final String s )
	{
		String s_name = "";
		String s_start = "1";
		String s_count = "0";

		String[] tokens1 = s.split( "\\:" );
		if ( tokens1.length == 1 )
		{
			if ( tokens1[ 0 ].startsWith( "@" ) )
			{
				s_name = tokens1[ 0 ].substring( 1 );
			}
			else
			{
				String[] tokens2 = tokens1[ 0 ].split( "\\." );
				if ( tokens2.length > 0 )
					s_start = tokens2[ 0 ];
				if ( tokens2.length > 1 )
					s_count = tokens2[ 1 ];
			}
		}
		else if ( tokens1.length > 1 )
		{
			if ( tokens1[ 0 ].startsWith( "@" ) )
			{
				s_name = tokens1[ 0 ].substring( 1 );
				String[] tokens2 = tokens1[ 1 ].split( "\\." );
				if ( tokens2.length > 0 )
					s_start = tokens2[ 0 ];
				if ( tokens2.length > 1 )
					s_count = tokens2[ 1 ];
			}
			else
			{
				String[] tokens2 = tokens1[ 0 ].split( "\\." );
				if ( tokens2.length > 0 )
					s_start = tokens2[ 0 ];
				if ( tokens2.length > 1 )
					s_count = tokens2[ 1 ];
			}
		}
		
		start = Long.parseLong( s_start );
		count = Long.parseLong( s_count );
		name = s_name;
	}

	public long get_start() { return start; }
	public long get_count() { return count; }
	public String get_name()  { return name; }
	public boolean has_name() { return !name.equals( "" ); }

    @Override public String toString()
    {
        StringBuffer sb = new StringBuffer();
		if ( name.length() > 0 )
		{
			sb.append( "@" );
			sb.append( name );
			sb.append( ":" );
		}
		sb.append( start );
		sb.append( "." );
		sb.append( count );
        return sb.toString();
    }
}
