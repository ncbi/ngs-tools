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

class interval_list_storage
{
	private int size, members;
    private long v_start[];
    private int v_len[];
	private boolean v_amb[];

	public void interval_list_storage( final int n )
	{
		size = n;
		members = 0;
		v_start = new long[ size ];
		v_len   = new int[ size ];
		v_amb   = new boolean[ size ];
	}
	
	public boolean idx_valid( final int idx ) { return ( idx >= 0 && idx < members ); }
	
	public long get_start( final int idx ) { if ( idx_valid( idx ) ) return v_start[ idx ]; else return 0; }
	public int get_len( final int idx ) { if ( idx_valid( idx ) ) return v_len[ idx ]; else return 0; }
	public long get_end( final int idx ) { if ( idx_valid( idx ) ) return v_start[ idx ] + v_len[ idx ] - 1; else return 0; }
	public boolean is_amb( final int idx ) { if ( idx_valid( idx ) ) return v_amb[ idx ]; else return false; }

	public void set_start( final int idx, final long value ) { if ( idx_valid( idx ) ) v_start[ idx ] = value; }
	public void set_len( final int idx, final int value ) { if ( idx_valid( idx ) && value >= 0 ) v_len[ idx ] = value; }
	public void set_end( final int idx, final long value ) { set_len( idx, (int)( value - get_start( idx ) + 1 ) ); }
	public void set_amb( final int idx, final boolean value ) { if ( idx_valid( idx ) ) v_amb[ idx ] = value; }

	public void set( final int idx, final long start, final int len, final boolean amb )
	{
		set_start( idx, start );
		set_len( idx, len );
		set_amb( idx, amb );
	}

	public int get_count() { return members; }
	public boolean is_empty() { return ( members == 0 ); }
	public void set_count( final int value ) { if ( value <= size && value >= 0 ) members = value; }
	public void clear() { set_count( 0 ); }
	
	public int get_size() { return size; }
	
	public void add( final long start, final int len, final boolean amb )
	{
		if ( members < size ) set( members++, start, len, amb );
	}

	public void remove( final int at )
	{
		if ( idx_valid( at ) )
		{
			int n = get_count() - 1;
			for ( int i = at; i < n; i++ )
				set( i, get_start( i + 1 ), get_len( i + 1 ), is_amb( i + 1 ) );
			set_count( n );
		}
	}

	public void insert( final int at, final long start, final int len, final boolean amb )
	{
		if ( idx_valid( at ) )
		{
			int n = get_count();
			if ( n < size )
			{
				set_count( n + 1 );
				for ( int i = n; i > at; i-- )
					set( i, get_start( i - 1 ), get_len( i - 1 ), is_amb( i - 1 ) );
				set( at, start, len, amb );
			}
		}
		else
			add( start, len, amb );
	}
	
	public boolean equals( final interval_list_storage other )
	{
		int n = get_count();
		boolean res = ( n == other.get_count() );
		for ( int i = 0; res && i < n; i++ )
		{
			res = ( get_start( i ) == other.get_start( i ) &&
					get_len( i ) == other.get_len( i ) &&
					is_amb( i ) == other.is_amb( i ) );
		}
		return res;
	}

    public void toStartLen( StringBuffer sb )
    {
		int n = get_count();
		for( int i = 0; i < n; i++ )
		{
            if ( i > 0 ) sb.append( ";" );
			if ( is_amb( i ) ) sb.append( "A" ); else sb.append( "N" );
            sb.append( get_start( i ) );
            sb.append( "." );
            sb.append( get_len( i ) );
		}
    }

	@Override public String toString()
	{
		StringBuffer sb = new StringBuffer();
		int n = get_count();
		for( int i = 0; i < n; i++ )
		{
            sb.append( "[" );
            if ( is_amb( i ) ) sb.append( "A" ); else sb.append( " " );
            sb.append( get_start( i ) );
            sb.append( "." );
            sb.append( get_len( i ) );
            sb.append( " .. " );
            sb.append( get_end( i ) );
            sb.append( "]" );
		}
		return sb.toString();
	}

}


class sortable_interval_list extends interval_list_storage
{
	private boolean sorted;
	
	public void sortable_interval_list( final int n )
	{
		interval_list_storage( n );
		sorted = true;
	}
	
	private void swap( final int idx1, final int idx2 )
	{
		long temp_start = get_start( idx1 );
		int temp_len = get_len( idx1 );
		boolean temp_amb = is_amb( idx1 );
		
		set( idx1, get_start( idx2 ), get_len( idx2 ), is_amb( idx2 ) );
		set( idx2, temp_start, temp_len, temp_amb );
	}

	public boolean is_sorted() { return sorted; }
	
	public void add( final long start, final int len, final boolean amb )
	{
		if ( is_empty() )
		{
			super.add( start, len, amb );
			sorted = true;
		}
		else if ( sorted )
		{
			long last_start = get_start( get_count() - 1 );
			super.add( start, len, amb );
			sorted = ( last_start <= start );
		}
	}
	
	public void sort()
	{
		if ( !sorted )
		{
			int n = get_count();
			for ( int c = 0; c < ( n - 1 ); c++ )
			{
				for ( int ofs = 0; ofs < n - c - 1; ofs++ )
				{
					
					if ( get_start( ofs ) > get_start( ofs + 1 ) )
						swap( ofs, ofs + 1 );
				}
			}
			sorted = true;
		}
	}
}


class merge_interval_list extends sortable_interval_list
{
	public void merge_interval_list( final int n )
	{
		sortable_interval_list( n );
	}

	// removes elements with len==0 ( created by merge_intervals() )
	private void compact()
	{
		int n = get_count();
		if ( n > 1 )
		{
			int dst = 0;
			int src = 0;
			
			while ( src < n )
			{
				if ( dst == src )
				{
					if ( get_len( dst ) > 0 ) dst++;
				}
				else
				{
					if ( get_len( src ) > 0 )
					{
						set( dst, get_start( src ), get_len( src ), is_amb( src ) );
						set_len( src, 0 );
						dst++;
					}
				}
				src++;
			}
			// reduce the count if we have removed elements
			if ( src > dst ) set_count( n - ( src - dst ) );
		}
	}

	public void merge_intervals()
	{
		int n = get_count();
		if ( n > 1 )
		{
			int last = 0;
			int curr = 1;

			sort();
			while ( curr < n )
			{
				long end_last = get_start( last ) + get_len( last );
				if ( end_last >= get_start( curr ) )
				{
					long end_curr = get_start( curr ) + get_len( curr );
					if ( end_curr > end_last )
                    {
                        long value = get_start( curr ) - get_start( last );
                        value += get_len( curr );
						set_len( last, (int)value );
                    }
					set_len( curr, 0 );
				}
				else
				{
					last++;
					while( get_len( last ) == 0 ) last++;
				}
				curr++;
			}
			compact();
		}
	}
}


class event_iter_list extends merge_interval_list
{
	private int interval_idx;
	private boolean request_start;
	long tmp[];
	
	public void event_iter_list( final int n )
	{
		merge_interval_list( n );
		tmp = new long[ 2 ];
		reset_event_iter();
	}

	public void reset_event_iter()
	{
		interval_idx = 0;
		request_start = true;
	}

	boolean start_event( int idx, long event[] )
	{
		boolean res = ( idx < get_count() );
		if ( res )
		{
			event[ 0 ] = get_start( idx );
			event[ 1 ] = is_amb( idx ) ? 2 : 1;
		}
		return res;
	}
	
	boolean end_event( int idx, long event[] )
	{
		boolean res = ( idx < get_count() );
		if ( res )
		{
			event[ 0 ] = get_end( idx );
			event[ 1 ] = is_amb( idx ) ? 4 : 3;
		}
		return res;
	}
	

	// event[ 0 ] ... position
	// event[ 1 ] ... 1=start, 2=start_amb, 3=end, 4=end_amb
	public boolean get_event( long event[] )
	{
		boolean res;
		if ( request_start )
		{
			res = start_event( interval_idx, event );
			if ( res ) request_start = false;
		}
		else
		{
			res = end_event( interval_idx, event );
			if ( res )
			{
				boolean has_nxt = start_event( interval_idx + 1, tmp );
				if ( has_nxt )
				{
					boolean follows = ( tmp[ 0 ] == ( event[ 0 ] + 1 ) );
					if ( follows )
					{
						event[ 0 ] = tmp[ 0 ];
						event[ 1 ] = tmp[ 1 ];
					}
					else
						request_start = true;
				}
				else
					request_start = true;
				interval_idx++;
			}
		}
		return res;
	}
}


public class interval_list extends event_iter_list
{
	/* constructor to create an empty interval-list with n entries */
	interval_list( final int n )
	{
		event_iter_list( n );
	}

	/* copy constructor to create a identical copy of another interval list */
	interval_list( final interval_list other )
	{
		event_iter_list( other.get_size() );
		int n = other.get_count();
		for( int i = 0; i < n; i++ )
			add( other.get_start( i ), other.get_len( i ), other.is_amb( i ) );
	}

	/* constructor to create an interval-list from a set of given intervals, START1, LEN1, START2, LEN2 */
	interval_list( final int n, final long...setOfInts )
	{
		event_iter_list( n );
		
		long i_start = 0;
		boolean has_start = false;
		for( long value : setOfInts )
		{
			if ( has_start )
				add( i_start, (int)value, false );
			else
				i_start = value;
			has_start = !has_start;
		}
	}

	/* constructor from a cigar-string */
	interval_list( final int n, final long start, final String cigar )
	{
		event_iter_list( n );
		from_cigar( start, cigar );
	}

	
	/* helper to fill a list with values from a cigar-string */
	void from_cigar( final long start, final String cigar )
	{
		String s_num = "";
		long ofs = 0;
		for ( char c : cigar.toCharArray() )
		{
			if ( Character.isDigit( c ) )
				s_num += c;
			else
			{
				int len = 1;
				if ( !s_num.isEmpty() ) len = Integer.parseInt( s_num );
				s_num = "";
				
				switch( c )
				{
					case 'D' : ;
					case 'N' : break;			
					case 'M' : ;
					case '=' : ;
					case 'X' : add( start + ofs, len, false ); break;
					default  : len = 0; break; // includes H,P,I,S
				}
				ofs += len;
			}
		}
		merge_intervals();
	}
	
	void add( final interval_list other )
	{
		int n = other.get_count();
		for ( int i = 0; i < n; i++ )
			add( other.get_start( i ), other.get_len( i ), other.is_amb( i ) );
	}
	

	void detect_amb( final interval_list other )
	{
		int other_idx = 0;
		int other_count = other.get_count();
		int this_idx = 0;
        long other_start = other.get_start( 0 );
		long other_end = other.get_end( 0 );
		long this_start = get_start( 0 );
		long this_end = get_end( 0 );
		
		while ( other_idx < other_count && this_idx < get_count() )
		{
			if ( other_start <= this_start )
			{
				// could be AMB_01, AMB_03 or AMB_04

				if ( other_end < this_start )
				{
					//  CASE AMB_01: no overlap between other and this
					//	other:	XXXXXXXXXXXXXXX
					//	this:                     XXXXXXXXXXXXXXXX
					other_start = other.get_start( ++other_idx );		// ---> next from other
					other_end = other.get_end( other_idx );
				}
				else if ( other_end < this_end )
				{
					// CASE AMB_04:	partial overlap between other and this
					//	other:	XXXXXXXXXXXXXXXX
					//	this:         XXXXXXXXXXXXXXXX
					//	-----------------------------------
					//	this:         AAAAAAAAAAXXXXXX
					if ( !is_amb( this_idx ) )
					{	// split the current interval in 2
                        long value = other_end - this_start + 1;
						insert( this_idx++, this_start, (int)value, true );
						this_start = other_end + 1;
                        value = this_end - other_end;
						set( this_idx, this_start, (int)value, false );
					}
					other_start = other.get_start( ++other_idx );		// ---> next from other
					other_end = other.get_end( other_idx );
				}
				else
				{
					//  CASE AMB_03: complete overlap between other and this
					//	other:	XXXXXXXXXXXXXXXXXXXXXXXXXXX
					//	this:         XXXXXXXXXXXXXXXX
					//	-----------------------------------
					//	this:         AAAAAAAAAAAAAAAA
					if ( !is_amb( this_idx ) )
						set_amb( this_idx, true );		// set the amb-flag for the whole interval
					if ( this_end == other_end )
					{
						other_start = other.get_start( ++other_idx );	// ( maybe also next from other )
						other_end = other.get_end( other_idx );
					}
					this_start = get_start( ++this_idx );		// ---> next from this
					this_end = get_end( this_idx );
				}
			}
			else
			{
				// could be AMB_02, AMB_05 or AMB_06
				if ( other_start > this_end )
				{
					//  CASE AMB_02: no overlap between other and this
					//	other:	                  XXXXXXXXXXXXXXX
					//	this:   XXXXXXXXXXXXXXXX
					this_start = get_start( ++this_idx );		// ---> next from this
					this_end = get_end( this_idx );
				}
				else if ( other_end < this_end )
				{
					//	CASE AMB_06: partial overlap between other and this
					//	other:	           XXXXX
					//	this:         XXXXXXXXXXXXXXXX
					//	-----------------------------------
					//	this:         XXXXXAAAAAXXXXXX
					if ( !is_amb( this_idx ) )
					{	// split the current interval in 3
                        long value = other_start - this_start;
						set_len( this_idx++, (int)value );
                        value = other_end - other_start + 1;
						insert( this_idx++, other_start, (int)value, true );
						this_start = other_end + 1;
                        value = this_end - other_end;
						insert( this_idx, this_start, (int)value, false );
					}
					other_start = other.get_start( ++other_idx );	// ---> next from other
					other_end = other.get_end( other_idx );
				}
				else
				{
					//  CASE AMB_05: partial overlap between other and this
					//	other:	             XXXXXXXXXXXXXXXX
					//	this:         XXXXXXXXXXXXXXXX
					//	-----------------------------------
					//	this:         XXXXXXXAAAAAAAAA
					if ( !is_amb( this_idx ) )
					{	// split the current interval in 2
                        long value = other_start - this_start;
						set_len( this_idx++, (int)value );
                        value = this_end - other_start + 1;
						insert( this_idx, other_start, (int)value, true );
					}
					this_start = get_start( ++this_idx );		// ---> next from this
					this_end = get_end( this_idx );
				}
			}
		}
	}

}
