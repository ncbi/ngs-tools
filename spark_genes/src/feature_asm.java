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

/*
    - the feature_asm wraps a gtf_reader to assmeble multiple gtf-lines into a feature
    - it can be reset by calling reset() ... calls reset on it's gtf_reader

    - it is an iterator, the user calls next() to get the next feature
    - next() returns null when the file is processed, if successful it returns a new gtf_feature
    - the iterator groups multiple gtf-lines based on having the same ID into one feature
    - the iterator sorts the start/len pairs by start-position before creating a new gtf_feature
    - the iterator merges overlapping start/len pairs before creating a new gtf_feature
*/

public class feature_asm
{
    private gtf_reader reader;
    private String[] assembled;
    private int[] start;
    private int[] len;
    private int ranges, features;
    private String[] current;
    private boolean new_feature;
    private boolean reader_ok;

    feature_asm( gtf_reader reader )
    {
        this.reader = reader;
        assembled = new String[ 4 ];
        start = new int[ 2048 ];
        len   = new int[ 2048 ];
        current = new String[ 5 ];
        ranges = 0;
        features = 0;
        new_feature = true;
        if ( reader != null )
            reader_ok = reader.next( current );
        else
            reader_ok = false;
    }

    void reset()
    {
        ranges = 0;
        features = 0;
        new_feature = true;
        if ( reader != null )
        {
            reader.reset();
            reader_ok = reader.next( current );
        }
        else
            reader_ok = false;
    }

    int get_features() { return features; }

	private void swap( final int idx1, final int idx2 )
	{
		int temp = start[ idx1 ];
        start[ idx1 ] = start[ idx2 ];
        start[ idx2 ] = temp;

		temp = len[ idx1 ];
        len[ idx1 ] = len[ idx2 ];
        len[ idx2 ] = temp;
	}

	private void sort()
	{
        for ( int c = 0; c < ( ranges - 1 ); c++ )
        {
            for ( int ofs = 0; ofs < ranges - c - 1; ofs++ )
            {
                if ( start[ ofs ] > start[ ofs + 1 ] )
                    swap( ofs, ofs + 1 );
            }
        }
	}

	private void merge()
    {
        int last = 0;
        int curr = 1;
        while ( curr < ranges )
        {
            int end_last = start[ last ] + len[ last ];
            if ( end_last >= start[ curr ] )
            {
                int end_curr = start[ curr ] + len[ curr ];
                if ( end_curr > end_last )
                    len[ last ] = len[ curr ] + ( start[ curr ] - start[ last ] );
                len[ curr ] = 0;
            }
            else
            {
                last++;
                while( len[ last ] == 0 ) last++;
            }
            curr++;
        }
    }

    private void store_ranges()
    {
        start[ ranges ] = Integer.parseInt( current[ 1 ] );
        len[ ranges ] = ( Integer.parseInt( current[ 2 ] ) - start[ ranges ] ) + 1;
        ranges++;
    }

    // 0...done, 1...skip, 2...valid feature assembled
    private int read_feature()
    {
        if ( reader_ok )
        {
            if ( new_feature )
            {
                assembled[ 0 ] = current[ 4 ];  // id
                assembled[ 1 ] = current[ 3 ];  // strand
                assembled[ 2 ] = current[ 0 ];  // chromosome
                ranges = 0;
                store_ranges();
                new_feature = false;
                return 1;
            }
            else
            {
                reader_ok = reader.next( current );
                if ( reader_ok )
                {
                    if ( !assembled[ 0 ].equals( current[ 4 ] ) )
                        return 2;
                    else
                        store_ranges();
                    return 1;
                }
            }
        }
        return 0;
    }

    gtf_feature next()
    {
        int rf;
        while ( ( rf = read_feature() ) != 0 )
        {
            if ( rf == 2 )
            {
                new_feature = true;
                features++;

                sort();
                merge();

                gtf_feature res = new gtf_feature( ranges, assembled[ 0 ], assembled[ 2 ], assembled[ 1 ].equals( "-" ) );
                for ( int src = 0; src < ranges; ++src )
                    res.add_seg( start[ src ], len[ src ] );
                return res;
            }
        }
        return null;
    }
    
}