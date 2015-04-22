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
import java.util.concurrent.*;

class feature_bin implements java.io.Serializable
{
    private ArrayList< gtf_feature > features;

    feature_bin( gtf_feature f )
    {
        features = new ArrayList< gtf_feature >();
        add( f );
    }

    public void add( gtf_feature f ) { features.add( f ); }

    public int lookup( long start, long end, Map< String, gtf_feature > collector )
    {
        int res = 0;
        for ( gtf_feature f : features )
        {
            if ( start <= f.get_end() && end >= f.get_start() )
            {
                if ( !collector.containsKey( f.get_id() ) )
                {
                    collector.put( f.get_id(), f );
                    res++;
                }
            }
        }
        return res;
    }

    public int get_count() { return features.size(); }
}


class ref_features implements java.io.Serializable
{
    private final String ref_name;
    private Map< Long, feature_bin > bins;
    private Map< String, gtf_feature > collector;
    private long features_start, features_end, features_max_len, bin_size;
    private int count;

    ref_features( gtf_feature f, long bin_size )
    {
        this.ref_name = f.get_ref();
        this.bin_size = bin_size;
        bins = new ConcurrentHashMap< Long, feature_bin >();
        collector = new ConcurrentHashMap< String, gtf_feature >();

        features_start = Long.MAX_VALUE;
        features_end = 0;
        features_max_len = 0;
        count = 0;

        add( f );
    }

    private void add_to_bin( long bin_id, gtf_feature f )
    {
        feature_bin bin = bins.get( bin_id );
        if ( bin != null )
            bin.add( f );
        else
            bins.put( bin_id, new feature_bin( f ) );
    }

    public void add( gtf_feature f )
    {
        long start = f.get_start();
        long end = f.get_end();
        long len = ( end - start ) + 1;

        if ( start < features_start ) features_start = start;
        if ( end > features_end ) features_end = end;
        if ( len > features_max_len ) features_max_len = len;

        long start_bin = start / bin_size;
        long end_bin = end / bin_size;
        long bin_id;
        for ( bin_id = start_bin; bin_id <= end_bin; ++ bin_id )
            add_to_bin( bin_id, f );
        count++;
    }

    public int lookup( long start, long len, List< gtf_feature > found )
    {
        int res = 0;
        long end = start + len - 1;
        if ( start <= features_end && end >= features_start )
        {
            // the requested rangs is not outside of what we have stored
            long start_bin = start / bin_size;
            long end_bin = end / bin_size;
            long bin_id;
            for ( bin_id = start_bin; bin_id <= end_bin; ++ bin_id )
            {
                collector.clear();
                feature_bin bin = bins.get( bin_id );
                if ( bin != null )
                    res += bin.lookup( start, end, collector );
            }
            for ( gtf_feature f : collector.values() )
                found.add( f );
        }
        return res;
    }

    public void report()
    {
        int min_per_bin = Integer.MAX_VALUE;
        int max_per_bin = 0;
        for ( feature_bin bin : bins.values() )
        {
            int n = bin.get_count();
            if ( n < min_per_bin ) min_per_bin = n;
            if ( n > max_per_bin ) max_per_bin = n;
        }

        System.out.println( ref_name + ": " + bins.size() + 
                            " bins [ " + features_start + " ... " + features_end + " ] maxlen = " +
                            features_max_len + " features = " + count + 
                            " ( min= " + min_per_bin + " max= " + max_per_bin + " )" );
    }
}

// the lookup-table that contains all the features...
public class gtf_lookup implements java.io.Serializable
{
    private Map< String, ref_features > lookup;
    private long bin_size;
	private ref_cmp_res refs;

    // the constructor takes a flat list of features and constructs the lookup-tree
    gtf_lookup( Iterable<gtf_feature> features, long bin_size, ref_cmp_res refs )
    {
        this.bin_size = bin_size;
		this.refs = refs;
        lookup = new ConcurrentHashMap< String, ref_features >();
        for ( gtf_feature f : features ) add( f );
    }
    
    // adding a feature splits into 2 cases, the reference has been seen already or not
    void add( gtf_feature f )
    {
		String gtf_ref = f.get_ref();
		if ( refs.is_in_gtf_and_acc( gtf_ref ) )
		{
			ref_features rf = lookup.get( gtf_ref );
			if ( rf != null )
			{
				// reference-features are already in lookup-map
				rf.add( f );
			}
			else
			{
				// reference-features have to be created...
				rf = new ref_features( f, bin_size );
				lookup.put( gtf_ref, rf );

				String acc_ref = refs.translate( gtf_ref );
				if ( acc_ref != null )
					lookup.put( acc_ref, rf );
			}
		}
    }

    // lookup puts the features intersecting the interval into the given list...
    int lookup( String ref, long start, long len, List< gtf_feature > found )
    {
        int res = 0;
        found.clear();
        ref_features rf = lookup.get( ref );
        if ( rf != null )
            res = rf.lookup( start, len, found );
        return res;
    }

    void report()
    {
        for ( ref_features rf : lookup.values() ) rf.report();
    }
}
