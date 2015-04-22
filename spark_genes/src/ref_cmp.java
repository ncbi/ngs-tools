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

public class ref_cmp
{
	private static Set< String > make_gtf_refs( final Iterable< gtf_feature > features )
	{
        Set< String >res = new TreeSet< String >();
        for ( gtf_feature f : features )
        {
			String ref = f.get_ref();
			if ( !res.contains( ref ) ) res.add( ref );
        }
		return res;
	}
	
	private static Iterable< StringPair > make_acc_refs( final String accession )
	{
        LinkedList< StringPair >res = new LinkedList< StringPair >();
        try
        {
            ReadCollection ngs_run = gov.nih.nlm.ncbi.ngs.NGS.openReadCollection( accession );
            ReferenceIterator ngs_refs = ngs_run.getReferences();
            while( ngs_refs.nextReference() )
				res.add( new StringPair( ngs_refs.getCanonicalName(), ngs_refs.getCommonName() ) );	
        }
        catch ( ErrorMsg err )
        {
            System.err.println( err.toString() );
        }
		return res;
	}
	
	private static boolean gtf_refs_contain( final StringPair sp, final Set< String >gtf )
	{
		boolean c1 = sp.get_S1() == null ? false : gtf.contains( sp.get_S1() );
		boolean c2 = sp.get_S2() == null ? false : gtf.contains( sp.get_S2() );
		return ( c1 || c2 );
	}
	
	private static boolean acc_refs_contain( final String s, final Iterable< StringPair >acc  )
	{
		boolean res = false;
		for ( StringPair sp_acc : acc )
		{
			boolean c1 = sp_acc.get_S1() == null ? false : s.equals( sp_acc.get_S1() );
			boolean c2 = sp_acc.get_S2() == null ? false : s.equals( sp_acc.get_S2() );
			if ( c1 || c2 ) res = true;
		}
		return res;
	}

	private static Iterable< StringPair > in_both( final Set< String >gtf, final Iterable< StringPair >acc )
	{
		LinkedList< StringPair > res = new LinkedList< StringPair >();
		for ( StringPair sp : acc )
		{
			if ( gtf_refs_contain( sp, gtf ) )
				res.add( sp );
		}
		return res;
	}

	private static Iterable< StringPair > only_in_acc( final Set< String >gtf, final Iterable< StringPair >acc )
	{
		LinkedList< StringPair > res = new LinkedList< StringPair >();
		for ( StringPair sp : acc )
		{
			if ( !gtf_refs_contain( sp, gtf ) )
				res.add( sp );
		}
		return res;
	}
	
	private static Iterable< String > only_in_gtf( final Set< String >gtf, final Iterable< StringPair >acc )
	{
		Set< String > res = new TreeSet< String >();
		for ( String s : gtf )
		{
			if ( !acc_refs_contain( s, acc ) )
				res.add( s );
		}
		return res;
	}
	
	static ref_cmp_res compare( final Iterable< gtf_feature > features, final String accession )
	{
		Set< String > gtf = make_gtf_refs( features );
		Iterable< StringPair > acc = make_acc_refs( accession );

		ref_cmp_res res = new ref_cmp_res();
		res.in_gtf_and_acc = in_both( gtf, acc );
		res.only_in_acc = only_in_acc( gtf, acc );
		res.only_in_gtf = only_in_gtf( gtf, acc );

		return res;
	}
	
}
