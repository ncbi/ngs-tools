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

public class ref_cmp_res implements java.io.Serializable
{
	public Iterable< StringPair > in_gtf_and_acc;
	public Iterable< StringPair > only_in_acc;
	public Iterable< String > only_in_gtf;

	int longest_string()
	{
		int res = 0;
		
		for ( StringPair sp : in_gtf_and_acc )
		{
			int l = sp.toString().length();
			if ( l > res ) res = l;
		}

		for ( StringPair sp : only_in_acc )
		{
			int l = sp.toString().length();
			if ( l > res ) res = l;
		}

		for ( String s : only_in_gtf )
		{
			int l = s.length();
			if ( l > res ) res = l;
		}
		
		return res;
	}
	
	public String translate( final String s )
	{
		String res = null;
		for ( StringPair sp : in_gtf_and_acc )
		{
			if ( res == null )
			{
				String S1 = sp.get_S1();
				String S2 = sp.get_S2();
				if ( S1 != null && S2 != null )
				{
					if ( s.equals( S1 ) )
						res = S2;
					else if ( s.equals( S2 ) )
						res = S1;
				}
			}
		}
		return res;
	}
	
	public boolean is_in_gtf_and_acc( final String s )
	{
		boolean res = false;
		for ( StringPair sp : in_gtf_and_acc )
		{
			if ( !res )
			{
				String S1 = sp.get_S1();
				String S2 = sp.get_S2();
				boolean c1 = ( S1 == null ) ? false : s.equals( S1 );
				boolean c2 = ( S2 == null ) ? false : s.equals( S2 );
				res = ( c1 || c2 );
			}
		}
		return res;
	}
}