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

public enum run_mode
{
	LOCAL_P     { public String toString() { return "LOCAL.PARALLEL"; } },
	LOCAL_S     { public String toString() { return "LOCAL.SERIAL"; } },
	SPARKED     { public String toString() { return "SPARKED"; } },
	PRE_GTF     { public String toString() { return "PRE_GTF"; } },
	REFS        { public String toString() { return "REFS"; } };
	
	public static run_mode from_string( String s )
	{
		if ( s.equals( "LOCAL_P" ) )
			return LOCAL_P;
		else if ( s.equals( "local_p" ) )
			return LOCAL_P;
		else if ( s.equals( "LOCAL_S" ) )
			return LOCAL_S;
		else if ( s.equals( "local_s" ) )
			return LOCAL_S;
		else if ( s.equals( "SPARKED" ) )
			return SPARKED;
		else if ( s.equals( "sparked" ) )
			return SPARKED;
		else if ( s.equals( "PRE_GTF" ) )
			return PRE_GTF;
		else if ( s.equals( "pre_gtf" ) )
			return PRE_GTF;
		else if ( s.equals( "REFS" ) )
			return REFS;
		else if ( s.equals( "refs" ) )
			return REFS;

		return LOCAL_P;
	}
	
}
