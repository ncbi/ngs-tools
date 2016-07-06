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

package Bio;

import ngs.*;

public class BioMaxNFilter implements BioFilter
{
    private final int max_N;
    
    private int count_N( final String s )
    {
        int res = 0;
        for( char c : s.toCharArray() )
        {
            if ( c == 'N' ) res++;
        }
        return res;
    }

    @Override  public boolean pass( Read obj, final String bases )
    {
        return ( count_N( bases ) < max_N );
    }

    @Override public boolean pass( Read read_obj, Fragment frag_obj, final String bases )
    {
        return ( count_N( bases ) < max_N );
    }

    @Override public boolean pass( Alignment obj, final String bases )
    {
        return ( count_N( bases ) < max_N );
    }

    @Override public boolean pass( ReferenceSequence obj, final String bases, long position )
    {
        return ( count_N( bases ) < max_N );
    }
    
    public BioMaxNFilter( final int max_N )
    {
        this.max_N = max_N;
    }
}
