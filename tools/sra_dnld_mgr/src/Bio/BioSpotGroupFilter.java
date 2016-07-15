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

public class BioSpotGroupFilter implements BioFilter
{
    private final String spotgroup;
    
    @Override public boolean pass( Read obj, final String bases )
    {
        boolean res;
        try
        {
            res = spotgroup.equals( obj.getReadGroup() );
        }
        catch ( ErrorMsg ex )
        {
            res  = false;
        }
        return res;
    }

    @Override public boolean pass( Read read_obj, Fragment frag_obj, final String bases )
    {
        boolean res;
        try
        {
            res = spotgroup.equals( read_obj.getReadGroup() );
        }
        catch ( ErrorMsg ex )
        {
            res  = false;
        }
        return res;
   }

    @Override public boolean pass( Alignment obj, final String bases )
    {
        boolean res;
        try
        {
            res = spotgroup.equals( obj.getReadGroup() );
        }
        catch ( ErrorMsg ex )
        {
            res  = false;
        }
        return res;
    }

    @Override public boolean pass( ReferenceSequence obj, final String bases, long position )
    {
        return true;
    }
    
    public BioSpotGroupFilter( final String spotgroup )
    {
        this.spotgroup = spotgroup;
    }

}
