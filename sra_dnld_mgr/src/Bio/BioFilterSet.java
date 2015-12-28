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

import java.util.ArrayList;
import ngs.*;

public class BioFilterSet
{
    private final ArrayList< BioFilter > filters;
            
    public boolean pass ( BioStats stats, Read obj, String bases )
    {
        boolean res = true;
        for ( int i = 0; res && i < filters.size(); i++ )
        {
            res = filters.get( i ).pass( obj , bases );
        }
        if ( !res && stats != null ) stats.number_rejected++;
        return res;
    }

    public boolean pass ( BioStats stats, Read read_obj, Fragment frag_obj, String bases )
    {
        boolean res = true;
        for ( int i = 0; res && i < filters.size(); i++ )
        {
            res = filters.get( i ).pass( read_obj , frag_obj , bases );
        }
        if ( !res && stats != null ) stats.number_rejected++;        
        return res;
    }

    public boolean pass ( BioStats stats, Alignment obj, String bases )
    {
        boolean res = true;
        for ( int i = 0; res && i < filters.size(); i++ )
        {
            res = filters.get( i ).pass( obj , bases );
        }
        if ( !res && stats != null ) stats.number_rejected++;        
        return res;
    }

    public boolean pass ( BioStats stats, ReferenceSequence obj, String bases, long position )
    {
        boolean res = true;
        for ( int i = 0; res && i < filters.size(); i++ )
        {
            res = filters.get( i ).pass( obj , bases, position );
        }
        if ( !res && stats != null ) stats.number_rejected++;        
        return res;
    }

    public int count()
    {
        return filters.size();
    }
    
    public void addFilter ( BioFilter filter )
    {
        if ( filter != null )
            filters.add( filter );
    }
    
    public BioFilterSet()
    {
        filters = new ArrayList<>();
    }
}
