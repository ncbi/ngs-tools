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

import java.util.LinkedList;
import java.util.Queue;

public class BioRead
{
    private final Queue< BioRecord > q;
    public long read_id;

    public boolean is_empty() { return q.isEmpty(); }
    public void put_record( final BioRecord rec ) { q.offer( rec ); }
    public BioRecord get_record() { return q.poll(); }
    public void clear() { q.clear(); read_id = 0; }
            
    public void parse_read_id( final String id )
    {
        // SRR000001.R.1
        String[] parts = id.split( "\\." );
        if ( parts.length > 2 ) read_id = Long.parseLong( parts[ 2 ] );
    }

    public BioRead()
    {
        q = new LinkedList<>(); 
        read_id = 0;
    }
}
