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

public interface BioFilter
{
    /**
     * pass
     * @param obj read iterator to inspect
     * @param bases bases to inspect
     * @return boolean if obj has passed this filter
     */
    boolean pass ( Read obj, String bases );

    /**
     * pass
     * @param frag_obj fragment iterator to inspect
     * @param bases bases to inspect
     * @return boolean if obj has passed this filter
     */
    //boolean pass ( Fragment obj, String bases );
    boolean pass ( Read read_obj, Fragment frag_obj, String bases );
    
    /**
     * pass
     * @param obj alignment iterator to inspect
     * @param bases bases to inspect
     * @return boolean if obj has passed this filter
     */
    boolean pass ( Alignment obj, String bases );
    
    /**
     * pass
     * @param obj RefSeq-obj to inspect
     * @param bases bases to inspect
     * @param position to inspect
     * @return boolean if obj has passed this filter
     */
    boolean pass ( ReferenceSequence obj, String bases, long position );
}
