/* ===========================================================================
#
#                            PUBLIC DOMAIN NOTICE
#               National Center for Biotechnology Information
#
#  This software/database is a "United States Government Work" under the
#  terms of the United States Copyright Act.  It was written as part of
#  the author's official duties as a United States Government employee and
#  thus cannot be copyrighted.  This software/database is freely available
#  to the public for use. The National Library of Medicine and the U.S.
#  Government have not placed any restriction on its use or reproduction.
#
#  Although all reasonable efforts have been taken to ensure the accuracy
#  and reliability of the software and data, the NLM and the U.S.
#  Government do not and cannot warrant the performance or results that
#  may be obtained by using this software or data. The NLM and the U.S.
#  Government disclaim all warranties, express or implied, including
#  warranties of performance, merchantability or fitness for any particular
#  purpose.
#
#  Please cite the author in any work or product based on this material.
#
=========================================================================== */
package Bio;

import data.CLogger;
import java.util.ArrayList;
import java.util.List;
import ngs.*;

public class BioVerifier
{
    private boolean supported_and_valid( final String source )
    {
        boolean res = gov.nih.nlm.ncbi.ngs.NGS.isSupported();
        if ( !res )
            CLogger.log( "verifier: NGS not supported" );
        if ( res )
        {
            res = gov.nih.nlm.ncbi.ngs.NGS.isValid( source );
            if ( !res )
                CLogger.logfmt( ">%s< is an invalid NGS accession" , source );
        }
        return res;
    }
    
    public BioSpec verify( final String source )
    {
        BioAccessionType biotype = BioAccessionType.INVALID;
        long count = 0;
        
        boolean res = supported_and_valid( source );
        if ( res )
        {
            try
            {
                ReadCollection run = gov.nih.nlm.ncbi.ngs.NGS.openReadCollection( source );
                count = run.getReadCount();
                if ( run.getAlignmentCount() > 0 )
                    biotype = BioAccessionType.READ_COLLECTION_ALIGNED;
                else
                    biotype = BioAccessionType.READ_COLLECTION_UNALIGNED;
            }
            catch ( ErrorMsg ex )
            {  }
            
            if ( biotype.equals( BioAccessionType.INVALID ) )
            {
                try
                {
                    ReferenceSequence ref = gov.nih.nlm.ncbi.ngs.NGS.openReferenceSequence( source );
                    count = ref.getLength();
                    biotype = BioAccessionType.REF_SEQUENCE;
                }
                catch ( ErrorMsg ex )
                {  }
            }
        }
        return new BioSpec( biotype, count, source );
    }
    
    public List<String> extract_spotgroups( final String source )
    {
        List<String> res = new ArrayList<String>();
        if ( supported_and_valid( source ) )
        {
            try
            {
                ReadCollection run = gov.nih.nlm.ncbi.ngs.NGS.openReadCollection( source );
                ReadGroupIterator iter = run.getReadGroups();
                while( iter.nextReadGroup() )
                {
                    if ( iter.getName().length() > 0 )
                        res.add( iter.getName() );
                }
            }
            catch ( ErrorMsg ex )
            {
                CLogger.logfmt( "errro %s opening %s", ex.toString(), source );
            }
        }
        return res;
    }

    public BioVerifier()
    {
    }
}
