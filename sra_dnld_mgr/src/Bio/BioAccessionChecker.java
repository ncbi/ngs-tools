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
import javax.swing.SwingWorker;
import ngs.ErrorMsg;
import ngs.ReadCollection;
import ngs.ReferenceSequence;

/* a SwingWorker that produces an BioSpec, and publishes Ints */
public class BioAccessionChecker extends SwingWorker< BioSpec, Integer >
{
    private final String accession;
    
    private boolean supported_and_valid()
    {
        boolean res = gov.nih.nlm.ncbi.ngs.NGS.isSupported();
        if ( !res )
            CLogger.log( "verifier: NGS not supported" );
        return res;
    }

    @Override protected BioSpec doInBackground()
    {
        BioSpec res = new BioSpec();
        res.set_accession( accession );
        if ( supported_and_valid() && !isCancelled() )
        {
            try
            {
                ReadCollection run = gov.nih.nlm.ncbi.ngs.NGS.openReadCollection( accession );
                if ( !isCancelled() )
                {
                    res.set_count( run.getReadCount() );
                    if ( !isCancelled() )
                    {
                        if ( run.getAlignmentCount() > 0 )
                            res.set_type( BioAccessionType.READ_COLLECTION_ALIGNED );
                        else
                            res.set_type( BioAccessionType.READ_COLLECTION_UNALIGNED );
                    }
                }
            }
            catch ( ErrorMsg ex )
            {  }
            
            if ( !isCancelled() && res.get_type().equals( BioAccessionType.INVALID ) )
            {
                try
                {
                    ReferenceSequence ref = gov.nih.nlm.ncbi.ngs.NGS.openReferenceSequence( accession );
                    if ( !isCancelled() )
                    {
                        res.set_count( ref.getLength() );
                        if ( !isCancelled() )
                            res.set_type( BioAccessionType.REF_SEQUENCE );
                    }
                }
                catch ( ErrorMsg ex )
                {  }
            }
        }
        return res;
    }

    public BioAccessionChecker( final String accession )
    {
        this.accession = accession;
    }
}
