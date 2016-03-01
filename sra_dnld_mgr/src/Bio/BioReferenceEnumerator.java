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

import GUI.BioRefSpecReceiver;
import data.CLogger;
import java.util.List;
import javax.swing.SwingWorker;
import ngs.ErrorMsg;
import ngs.ReadCollection;
import ngs.ReferenceIterator;

public class BioReferenceEnumerator extends SwingWorker< Integer, BioRefSpec >
{
    private final String accession;
    private final BioRefSpecReceiver receiver;
    private int progress;

    public String get_accession() { return accession; }
    
    private boolean supported_and_valid()
    {
        boolean res = gov.nih.nlm.ncbi.ngs.NGS.isSupported();
        if ( !res )
            CLogger.log( "BioReadGroupEnumerator: NGS not supported" );
        return res;
    }

    private void publish_item( final BioRefSpec spec )
    {
        //CLogger.logfmt( "publish reference: %s", spec.name() );
        publish( spec );
        setProgress( ++progress );
    }

    /* this runs in the background in a hidden thread */
    @Override protected Integer doInBackground()
    {
        if ( supported_and_valid() )
        {
            try
            {
                ReadCollection run = gov.nih.nlm.ncbi.ngs.NGS.openReadCollection( accession );
                ReferenceIterator iter = run.getReferences();
                while ( iter.nextReference() )
                {
                    publish_item( new BioRefSpec(
                        iter.getCanonicalName(),
                        iter.getCommonName(),
                        iter.getIsCircular(),
                        iter.getLength() ) );
                    progress++;
                }
            }
            catch ( ErrorMsg ex )
            {
                CLogger.logfmt( "errro %s opening %s", ex.toString(), accession );
            }
        }
        return progress;
    }

    /* this runs in the event-loop and has access to Swing-components */
    @Override protected void process( final List< BioRefSpec > chunks )
    {
        for ( final BioRefSpec spec : chunks ) receiver.on_received_spec( spec );   
    }

    public BioReferenceEnumerator( final BioRefSpecReceiver receiver, final String accession )
    {
        this.accession = accession;
        this.receiver = receiver;
        progress = 0;
    }
}
