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

import java.util.List;
import javax.swing.JTextArea;
import javax.swing.SwingWorker;
import job.JobData;

public class BioDumperPreview  extends SwingWorker< Integer, String >
{
    private final JobData job;
    private final BioRead read;
    private BioRecord rec;
    private final JTextArea txt;
    private boolean running;
    private int progress;
    private final int rows;
    
    private long start_value()
    {
        long res = 0;
        BioAccessionType t = job.get_bio_type();
        switch( t )
        {
            case READ_COLLECTION_UNALIGNED :
            case READ_COLLECTION_ALIGNED   : res = 1;
        }
        return res;
    }
    
    private int produce( final BioDumper dumper )
    {
        int res = 0;
        running = dumper.next();
        if ( running )
        {
            read.clear();
            running = dumper.produce( read );
            while ( !read.is_empty() )
            {
                rec = read.get_record();
                if ( rec != null )
                {
                    if ( !rec.is_empty() )
                    {
                        publish( rec.get() );
                        res++;
                    }
                }
            }
        }
        return res;
    }
    
    @Override protected Integer doInBackground()
    {
        final BioDumper dumper = BioDumperFactory.make_dumper( job, start_value() );
        progress = 0;
        if ( dumper != null )
        {
            running = true;
            while( running )
            {
                if ( produce( dumper ) > 0 )
                    setProgress( ++progress );
                running = ( progress < rows );
            }
        }
        return progress;
    }
    
    /* this runs in the event-loop and has access to Swing-components */
    @Override protected void process( final List< String > chunks )
    {
        for ( final String s : chunks ) txt.append( s );   
    }
    
    public BioDumperPreview( final JobData job,
            final JTextArea txt, final int rows )
    {
        this.job = job;
        this.rows = rows;
        read = new BioRead();
        rec = null;
        this.txt = txt;
        progress = 0;
    }

}
