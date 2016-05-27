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

package job;

import Bio.BioRead;
import Bio.BioRecord;
import data.QConnect;
import java.util.concurrent.ConcurrentLinkedQueue;

public class BioRunState extends QConnect< BioRead >
{
    private final ConcurrentLinkedQueue< BioRecord > empty_records;
    private boolean running;
    
    public synchronized boolean is_running() { return running; }
    public synchronized void set_running( final boolean value ) { running = value; }
    public synchronized void put_record_back( final BioRecord rec ) { empty_records.offer( rec ); }
    
    @Override public BioRead make_entry() { return new BioRead(); }
    
    public synchronized BioRecord get_empty_record()
    {
        BioRecord res;
        if ( empty_records.isEmpty() )
        {
            res = new BioRecord();
        }
        else
        {
            res = empty_records.poll();
            res.clear();
        }
        return res;
    }
   
    @Override public synchronized void clear()
    {
        super.clear();
        empty_records.clear();
    }
    
    public BioRunState( int max_created )
    {
        super( max_created );
        empty_records = new ConcurrentLinkedQueue<BioRecord>();
        running = false;
    }
    
}
