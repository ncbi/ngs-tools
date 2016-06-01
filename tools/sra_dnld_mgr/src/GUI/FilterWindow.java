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

package GUI;

import job.JobData;
import java.awt.*;
import javax.swing.*;

public class FilterWindow extends DlgWithMaxSize
{
    static final long serialVersionUID = 1;
    
    private static FilterWindow INSTANCE = null;
    public static FilterWindow getInstance() { return INSTANCE; }
    
    public static void make_instance( final MainWindow parent )
    {
        if ( INSTANCE == null )
            INSTANCE = new FilterWindow( parent );
    }

    public static boolean edit( final JobData job )
    {
        boolean res = false;
        if ( INSTANCE != null )
        {
            INSTANCE.setTitle( "filter" );
            res = INSTANCE.set_show_and_get( job );        
        }
        return res;
    }
            
    private final IntInputPanel min_read_filter;
    private final IntInputPanel max_N_filter;
    private final IntInputPanel end_N_filter;
    private final SpotGroupFilterPanel spotgroup;
    private final JobReadTypePanel read_type;
            
    private boolean set_show_and_get( final JobData job )
    {
        min_read_filter.set_value( job.get_min_read_len() );
        min_read_filter.set_editable( job.get_use_min_read_len() );
        max_N_filter.set_value( job.get_max_N() );
        max_N_filter.set_editable( job.get_use_max_N() );
        end_N_filter.set_value( job.get_end_N() );
        end_N_filter.set_editable( job.get_use_end_N() );
        spotgroup.set_text( job.get_spotgroup(), job.get_full_source() );
        spotgroup.set_editable( job.get_use_spotgroup() );
        read_type.set_value( job.get_bio_read_type() );
        read_type.set_editable( job.get_use_bio_read_type() );
                
        boolean res  = show_dialog(); /* from DlgWidthMaxSize.java */
        if ( res )
        {
            job.set_min_read_len( min_read_filter.get_value() );
            job.set_use_min_read_len( min_read_filter.is_editable() );
            job.set_max_N( max_N_filter.get_value() );
            job.set_use_max_N( max_N_filter.is_editable() );
            job.set_end_N( end_N_filter.get_value() );
            job.set_use_end_N( end_N_filter.is_editable() );
            job.set_spotgroup( spotgroup.get_text() );
            job.set_use_spotgroup( spotgroup.get_editable() );
            job.set_bio_read_type( read_type.get_value() );
            job.set_use_bio_read_type( read_type.get_editable() );
            /* do not store here, because the client has to verify
               the job before actually storing it */
        }
        return res;
    }

    public FilterWindow( final MainWindow parent )
    {
        super( parent, "", new Dimension( 500, 100 ) );

        Container pane = getContentPane();
        pane.setLayout( new BoxLayout( pane, BoxLayout.PAGE_AXIS ) );
        
        min_read_filter = new IntInputPanel( "min read", "length", true, true );
        pane.add( min_read_filter );

        max_N_filter = new IntInputPanel( "max N", "count", true, true );
        pane.add( max_N_filter );

        end_N_filter = new IntInputPanel( "end N", "count", true, true );
        pane.add( end_N_filter );

        spotgroup = new SpotGroupFilterPanel( "spotgroup" );
        pane.add( spotgroup );

        read_type = new JobReadTypePanel( "read type" );
        pane.add( read_type );
        
        resize_labels( pane );
        add_save_cancel_panel( pane );
        adjust_height( 35 );
    }
}
