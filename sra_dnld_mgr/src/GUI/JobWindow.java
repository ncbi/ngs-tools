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
import java.awt.event.*;
import javax.swing.*;
import job.JobFormat;

public class JobWindow extends DlgWithMaxSize
    implements ItemListener
{
    private static JobWindow INSTANCE = null;
    public static JobWindow getInstance() { return INSTANCE; }
    
    public static void make_instance( final MainWindow parent )
    {
        if ( INSTANCE == null )
            INSTANCE = new JobWindow( parent );
    }
    
    public static boolean edit( JobData job )
    {
        boolean res = false;
        if ( INSTANCE != null ) res = INSTANCE.edit_job( job );
        return res;
    }
    
    private final TextInputPanel source;
    private final TextInputPanel bio_type;
    private final PathChooserPanel export_path;
    private final TextInputPanel status;
    private final JobFormatPanel format;
    private final LineEndingsPanel line_endings;
    private final IntInputPanel line_wrap;
    private final IntInputPanel fixed_qual;
    private final BoolSettingsPanel fastq_dump_name;
    private final TextInputPanel md5;
    private final TextInputPanel rejected;
    private final Save_Cancel_Filter_Panel save_cancel_filter;
    
    private JobData current_job;
    
    private void setup_by_format()
    {
        boolean is_dnld = current_job.is_dnld();
        
        format.set_sub_format_enabled( !is_dnld );
        line_endings.set_enabled( !is_dnld );
        line_wrap.set_enabled( !is_dnld );
        fixed_qual.set_enabled( !is_dnld );
        fastq_dump_name.set_enabled( !is_dnld );
        save_cancel_filter.set_filter_btn_enabled( !is_dnld );
        save_cancel_filter.set_preview_btn_enabled( !is_dnld );
        export_path.set_text( is_dnld ? 
                                current_job.get_downloadpath() :
                                current_job.get_exportpath() );
    }
    
    @Override public void itemStateChanged( ItemEvent e )
    {
        current_job.set_format( format.get_format() );
        setup_by_format();
    }

    private void populate_window()
    {
        source.set_text( current_job.get_full_source() );
        bio_type.set_text( current_job.get_bio_type().toString() );
        status.set_text( current_job.get_state().toString() );
        format.set_format( current_job.get_format() );
        format.set_sub_format( current_job.get_subformat() );
        line_endings.set_value( current_job.get_line_ending() );
        line_wrap.set_value( current_job.get_line_wrap() );
        line_wrap.set_editable( current_job.get_use_line_wrap() );
        fixed_qual.set_value( current_job.get_fixed_qual() );
        fixed_qual.set_editable( current_job.get_use_fixed_qual() );
        fastq_dump_name.set_value( current_job.get_fastq_dump_name() );
        md5.set_text( current_job.get_md5() );
        rejected.set_text( String.format( "%d" , current_job.get_rejected() ) );
        setup_by_format();
    }

    private void extract_values( final JobData job )
    {
        if ( current_job.is_dnld() )
            job.set_downloadpath( export_path.get_text() );
        else
            job.set_exportpath( export_path.get_text() );
        
        job.set_format( format.get_format() );
        job.set_subformat( format.get_sub_format() );
        job.set_line_ending( line_endings.get_value() );
        job.set_use_line_wrap( line_wrap.is_editable() );
        job.set_line_wrap( line_wrap.get_value() );
        job.set_use_fixed_qual( fixed_qual.is_editable() );
        job.set_fixed_qual( fixed_qual.get_value() );
        job.set_fastq_dump_name( fastq_dump_name.get_value() );
    }
    
    private boolean set_show_and_get( final JobData job )
    {
        current_job = job;
        populate_window();
        boolean res  = show_dialog(); /* from DlgWidthMaxSize.java */
        if ( res )
        {
            extract_values( job );
            /* do not store here, because the client has to verify
               the job before actually storing it */
        }
        current_job = null;
        return res;
    }
    
    /* transfer the values from the job into the dialog-fields
       and show the dialog, if OK pressed transfer the values
       back into the job-instance */
    private boolean edit_job( JobData job )
    {
        this.setTitle( job.get_short_source() );
        return set_show_and_get( job );        
    }

    @Override public void on_save_cancel_filter_event( final SaveCancelFilterEventType event_type )
    {
        super.on_save_cancel_filter_event( event_type );
        switch( event_type )
        {
            case FILTER     : show_filter_window(); break;
            case PREVIEW    : show_preview(); break;
        }
    }
    
    public void show_filter_window()
    {
        if ( current_job != null )
        {
            MyKeyEventListener.set_receiver( null );
            this.setVisible( false );
            FilterWindow.edit( current_job );
            MyKeyEventListener.set_receiver( this );
            this.setVisible( true );
        }
    }

    public void show_preview()
    {
        if ( current_job != null )
        {
            extract_values( current_job );
            MyKeyEventListener.set_receiver( null );
            this.setVisible( false );
            PreviewWindow.show_preview( current_job );
            MyKeyEventListener.set_receiver( this );
            this.setVisible( true );
        }
    }

    /* make the job-dialog, but do not show it */
    public JobWindow( final MainWindow parent )
    {
        super( parent, "", new Dimension( 500, 100 ) );
        current_job = null;
        
        Container pane = getContentPane();
        pane.setLayout( new BoxLayout( pane, BoxLayout.PAGE_AXIS ) );
        
        source = new TextInputPanel( "source", false );
        pane.add( source );
        
        bio_type = new TextInputPanel( "type", false );
        pane.add( bio_type );
        
        export_path = new PathChooserPanel( "export path" );
        pane.add( export_path );
        
        status = new TextInputPanel( "state", false );
        pane.add( status );

        format = new JobFormatPanel( "format", this );
        pane.add( format );

        line_endings = new LineEndingsPanel( "line ending" );
        pane.add( line_endings );

        line_wrap = new IntInputPanel( "line wrap", "chars", true, true );
        pane.add( line_wrap );

        fixed_qual = new IntInputPanel( "fixed quality", "pthread", false, true );
        pane.add( fixed_qual );

        fastq_dump_name = new BoolSettingsPanel( "fastq-dump-style", false );
        pane.add( fastq_dump_name );
        
        md5 = new TextInputPanel( "md5-sum", false );
        pane.add( md5 );

        rejected = new TextInputPanel( "rejected", false );
        pane.add( rejected );

        resize_labels( pane );
        
        save_cancel_filter = add_save_cancel_filter_panel( pane );
        
        adjust_height( 35 );
    }
}
