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

public class JobWindow extends DlgWithMaxSize
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
    
    private JobData current_job;
    
    private final ActionListener on_preview = new on_preview_event();
    private final ActionListener on_filter = new on_filter_event();
    private final ItemListener on_format = new on_format_event();
    private final ActionListener on_references = new on_references_event();
    
    private final TextInputPanel source = new TextInputPanel( "source", false );
    private final TextInputPanel bio_type = new TextInputPanel( "type", false );
    private final PathChooserPanel export_path = new PathChooserPanel( "export path" );
    private final TextInputPanel status = new TextInputPanel( "state", false );
    private final JobFormatPanel format = new JobFormatPanel( "format", on_format );
    private final LineEndingsPanel line_endings = new LineEndingsPanel( "line ending" );
    private final IntInputPanel line_wrap = new IntInputPanel( "line wrap", "chars", true, true );
    private final IntInputPanel fixed_qual = new IntInputPanel( "fixed quality", "pthread", false, true );
    private final BoolSettingsPanel fastq_dump_name = new BoolSettingsPanel( "fastq-dump-style", false );
    private final TextInputPanel md5 = new TextInputPanel( "md5-sum", false );
    private final TextInputPanel rejected = new TextInputPanel( "rejected", false );
    private final ButtonPanel preview =  new ButtonPanel( "preview", "preview FASTQ/FASTA",
                    on_preview, ResourceImages.get_preview_img() );
    private final ButtonPanel filter = new ButtonPanel( "preview", "filter FASTQ/FASTA",
                    on_filter, ResourceImages.get_filter_img() );
    private final ButtonPanel refs = new ButtonPanel( "references", "show used references",
                    on_references, null );

    private class on_preview_event implements ActionListener
    {
        @Override public void actionPerformed( ActionEvent ae ) { show_preview(); }
    }

    private class on_filter_event implements ActionListener
    {
        @Override public void actionPerformed( ActionEvent ae ) { show_filter_window(); }
    }

    private class on_format_event implements ItemListener
    {
        @Override  public void itemStateChanged( ItemEvent ie )
        {
            current_job.set_format( format.get_format() );
            setup_by_format();
        }
    }

    private class on_references_event implements ActionListener
    {
        @Override public void actionPerformed( ActionEvent ae ) { show_references(); }
    }

    private void setup_by_format()
    {
        boolean is_dnld = current_job.is_dnld();
        
        format.set_sub_format_enabled( !is_dnld );
        line_endings.set_enabled( !is_dnld );
        line_wrap.set_enabled( !is_dnld );
        fixed_qual.set_enabled( !is_dnld );
        fastq_dump_name.set_enabled( !is_dnld );
        preview.set_enabled( !is_dnld );
        filter.set_enabled( !is_dnld );
        export_path.set_text( is_dnld ? 
                                current_job.get_downloadpath() :
                                current_job.get_exportpath() );
        refs.set_enabled( current_job.is_csra() );
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
    
    /* transfer the values from the job into the dialog-fields
       and show the dialog, if OK pressed transfer the values
       back into the job-instance */
    private boolean edit_job( JobData job )
    {
        this.setTitle( job.get_short_source() );
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

    public void show_references()
    {
        if ( current_job != null )
        {
            extract_values( current_job );
            MyKeyEventListener.set_receiver( null );
            this.setVisible( false );
            ReferenceWindow.show_refs( current_job );
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
        
        pane.add( source );
        pane.add( bio_type );
        pane.add( export_path );
        pane.add( status );
        pane.add( format );
        pane.add( line_endings );
        pane.add( line_wrap );
        pane.add( fixed_qual );
        pane.add( fastq_dump_name );
        pane.add( md5 );
        pane.add( rejected );
        pane.add( preview );
        pane.add( filter );
        pane.add( refs );
        
        resize_labels( pane );
        add_save_cancel_panel( pane );
        adjust_height( 35 );
    }
}
