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

import data.CLogger;
import data.Settings;
import java.awt.*;
import javax.swing.*;

public class SettingsWindow extends DlgWithMaxSize
{
    static final long serialVersionUID = 1;

    private static SettingsWindow INSTANCE = null;
    public static SettingsWindow getInstance() { return INSTANCE; }
    
    public static void make_instance( final MainWindow parent )
    {
        if ( INSTANCE == null )
            INSTANCE = new SettingsWindow( parent );
    }
    
    public static boolean edit()
    {
        boolean res = false;
        if ( INSTANCE != null ) res = INSTANCE.set_and_get();
        return res;
    }
    
    private final PathChooserPanel job_path = new PathChooserPanel( "job path" );
    private final PathChooserPanel export_path = new PathChooserPanel( "export path" );
    private final IntSettingsPanel max_dnld = new IntSettingsPanel( "max. dnld", 1, 5 );
    private final BoolSettingsPanel autostart = new BoolSettingsPanel( "autoresume", false );
    private final BoolSettingsPanel log_to_file = new BoolSettingsPanel( "log to file", false );
    private final BoolSettingsPanel log_to_cons = new BoolSettingsPanel( "log to console", false );
    private final LineEndingsPanel line_endings = new LineEndingsPanel( "line ending" );
    private final IntInputPanel line_wrap = new IntInputPanel( "line wrap", "chars", true, true );
    private final IntInputPanel fixed_qual = new IntInputPanel( "fixed qualitiy", "pthread", false, true );
    private final IntInputPanel preview_rows = new IntInputPanel( "preview rows", "rows", true, false );
    private final IntInputPanel conn_timeout = new IntInputPanel( "conn. timeout", "ms", true, false );
    private final IntInputPanel read_timeout = new IntInputPanel( "read timeout", "ms", true, false );
    
    private boolean set_and_get()
    {
        Settings settings = Settings.getInstance();
        job_path.set_text( settings.get_jobpath() );
        export_path.set_text( settings.get_exportpath() );
        max_dnld.set_value( settings.get_maxdownloads() );
        autostart.set_value( settings.get_autostart() );
        log_to_file.set_value( settings.get_log_to_file() );
        log_to_cons.set_value( settings.get_log_to_cons() );
        line_endings.set_value( settings.get_line_ending() );
        line_wrap.set_value( settings.get_line_wrap() );
        line_wrap.set_editable( settings.get_use_line_wrap() );
        fixed_qual.set_value( settings.get_fixed_qual() );
        fixed_qual.set_editable( settings.get_use_fixed_qual() );
        preview_rows.set_value( settings.get_preview_rows() );
        conn_timeout.set_value( settings.get_conn_timeout() );
        read_timeout.set_value( settings.get_read_timeout() );
        
        boolean res = show_dialog();
        if ( res )
        {
            settings.set_jobpath( job_path.get_text() );
            settings.set_exportpath( export_path.get_text() );
            settings.set_maxdownloads( max_dnld.get_value() );
            settings.set_autostart( autostart.get_value() );

            settings.set_log_to_file( log_to_file.get_value() );
            if ( log_to_file.has_changed() )
                CLogger.set_file_logging( log_to_file.get_value() );
            
            settings.set_log_to_cons( log_to_cons.get_value() );
            if ( log_to_cons.has_changed() )
                CLogger.set_console_logging( log_to_cons.get_value() );
            
            settings.set_line_ending(  line_endings.get_value() );
            settings.set_line_wrap( line_wrap.get_value() );
            settings.set_use_line_wrap( line_wrap.is_editable() );
            settings.set_fixed_qual( fixed_qual.get_value() );
            settings.set_use_fixed_qual( fixed_qual.is_editable() );
            settings.set_preview_rows( preview_rows.get_value() );
            settings.set_conn_timeout( conn_timeout.get_value() );
            settings.set_read_timeout( read_timeout.get_value() );
            
            settings.store();
        }
        return res;
    }
    
    public SettingsWindow( final MainWindow parent )
    {
        super( parent, "global settings", new Dimension( 500, 100 ) );
       
        Container pane = getContentPane();
        pane.setLayout( new BoxLayout( pane, BoxLayout.PAGE_AXIS ) );

        pane.add( job_path );
        pane.add( export_path );
        pane.add( max_dnld );
        pane.add( autostart );
        pane.add( log_to_file );
        pane.add( log_to_cons );
        pane.add( line_endings );
        pane.add( line_wrap );
        pane.add( fixed_qual );
        pane.add( preview_rows );
        pane.add( conn_timeout );
        pane.add( read_timeout );

        resize_labels( pane );
        add_save_cancel_panel( pane );
        
        adjust_height( 35 );
    }
}
