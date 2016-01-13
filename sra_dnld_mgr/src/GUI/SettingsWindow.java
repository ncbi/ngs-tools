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
    
    private final PathChooserPanel job_path;
    private final PathChooserPanel export_path;
    private final IntSettingsPanel max_dnld;
    private final BoolSettingsPanel autostart;
    private final LineEndingsPanel line_endings;
    private final IntInputPanel line_wrap;
    private final IntInputPanel fixed_qual;
    private final IntInputPanel preview_rows;
    
    private boolean set_and_get()
    {
        Settings settings = Settings.getInstance();
        job_path.set_text( settings.get_jobpath() );
        export_path.set_text( settings.get_exportpath() );
        max_dnld.set_value( settings.get_maxdownloads() );
        autostart.set_value( settings.get_autostart() );
        line_endings.set_value( settings.get_line_ending() );
        line_wrap.set_value( settings.get_line_wrap() );
        line_wrap.set_editable( settings.get_use_line_wrap() );
        fixed_qual.set_value( settings.get_fixed_qual() );
        fixed_qual.set_editable( settings.get_use_fixed_qual() );
        preview_rows.set_value( settings.get_preview_rows() );
        boolean res = show_dialog();
        if ( res )
        {
            settings.set_jobpath( job_path.get_text() );
            settings.set_exportpath( export_path.get_text() );
            settings.set_maxdownloads( max_dnld.get_value() );
            settings.set_autostart( autostart.get_value() );
            settings.set_line_ending(  line_endings.get_value() );
            settings.set_line_wrap( line_wrap.get_value() );
            settings.set_use_line_wrap( line_wrap.is_editable() );
            settings.set_fixed_qual( fixed_qual.get_value() );
            settings.set_use_fixed_qual( fixed_qual.is_editable() );
            settings.set_preview_rows( preview_rows.get_value() );
            settings.store();
        }
        return res;
    }
    
    public SettingsWindow( final MainWindow parent )
    {
        super( parent, "global settings", new Dimension( 500, 100 ) );
       
        Container pane = getContentPane();
        pane.setLayout( new BoxLayout( pane, BoxLayout.PAGE_AXIS ) );

        job_path = new PathChooserPanel( "job path" );
        pane.add( job_path );

        export_path = new PathChooserPanel( "export path" );
        pane.add( export_path );
        
        max_dnld = new IntSettingsPanel( "max. dnld", 1, 5 );
        pane.add( max_dnld );
        
        autostart = new BoolSettingsPanel( "autoresume", false );
        pane.add( autostart );
        
        line_endings = new LineEndingsPanel( "line ending" );
        pane.add( line_endings );

        line_wrap = new IntInputPanel( "line wrap", "chars", true, true );
        pane.add( line_wrap );

        fixed_qual = new IntInputPanel( "fixed qualitiy", "pthread", false, true );
        pane.add( fixed_qual );

        preview_rows = new IntInputPanel( "preview rows", "rows", true, false );
        pane.add( preview_rows );

        resize_labels( pane );
        add_save_cancel_panel( pane );
        
        adjust_height( 35 );
    }
}
