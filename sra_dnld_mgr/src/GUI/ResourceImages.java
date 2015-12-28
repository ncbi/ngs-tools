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

import java.awt.event.ActionListener;
import java.net.URL;
import javax.swing.BorderFactory;
import javax.swing.ImageIcon;
import javax.swing.JButton;

public class ResourceImages
{
    private static final ResourceImages INSTANCE = new ResourceImages();
    public static ResourceImages getInstance() { return INSTANCE; }
    
    static final String START = "Play";
    static final String PAUSE = "Pause";
    static final String STOP = "Stop";
    static final String SETTINGS = "Settings";
    static final String EXIT = "Exit";
    static final String ADD = "Add";
    static final String CHECK = "Check";
    static final String DEL = "Delete";
    static final String CLEAR = "Clear";
    static final String LOGO = "sra_logo";
    static final String FILTER = "filter";
    
    private final ImageIcon StartImage;
    private final ImageIcon PauseImage;
    private final ImageIcon StopImage;
    private final ImageIcon SettingsImage;
    private final ImageIcon ExitImage;
    private final ImageIcon AddImage;
    private final ImageIcon CheckImage;
    private final ImageIcon DelImage;
    private final ImageIcon ClearImage;
    private final ImageIcon LogoImage;
    private final ImageIcon FilterImage;
    
    private ImageIcon make_img( String key )
    {
        ImageIcon res = null;
        String path = String.format( "resources/%s.png", key );
        URL url = getClass().getClassLoader().getResource( path );
        if ( url != null ) res = new ImageIcon( url );
        return res;
    }
    
    public static ImageIcon get_start_img() { return INSTANCE.StartImage; }
    public static ImageIcon get_pause_img() { return INSTANCE.PauseImage; }
    public static ImageIcon get_stop_img() { return INSTANCE.StopImage; }
    public static ImageIcon get_settings_img() { return INSTANCE.SettingsImage; }
    public static ImageIcon get_exit_img() { return INSTANCE.ExitImage; }
    public static ImageIcon get_add_img() { return INSTANCE.AddImage; }
    public static ImageIcon get_check_img() { return INSTANCE.CheckImage; }
    public static ImageIcon get_del_img() { return INSTANCE.DelImage; }
    public static ImageIcon get_clear_img() { return INSTANCE.ClearImage; }
    public static ImageIcon get_logo_img() { return INSTANCE.LogoImage; }
    public static ImageIcon get_filter_img() { return INSTANCE.FilterImage; }
    
    static public JButton make_img_button( ImageIcon img,
            String txt, int border, ActionListener action )
    {
        JButton res;
        if ( img != null )
            res = new JButton( img );
        else
            res = new JButton( txt );

        res.setBorder( BorderFactory.createEmptyBorder( border, border, border, border ) );
        res.setContentAreaFilled( false );
        res.setFocusPainted( false );
        res.setToolTipText( txt );
        if ( action != null ) res.addActionListener( action );
        return res;
    }
            
    public ResourceImages()
    {
        StartImage = make_img( START );
        PauseImage = make_img( PAUSE );
        StopImage = make_img( STOP );
        SettingsImage = make_img( SETTINGS );
        ExitImage = make_img( EXIT );
        AddImage = make_img( ADD );
        CheckImage = make_img( CHECK );
        DelImage = make_img( DEL );
        ClearImage = make_img( CLEAR );
        LogoImage = make_img( LOGO );
        FilterImage = make_img( FILTER );
    }
}
