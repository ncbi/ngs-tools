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
package RobotScript;

import java.awt.AWTException;
import java.awt.MouseInfo;
import java.awt.Point;
import java.awt.Robot;
import java.awt.event.InputEvent;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

public class RobotScript
{
    private Robot robot;
    
    void delay( final String words[] )
    {
        int time = Integer.parseInt( words[ 1 ] );
        if ( time > 0 ) robot.delay( time );
    }
    
    void moveto( final String words[] )
    {
        if ( words.length > 2 )
        {
            int x = Integer.parseInt( words[ 1 ] );
            int y = Integer.parseInt( words[ 2 ] );
            robot.mouseMove( x, y );
        }
    }

    void moveby( final String words[] )
    {
        if ( words.length > 2 )
        {
            int dx = Integer.parseInt( words[ 1 ] );
            int dy = Integer.parseInt( words[ 2 ] );
            Point p = MouseInfo.getPointerInfo().getLocation();
            robot.mouseMove( p.x + dx, p.y + dy );
        }
    }

    void click( final String words[] )
    {
        String how = words[ 1 ];
        int btn = InputEvent.BUTTON1_MASK;
        if ( how.startsWith( "left" ) )
            btn = InputEvent.BUTTON1_MASK;
        else if ( how.startsWith( "middle" ) )
            btn = InputEvent.BUTTON2_MASK;
        else if ( how.startsWith( "right" ) )
            btn = InputEvent.BUTTON3_MASK;
        
        robot.mousePress( btn );
        robot.delay( 200 );
        robot.mouseRelease( btn );
        robot.delay( 200 );
    }

    void type( int code )
    {
        robot.delay( 40 );
        robot.keyPress( code );
        robot.keyRelease( code );
    }
    
    void type( final String s )
    {
        byte[] bytes = s.getBytes();
        for ( byte b : bytes )
        {
            int code = b;
            // keycode only handles [A-Z] (which is ASCII decimal [65-90])
            if (code > 96 && code < 123) code = code - 32;
            type( code );
        }        
    }

    void type( final String words[] )
    {
        for ( int i = 1; i < words.length; ++i )
            type( words[ i ] );
    }

    void command( final String words[] )
    {
        String cmd = words[ 0 ];
        if ( cmd.startsWith( "delay" ) )
        {
            delay( words );
        }
        else if ( cmd.startsWith( "moveto" ) )
        {
            moveto( words );
        }
        else if ( cmd.startsWith( "moveby" ) )
        {
            moveby( words );
        }
        else if ( cmd.startsWith( "click" ) )
        {
            click( words );
        }
        else if ( cmd.startsWith( "type" ) )
        {
            type( words );
        }
    }
    
    void play( final String filename )
    {
        if ( robot != null )
        {
            try
            {
                BufferedReader br = new BufferedReader( new FileReader( filename ) );
                for ( String line; ( line = br.readLine() ) != null; )
                {
                    if ( !line.isEmpty() && !line.startsWith( "#" ) )
                    {
                        String words[] = line.split( " " );
                        if ( words.length > 1 )
                            command( words );
                    }
                }
            }
            catch ( IOException ex ) { ; }
        }
    }
    
    RobotScript()
    {
        try
        {
            robot = new Robot();
            robot.setAutoDelay( 40 );
            robot.setAutoWaitForIdle( true );
        }
        catch ( AWTException ex )
        {
            robot = null;
        }
    }
}
