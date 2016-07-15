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

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

class Action_With_Parent implements ActionListener
{
    MainWindow parent;
    public Action_With_Parent( MainWindow parent ) { super(); this.parent = parent; }
    @Override public void actionPerformed( ActionEvent e ) {  }
}

class ExitAction extends Action_With_Parent
{
    public ExitAction( MainWindow parent ) { super( parent ); }
    @Override public void actionPerformed( ActionEvent e ) { parent.on_exit( 0 ); }
}

class SettingsAction extends Action_With_Parent
{
    public SettingsAction( MainWindow parent ) { super( parent ); }
    @Override public void actionPerformed( ActionEvent e ) { SettingsWindow.edit();  }
}

class AddJobAction extends Action_With_Parent
{
    public AddJobAction( MainWindow parent ) { super( parent ); }
    @Override public void actionPerformed( ActionEvent e ) { parent.NewJob(); }
}

/* ------------------------------------------------------ */

class JobsAdd extends JMenuItem
{
    static final long serialVersionUID = 1;
    
    public JobsAdd( MainWindow parent )
    {
        super( "Add Job" );
        setMnemonic( KeyEvent.VK_A );
        setAccelerator( KeyStroke.getKeyStroke( KeyEvent.VK_A, 0 ) );
        setToolTipText( "add a new job" );
        addActionListener( new AddJobAction( parent ) ); 
    }
}

class JobsDel extends JMenuItem
{
    static final long serialVersionUID = 1;
    
    public JobsDel( MainWindow parent )
    {
        super( "Delete Job" );
        setMnemonic( KeyEvent.VK_D );
        setAccelerator( KeyStroke.getKeyStroke( KeyEvent.VK_D, 0 ) );
        setToolTipText( "delete the selected job" );
    }
}

class JobsMenu extends JMenu
{
    static final long serialVersionUID = 1;
    
    public JobsMenu( MainWindow parent )
    {
        super( "Jobs" );
        setMnemonic( KeyEvent.VK_J );
        add( new JobsAdd( parent ) );
        add( new JobsDel( parent ) );
        getAccessibleContext().setAccessibleDescription( "Operations on jobs" );
        setToolTipText( "Operations on jobs" );
    }
}

class MainMenuSettings extends JMenuItem
{
    static final long serialVersionUID = 1;
    
    public MainMenuSettings( MainWindow parent )
    {
        super( "Settings" );
        setMnemonic( KeyEvent.VK_S );
        setAccelerator( KeyStroke.getKeyStroke( KeyEvent.VK_S, 0 ) );
        addActionListener( new SettingsAction( parent ) ); 
        setToolTipText( "global settings" );
    }
}

class MainMenuExit extends JMenuItem
{
    static final long serialVersionUID = 1;
    
    public MainMenuExit( MainWindow parent )
    {
        super( "Exit" );
        setMnemonic( KeyEvent.VK_X );
        setAccelerator( KeyStroke.getKeyStroke( KeyEvent.VK_X, 0 ) );
        addActionListener( new ExitAction( parent ) );
        setToolTipText( "exit the tool" );
    }
}

class MainMenu extends JMenu
{
    static final long serialVersionUID = 1;
    
    public MainMenu( MainWindow parent )
    {
        super( "Main" );
        setMnemonic( KeyEvent.VK_M );
        add( new MainMenuSettings( parent ) );
        addSeparator();
        add( new MainMenuExit( parent ) );
        getAccessibleContext().setAccessibleDescription( "The main menu" );
        setToolTipText( "The main menu" );
    }
}

public class TheMenuBar extends JMenuBar
{
    static final long serialVersionUID = 1;
    
    public TheMenuBar( MainWindow parent )
    {
        super();
        setOpaque( true );
        setBackground( new Color( 154, 165, 127 ) );
        setPreferredSize( new Dimension( 400, 25 ) );
        add( new MainMenu( parent ) );
        add( new JobsMenu( parent ) );
    }
}
