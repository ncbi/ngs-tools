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

import Bio.BioAccessionChecker;
import Bio.BioSpec;
import Bio.BioAccessionType;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.concurrent.CancellationException;
import java.util.concurrent.ExecutionException;
import javax.swing.BoxLayout;
import javax.swing.SwingWorker;

public class AccessionWindow extends DlgWithMaxSize
    implements ActionListener , PropertyChangeListener
{
    static final long serialVersionUID = 1;

    private static AccessionWindow INSTANCE = null;
    public static AccessionWindow getInstance() { return INSTANCE; }
    
    public static void make_instance( final MainWindow parent )
    {
        if ( INSTANCE == null )
            INSTANCE = new AccessionWindow( parent );
    }
    
    public static boolean select( final BioSpec spec )
    {
        boolean res = false;
        if ( INSTANCE != null ) res = INSTANCE.select_accession( spec );
        return res;
    }
            
    private final TextInputPanel source;
    private final TextInputPanel bio_type;
    private final TextInputPanel row_count;
    private final Save_Cancel_Panel save_cancel;
    private final BioSpec f_spec;
    private BioAccessionChecker checker;
    
    private void show_spec( final BioSpec spec )
    {
        bio_type.set_text( spec.get_type().toString() );
        row_count.set_text( String.format( "%,d rows", spec.get_count() ) );
        save_cancel.set_save_btn_status( spec.get_type() != BioAccessionType.INVALID );        
    }
    
    private void get_spec( final BioSpec spec )
    {
        if ( checker != null )
        {
            try
            {
                spec.copy( checker.get() );
            }
            catch ( InterruptedException | ExecutionException | CancellationException ex )
            {
                spec.set_type( BioAccessionType.INVALID );
            }
        }
        else
            spec.set_type( BioAccessionType.INVALID );
    }
    
    @Override public void propertyChange( final PropertyChangeEvent event )
    {
        final String propname = event.getPropertyName();
        if ( propname.equals( "state" ) )
        {
            switch ( ( SwingWorker.StateValue ) event.getNewValue() )
            {
                case DONE    : get_spec( f_spec );
                               show_spec( f_spec );
                               source.blink( false );
                               //source.set_editable( true );
                               break;
                    
                case STARTED : source.blink( true ); break;
                case PENDING : break;
            }
        }
    }
    
    private void new_check()
    {
        //source.set_editable( false );
        checker = new BioAccessionChecker( f_spec.get_accession() );
        checker.addPropertyChangeListener( this );
        checker.execute();
    }
    
    @Override public void actionPerformed( ActionEvent ae )
    {
        super.actionPerformed( ae );
        if ( isShowing() )
        {
            f_spec.clear();
            f_spec.set_accession( source.get_text() );

            show_spec( f_spec );
            save_cancel.set_save_btn_status( false );
            if ( f_spec.get_accession().length() > 8 )
            {
                if ( checker != null )
                {
                    if ( checker.isDone() )
                        new_check();
                    else
                    {
                        if ( !checker.isCancelled() )
                        {
                            checker.cancel( true );
                            new_check();
                        }
                    }
                }
                else
                    new_check();
            }
        }
    }

    private boolean select_accession( final BioSpec spec )
    {
        f_spec.clear();
        source.set_text( f_spec.get_accession() );
        //source.set_editable( true );
        show_spec( f_spec );
        boolean res  = show_dialog(); /* from DlgWidthMaxSize.java */
        if ( res )
            spec.copy( f_spec );
        else
            spec.set_type( BioAccessionType.INVALID );
        source.blink( false );
        return res;
    }
    
    /* make the job-dialog, but do not show it */
    public AccessionWindow( final MainWindow parent )
    {
        super( parent, "", new Dimension( 500, 100 ) );

        f_spec = new BioSpec();
        checker = null;
        
        Container pane = getContentPane();
        pane.setLayout( new BoxLayout( pane, BoxLayout.PAGE_AXIS ) );
        
        source = new TextInputPanel( "accession", true );
        source.add_listener( this );
        pane.add( source );

        bio_type = new TextInputPanel( "type", false );
        pane.add( bio_type );

        row_count = new TextInputPanel( "row count", false );
        pane.add( row_count );

        resize_labels( pane );
        
        save_cancel = add_save_cancel_panel( pane );
        save_cancel.set_save_btn_status( false );
        adjust_height( 35 );

    }
}
