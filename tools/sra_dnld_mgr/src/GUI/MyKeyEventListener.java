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

import java.awt.KeyEventPostProcessor;
import java.awt.KeyboardFocusManager;
import java.awt.event.KeyEvent;

public class MyKeyEventListener implements KeyEventPostProcessor
{
    private static final MyKeyEventListener INSTANCE = new MyKeyEventListener();
    
    public static void set_receiver( final MyKeyEventReceiver receiver )
    {
        if ( INSTANCE != null ) INSTANCE.receiver = receiver;
    }
            
    private MyKeyEventReceiver receiver;
    
    private MyKeyEventListener()
    {
        this.receiver = null;
        KeyboardFocusManager kfm = KeyboardFocusManager.getCurrentKeyboardFocusManager();
        kfm.addKeyEventPostProcessor( this );
    }

    @Override public boolean postProcessKeyEvent( KeyEvent e )
    {
        if ( e.getID() == KeyEvent.KEY_PRESSED && receiver != null )
            return receiver.on_key( e );
        return false;
    }

}
