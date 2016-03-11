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
package data;

import java.io.IOException;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.logging.*;

class MyFormatter extends Formatter
{
    private final StringBuilder sb;
    
    @Override public String format( LogRecord rec )
    {
        sb.setLength( 0 );
        sb.append( new java.util.Date() );
        sb.append( ' ' );
        sb.append( formatMessage( rec ) );
        sb.append( '\n' );                        
        return sb.toString();
    }
    
    public MyFormatter() { sb = new StringBuilder(); }
}

enum LogHandlerRequest { NONE, ENABLE, DISABLE; }

class LogTask implements Runnable
{
    private final String filename;
    private FileHandler file_handler;
    private boolean file_handler_used;
    private ConsoleHandler console_handler;
    private boolean console_handler_used;
    private Logger my_logger;
    private final boolean valid;
    private BlockingQueue<String> log_lines_queue;
    private boolean done;
    private LogHandlerRequest console_request, file_request;
    
    private void write( final String s )
    {
        if ( valid ) my_logger.log( Level.INFO, s );
    }
    
    private void enable_console_handler()
    {
        if ( !console_handler_used && console_handler != null )
        {
            my_logger.addHandler( console_handler );
            console_handler_used = true;
        }
    }

    private void disabled_console_handler()
    {
        if ( console_handler_used && console_handler != null )
        {
            my_logger.removeHandler( console_handler );
            console_handler_used = false;
        }
    }

    private void enable_file_handler()
    {
        if ( !file_handler_used && file_handler != null )
        {
            my_logger.addHandler( file_handler );
            file_handler_used = true;
        }
    }

    private void disabled_file_handler()
    {
        if ( file_handler_used && file_handler != null )
        {
            my_logger.removeHandler( file_handler );
            file_handler_used = false;
        }
    }

    @Override public void run()
    {
        try
        {
            do
            {
                write( log_lines_queue.take() );
                
                switch( console_request )
                {
                    case ENABLE  : enable_console_handler();
                                   console_request = LogHandlerRequest.NONE;
                                   break;
                    case DISABLE : disabled_console_handler();
                                   console_request = LogHandlerRequest.NONE;
                                   break;
                }

                switch( file_request )
                {
                    case ENABLE  : enable_file_handler();
                                   file_request = LogHandlerRequest.NONE;
                                   break;
                    case DISABLE : disabled_file_handler();
                                   file_request = LogHandlerRequest.NONE;
                                   break;
                }
                
            } while ( !done );
            /* empty the queue after the done flag is set... */
            while ( !log_lines_queue.isEmpty() )
            {
                write( log_lines_queue.poll() );
            }
        }
        catch ( InterruptedException iex )
        {
        }
    }
    
    public void stop() { done = true; }
    
    public void set_console_state( boolean value )
    {
        console_request = value ? LogHandlerRequest.ENABLE : LogHandlerRequest.DISABLE;
        try
        {
            log_lines_queue.put( String.format( "log.console = %s", value ? "enabled" : "disabled" ) );
        } catch ( InterruptedException ex ) { }
    }

    public void set_file_state( boolean value )
    {
        if ( file_handler == null && value )
        {
            MyFormatter f = new MyFormatter();
            try
            {
                file_handler = new FileHandler( filename, true );
                file_handler.setFormatter( f );
            }
            catch( IOException io_e )
            {
                System.out.println( "cannot create Log-Filehandler: " + io_e );
            }
        }
        
        file_request = value ? LogHandlerRequest.ENABLE : LogHandlerRequest.DISABLE;
        try
        {
            log_lines_queue.put( String.format( "log.file = %s", value ? "enabled" : "disabled" ) );
        } catch ( InterruptedException ex ) { }
    }
    
    public LogTask( final String filename,
                    final String info,
                    BlockingQueue<String> q,
                    final boolean log_to_file,
                    final boolean log_to_cons )
    {
        this.filename = filename;
        my_logger = null;
        log_lines_queue = q;
        boolean success = false;

        file_handler = null;
        console_handler = null;
        file_handler_used = log_to_file;
        console_handler_used = log_to_cons;
        console_request = LogHandlerRequest.NONE;
        file_request = LogHandlerRequest.NONE;
                
        try
        {
            MyFormatter f = new MyFormatter();

            my_logger = Logger.getLogger( info );
            my_logger.setLevel( Level.INFO );

            if ( filename != null && !filename.isEmpty() && file_handler_used )
            {
                file_handler = new FileHandler( filename, true );
                file_handler.setFormatter( f );
                my_logger.addHandler( file_handler );
            }
           
            console_handler = new ConsoleHandler();
            console_handler.setFormatter( f );
            if ( console_handler_used )
                my_logger.addHandler( console_handler );
        
            LogManager my_log_manager = LogManager.getLogManager();
            my_log_manager.addLogger( my_logger );
            
            success = true;
        }   
        catch( IOException io_e )
        {
            System.out.println( "cannot create Log-Filehandler: " + io_e );
        }
        valid = success;
        done = false;
    }
}
        
public class CLogger
{
    private final BlockingQueue<String> log_lines_queue;
    private LogTask my_task;
    private Thread my_thread;

    private static CLogger INSTANCE = null;
    public static CLogger getInstance()
    {
        if ( INSTANCE == null ) INSTANCE = new CLogger();
        return INSTANCE;
    }

    /**
     * the main functionality: log this string...
     * 
     * @param s ... string to write to log-file ( can be null )
     */
    public void log_str( final String s )
    {
        if ( s != null && !s.isEmpty() )
        {
            try
            {
                log_lines_queue.put( s );
            }
            catch( InterruptedException iex )
            {
                // we simply do nothing then...
            }
        }
    }

    public static void log( final String s )
    {
        if ( INSTANCE != null )
            INSTANCE.log_str( s );
    }

    public static void logfmt( String fmt, Object... arguments )
    {
        if ( INSTANCE != null )
            INSTANCE.log_str( String.format( fmt, arguments ) );
    }
    
    private void start_log( final String filename, final String info,
            final boolean log_to_file, final boolean log_to_cons )
    {
        my_task = new LogTask( filename, info, log_lines_queue, log_to_file, log_to_cons );
        my_thread = new Thread( my_task );
        my_thread.start();
    }
            
    /**
     * starts the logging task ( constructor does not! )
     * @param filename .... file to write log into
     * @param info ........ info string to be written into log
     * 
     */
    public static void start( final String filename, final String info,
            final boolean log_to_file, final boolean log_to_cons )
    {
        if ( INSTANCE != null )
            INSTANCE.start_log( filename, info, log_to_file, log_to_cons );
    }
    
    /**
     * stops the logging task
     * 
     * @param s ... optional final string to write
     */
    public static void stop( final String s )
    {
        if ( INSTANCE != null )
        {
            if ( INSTANCE.my_task != null )
            {
                if ( s != null ) log( s );
                INSTANCE.my_task.stop();
            }
        }
    }

    public static void set_file_logging( boolean value )
    {
        if ( INSTANCE != null )
        {
            if ( INSTANCE.my_task != null )
                INSTANCE.my_task.set_file_state( value );
        }
    }

    public static void set_console_logging( boolean value )
    {
        if ( INSTANCE != null )
        {
            if ( INSTANCE.my_task != null )
                INSTANCE.my_task.set_console_state( value );
        }
    }

    /**
     * custom logger
     */
    public CLogger()
    {
        Logger global_logger = Logger.getGlobal();
        Handler[] h = global_logger.getParent().getHandlers();
        for ( int i = 0; i < h.length; ++i )
            global_logger.getParent().removeHandler( h[ i ] );
        
        log_lines_queue = new LinkedBlockingQueue<String>();
        my_task = null;
        my_thread = null;
    }
}
