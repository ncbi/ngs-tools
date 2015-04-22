/*===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
*  Please cite the author in any work or product based on this material.
*
* ===========================================================================
*
*/

import java.util.*;
import java.io.*;

enum opt_sel
{
	OUT, TEMP, SQL, MASTER, SLICES, FORMAT, EXECMEM, DRVMEM, MAXRES, LOCNODES, RESREQ, NONE;

	public static opt_sel from_string( final String s )
	{
		try
		{
			return Enum.valueOf( opt_sel.class, s.trim().toUpperCase() );
		}
        catch ( IllegalArgumentException e )
        {
        }
		return NONE;
	}
}

public class NGS_SQL_Opt
{
	private List< NGS_TAB_Opt > tab_options;
	private List< String > invalid_lines;
	
    private String out, temp, sql, master, exec_mem, drv_mem, max_res;
    private int slices, local_nodes, res_requests;
	private out_fmt format;
	
    public NGS_SQL_Opt( final String filename )
    {
		master			= "";
        slices 			= 8;
		local_nodes		= 4;
		res_requests	= 4;
		tab_options 	= new ArrayList< NGS_TAB_Opt >();
		invalid_lines 	= new ArrayList< String >();
		out				= "result.txt";
		temp			= "result.temp";
		format 			= out_fmt.TAB;
		exec_mem		= "";
		drv_mem			= "";
		max_res			= "";		
		read_option_file( filename );
    }

	public NGS_TAB_Opt get( int idx )
	{
		NGS_TAB_Opt res = null;
		if ( idx < tab_options.size() ) res = tab_options.get( idx );
		return res;
	}
	
	private void define_line( final String key, final String value )
	{
		NGS_TAB_Opt tab_opt = new NGS_TAB_Opt( tab_options.size(), key, value );
		if ( tab_opt.valid() )
			tab_options.add( tab_opt );
		else
			invalid_lines.add( key + " = " + value );
	}
	
	private boolean read_option_file( final String filename )
	{
		try
		{
			BufferedReader buf_reader = new BufferedReader( new InputStreamReader( new FileInputStream( filename ) ) );
			String line;
			do
			{
				line = buf_reader.readLine();
				if ( line != null && !line.startsWith( "#" ) )
				{
					int equalsign = line.indexOf( '=' );
					if ( equalsign > 0 )
					{
						String key = line.substring( 0, equalsign ).trim();
						String value = line.substring( equalsign + 1 ).trim();
						switch( opt_sel.from_string( key ) )
						{
							case OUT		: out = value; break;
							case SQL		: sql = value; break;
							case MASTER		: master = value; break;
							case TEMP		: temp = value; break;
							case EXECMEM	: exec_mem = value; break;
							case DRVMEM		: drv_mem = value; break;
							case MAXRES		: max_res = value; break;
							case SLICES		: slices = Integer.parseInt( value ); break;
							case LOCNODES	: local_nodes = Integer.parseInt( value ); break;
							case RESREQ		: res_requests = Integer.parseInt( value ); break;
							case FORMAT		: format = out_fmt.from_string( value ); break;
							case NONE		: define_line( key, value ); break;
						}
					}
				}
			} while ( line != null );
		}
        catch ( IOException e )
        {
            System.out.printf( "cannot read file '%s'%n", filename );
        }
	
		return ( invalid_lines.size() == 0 );
	}
	
	public boolean valid()
	{
		boolean res = ( out != null && temp != null && sql != null && invalid_lines.size() == 0 );
		for ( NGS_TAB_Opt o : tab_options ) res = ( res && o.valid() );
		return res;
	}

	public String get_Master()
	{
		String res = master;
		if ( res.length() < 1 && local_nodes > 0 )
			res = String.format( "local[%d]", local_nodes );
		return res;
	}
	
	public String get_Out() 		{ return out; }
	public String get_Temp() 		{ return temp; }
	public String get_Sql() 		{ return sql; }
	public String get_Exec_Mem()	{ return exec_mem; }
	public String get_Drv_Mem()		{ return drv_mem; }
	public String get_Max_Res()		{ return max_res; }
	
	public int get_table_count() 	{ return tab_options.size(); }
    public int get_Slices() 		{ return slices; }
	public int get_LocalNodes() 	{ return local_nodes; }
	public int get_ResRequests() 	{ return res_requests; }
	
	public out_fmt get_format()	{ return format; }
	
	public void report()
	{
		System.out.printf( "OUT      =%s%n", get_Out() );
		System.out.printf( "TEMP     =%s%n", get_Temp() );
		System.out.printf( "SQL      =%s%n", get_Sql() );
		System.out.printf( "MASTER   =%s%n", get_Master() );
		System.out.printf( "FORMAT   =%s%n", format.toString() );
		System.out.printf( "EXECMEM  =%s%n", get_Exec_Mem() );
		System.out.printf( "DRVMEM   =%s%n", get_Drv_Mem() );
		System.out.printf( "MAXRES   =%s%n", get_Max_Res() );

		System.out.printf( "SLICES   =%d%n", get_Slices() );
		System.out.printf( "LOCNODES =%d%n", get_LocalNodes() );
		System.out.printf( "RESREQ   =%d%n", get_ResRequests() );
		
		for ( NGS_TAB_Opt o : tab_options ) o.report();
		
		if ( invalid_lines.size() > 0 )
		{
			System.out.printf( "%nINVALID:%n" );
			for ( String s : invalid_lines ) System.out.printf( "%s%n", s );
		}
	}
}
