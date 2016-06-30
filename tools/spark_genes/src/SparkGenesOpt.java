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
	/* String Options*/
	ACC					( 'a', true, "sra-accession to process arg=SRR123456" ),
	GTF					( 'f', true, "name of gtf-file to process arg=file.gtf" ),
	OUT					( 'o', true, "where to write the result" ),
	TRANS				( 'x', true, "translate references ( from GTF into Accession )" ),
	FEATURE_TYPE		( 't', true, "which feature_type to filter (exon)" ),
	FEATURE_ID			( 'i', true, "which feature_id to use (gene_id)" ),
	MASTER				( 'd', true, "where the spark-master is ( in case of run_mode = SPARKED )" ),

	/* Enum Options*/
	COUNT_MODE			( 'c', true, "how to count alignments UNION/STRICT/NONEMPTY" ),
	RUN_MODE			( 'r', true, "how to run local_s/local_p/sparked" ),

	/* Integer Options*/
	MIN_MAPQ			( 'm', true, "minimal mapping quality" ),
	SLICES				( 's', true, "number of slices" ),
	BINSIZE				( 'j', true, "how big the gtf-lookup bins are (50000)"),
	
	/* boolean Options */
	PROGRESS			( 'p', false, "show progress while processing gtf-file" ),
	PRESCAN				( 'b', false, "prescan gtf-file for unsorted references" ),
	IGNORE_ORIENTATION	( 'g', false, "ignore orientation of alignment" ),
	MEASURE_TIME		( 'e', false, "measure time to execute function" ),
	COUNT_UNALIGNED		( 'u', false, "count unaligned reads" ),
	HELP				( 'h', false, "print help" ),

	NONE				( '.', false, "none" );
	
	// private members of each enum-value
	private char priv_short_opt;
	private boolean priv_has_value;
	private String priv_help_line;
	
	// private constructor of each enum-value
	private opt_sel( char short_opt, boolean has_value, final String help_line )
	{
		priv_short_opt = short_opt;
		priv_has_value = has_value;
		priv_help_line = help_line;
	}

	public String help() { return priv_help_line; }
	public boolean has_value() { return priv_has_value; }
	
	public static opt_sel from_string( final String s )
	{
		opt_sel res = NONE;
		try
		{
			res = Enum.valueOf( opt_sel.class, s.trim().toUpperCase() );
		}
        catch ( IllegalArgumentException e )
        {
        }

		if ( res == NONE && s.length() == 1 )
			res = from_char( s.charAt( 0 ) );
		return res;
	}

	public static opt_sel from_char( final char c )
	{
		for ( opt_sel o : opt_sel.values() )
		{
		  if ( o.priv_short_opt == c ) return o;
		}
		return NONE;
	}
	
	public static String help_string( opt_sel[] options )
	{
		StringBuffer sb = new StringBuffer();
		for ( opt_sel o : options )
		{
			if ( o != NONE )
			{
				String s = String.format( "  -%c|--%s", o.priv_short_opt, o.toString().toLowerCase() );
				sb.append( String.format( "%-30s %s%n", s, o.priv_help_line ) );
			}
		}
		return sb.toString();
	}
	
}

public class SparkGenesOpt
{
    private String Accession;
    private String GTFfile;
    private String Outputfile;
    private String RefTranslation;
    private String FeatureType;
    private String FeatureID;
    private String SparkMaster;

    private count_mode CountMode;
    private run_mode RunMode;

    private int MinMapq;
    private int NumSlices;
    private int LookupBinSize;
	
    private boolean ShowProgress;
    private boolean PreScanGtf;
    private boolean IgnoreOrientation;
    private boolean MeasureTime;
    private boolean CountUnaligned;
    private boolean HelpRequested;
	
    public SparkGenesOpt( String[] args )
    {
        Accession = null;
        GTFfile = null;
        Outputfile = null;
        RefTranslation = "";
        FeatureType = "exon";
        FeatureID = "gene_id";
        SparkMaster = null;
		
        CountMode = count_mode.UNION;
        RunMode = run_mode.LOCAL_P;
		
        MinMapq = 0;
        NumSlices = 8;
        LookupBinSize = 500000;
		
        ShowProgress = false;
        PreScanGtf = false;
        IgnoreOrientation = false;
        MeasureTime = false;
        CountUnaligned = false;
        HelpRequested = false;
		
		from_args( args );
    }

	private void set_by_enum_and_value( opt_sel sel, final String value )
	{
		switch( sel )
		{
			case ACC					: Accession = value; break;
			case GTF					: GTFfile = value; break;
			case OUT					: Outputfile = value; break;
			case TRANS					: RefTranslation = value; break;
			case FEATURE_TYPE			: FeatureType = value; break;
			case FEATURE_ID				: FeatureID = value; break;
			case MASTER					: SparkMaster = value; break;

			case COUNT_MODE				: CountMode = count_mode.from_string( value ); break;
			case RUN_MODE				: RunMode = run_mode.from_string( value ); break;

			case MIN_MAPQ				: MinMapq = Integer.parseInt( value ); break;
			case SLICES					: NumSlices = Integer.parseInt( value ); break;
			case BINSIZE				: LookupBinSize = Integer.parseInt( value ); break;
			
			case PROGRESS				: ShowProgress = Boolean.parseBoolean( value ); break;
			case PRESCAN				: PreScanGtf = Boolean.parseBoolean( value ); break;
			case IGNORE_ORIENTATION		: IgnoreOrientation = Boolean.parseBoolean( value ); break;
			case MEASURE_TIME			: MeasureTime = Boolean.parseBoolean( value ); break;
			case COUNT_UNALIGNED		: CountUnaligned = Boolean.parseBoolean( value ); break;
			case HELP					: HelpRequested = Boolean.parseBoolean( value ); break;
		}
	}
	
	private void read_from_file( final String filename )
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
						set_by_enum_and_value( opt_sel.from_string( key ), value );
					}
					else
						set_by_enum_and_value( opt_sel.from_string( line ), "true" );
				}
			} while ( line != null );
		}
        catch ( IOException e )
        {
            System.out.printf( "cannot read file '%s'%n", filename );
        }
	}

    private void from_args( String[] args )
	{
		opt_sel sel = opt_sel.NONE;
		String value = "";
		for ( String arg : args )
		{
			if ( sel == opt_sel.NONE )
			{
				if ( arg.startsWith( "--" ) )
					sel = opt_sel.from_string( arg.substring( 2 ) );
				else if ( arg.startsWith( "-" ) )
				{
					sel = opt_sel.from_string( arg.substring( 1 ) );
					if ( !sel.has_value() )
					{
						set_by_enum_and_value( sel, "true" );
						sel = opt_sel.NONE;
					}
				}
			}
			else
			{
				set_by_enum_and_value( sel, arg );
				sel = opt_sel.NONE;
			}
		}	
	}

    public void print_help()
    {
		System.out.println( opt_sel.help_string( opt_sel.values() ) );
    }

	public String get_Raw_Outputfile() { return Outputfile; }
		
    public String get_Outputfile()
    {
        if ( Outputfile != null )
            return Outputfile;
        else if ( Accession != null )
            return String.format( "%s.counts.txt", Accession );
        else
            return "counts.txt";
    }

    public String get_Accession() { return Accession; }
    public String get_GTFfile() { return GTFfile; }
    public String get_RefTranslation() { return RefTranslation; }
    public boolean has_RefTranslation() { return ( RefTranslation != "" ); }
    public String get_FeatureType() { return FeatureType; }
    public String get_FeatureID() { return FeatureID; }
    public String get_SparkMaster() { return SparkMaster; }

    public count_mode get_CountMode() { return CountMode; }
	public run_mode get_RunMode() { return RunMode;}
	
    public int get_MinMapq() { return MinMapq; }
    public int get_NumSlices() { return NumSlices; }
    public int get_LookupBinSize() { return LookupBinSize; }

    public boolean get_ShowProgress() { return ShowProgress; }
    public boolean get_PreScanGtf() { return PreScanGtf; }
    public boolean get_IgnoreOrientation() { return IgnoreOrientation; }
    public boolean get_MeasureTime() { return MeasureTime; }
    public boolean get_CountUnaligned() { return CountUnaligned; }
    public boolean get_HelpRequested() { return HelpRequested; }
}
