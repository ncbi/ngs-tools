package Bio;

import data.Settings;
import java.io.File;
import job.JobData;
import job.JobFormat;
import job.JobList;
import job.JobState;

public class BioRefSpec
{
    private final String canonical_name;
    private final String common_name;
    private final Boolean circular;
    private final Long length;
    private final Boolean downloaded;
    
    public String name() { return canonical_name; }
    public Boolean get_downloaded() { return downloaded; }
    
    public JobData make_job( JobList jl )
    {
        JobData res = null;
        if ( !jl.contains_job( String.format( "%s.JOB", canonical_name ) ) )
        {
            BioSpec spec = new BioSpec( BioAccessionType.REF_SEQUENCE, 0, canonical_name );
            res = new JobData( spec );
            res.set_state( JobState.READY, false );
            res.set_format( JobFormat.DOWNLOAD );
            res.set_downloadpath( get_dnload_path() );
            res.set_valid( true );
            if ( !res.save( Settings.getInstance().get_jobpath() ) )
                res = null;
        }
        return res;
    }
    
    public String desc()
    {
        if ( circular )
            return String.format( "\'%s\' : %d bases, circular",
                    common_name, length );
        else
            return String.format( "\'%s\' : %d bases",
                    common_name, length );
    }
    
    private String get_dnload_path()
    {
        String sep = File.separator;
        String home = System.getProperty( "user.home" );
        return String.format( "%s%sncbi%spublic%srefseq",
                home, sep, sep, sep );
    }
    
    private String get_local_path( final String c_name )
    {
        String sep = File.separator;
        return String.format( "%s%s%s",
                get_dnload_path(), sep, c_name );
    }
    
    private Boolean get_already_downloaded( final String path )
    {
        File f = new File( path );
        return ( f.exists() && !f.isDirectory() );
    }
    
    public BioRefSpec( final String can, final String com,
            final Boolean circ, final Long len )
    {
        canonical_name = can;
        common_name = com;
        circular = circ;
        length = len;
        downloaded = get_already_downloaded( get_local_path( can ) );
    }
}
