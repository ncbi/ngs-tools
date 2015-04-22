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
import java.util.HashSet;

public class gtf_prescan
{
    private feature_asm ft_asm;
    private String chromosome;
    private int features;
    private long last_start;
    private boolean sorted;
    private boolean show_progress;
    private HashSet<String> sorted_chromosomes;
    private gtf_feature feature;

	gtf_prescan ( feature_asm assembler, final boolean show_progress )
	{
        this.show_progress = show_progress;
        ft_asm = assembler;
        sorted_chromosomes = new HashSet<String>( 32 );
        chromosome = "";
        features = 0;
        last_start = 0;
	}

    void report_chromosome()
    {
        if ( !chromosome.equals( "" ) )
        {
            if ( sorted )
            {
                sorted_chromosomes.add( chromosome );
                if ( show_progress )
                    System.out.println( "\n" + chromosome + " : " + features + " features, sorted" );
            }
            else
            {
                if ( show_progress )
                    System.out.println(  "\n" + chromosome + " : " + features + " features, not sorted" );
            }
        }
    }

    boolean is_chromosome_sorted( final String chr )
    {
        return sorted_chromosomes.contains( chr );
    }

    void run()
    {
        if ( ft_asm != null )
        {
            while ( ( feature = ft_asm.next() ) != null )
            {
                if ( !feature.has_ref( chromosome ) )
                {
                    report_chromosome();
                    chromosome = feature.get_ref();
                    features = 0;
                    last_start = 0;
                    sorted = true;
                }
                features++;
                if ( feature.get_start() < last_start )
                    sorted = false;
                last_start = feature.get_start();
            }
            report_chromosome();
        }
    }
}
