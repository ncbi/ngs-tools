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

enum State
{
	EMPTY   { public String toString() { return "EMPTY"; } },
	GEN     { public String toString() { return "GEN"; } },
	AMB     { public String toString() { return "AMB"; } },
	ALI     { public String toString() { return "ALI"; } },
	GEN_ALI { public String toString() { return "GEN+ALI"; } },
	AMB_ALI { public String toString() { return "AMB+ALI"; } },
}

enum Event
{
    NONE,
    START_ALI,
    END_ALI,

    START_GEN,
    END_GEN,
    START_GEN_START_ALI,
    START_GEN_END_ALI,
    END_GEN_START_ALI,
    END_GEN_END_ALI,

    START_AMB,
    END_AMB,
    START_AMB_START_ALI,
    START_AMB_END_ALI,
    END_AMB_START_ALI,
    END_AMB_END_ALI
}

public class interval_list_cmp
{
	private interval_list gene;
	private interval_list alig;
	private long gene_event[];
	private long alig_event[];
    private State state;
    private int gen_ali, amb_ali, ali;
    private boolean debug;

	public interval_list_cmp( final boolean debug )
	{
        this.debug = debug;
		this.gene = null;
		this.alig = null;
		gene_event = new long[ 2 ];
		alig_event = new long[ 2 ];
	}

	private Event get_alig_event()
	{
		Event res = Event.NONE;
		if ( alig_event[ 1 ] != 0 )
		{
			// but we have events from the alignment
			switch ( (int)alig_event[ 1 ] )
			{
				case 1 : res = Event.START_ALI; break; // start alignment
				case 3 : res = Event.END_ALI; break; // end alignment
			}
			// get the next event from the alignment
			if ( !alig.get_event( alig_event ) ) alig_event[ 1 ] = 0;
		}
		return res;
	}

	private Event get_gene_event()
	{
		Event res = Event.NONE;
		if ( gene_event[ 1 ] != 0 )
		{
			// but we have events from the alignment
			switch( (int)gene_event[ 1 ] )
			{
				case 1 : res = Event.START_GEN; break;  // start gene ( not amb. )
				case 2 : res = Event.START_AMB; break;  // start gene ( amb. )
				case 3 : res = Event.END_GEN; break;    // end gene ( not amb. )
				case 4 : res = Event.END_AMB; break;    // end gene ( amb. )
			}
			// get the next event from the gene
			if ( !gene.get_event( gene_event ) ) gene_event[ 1 ] = 0;		
		}
		return res;
	}

    private Event combine( final Event ae, final Event ge )
    {
        switch( ae )
        {
            case NONE       : return ge;

            case START_ALI  : switch( ge )
                                {
                                    case NONE       : return ae;
                                    case START_GEN  : return Event.START_GEN_START_ALI;
                                    case START_AMB  : return Event.START_AMB_START_ALI;
                                    case END_GEN    : return Event.END_GEN_START_ALI;
                                    case END_AMB    : return Event.END_AMB_START_ALI;
                                }
                            break;

            case END_ALI    : switch( ge )
                                {
                                    case NONE       : return ae;
                                    case START_GEN  : return Event.START_GEN_END_ALI;
                                    case START_AMB  : return Event.START_AMB_END_ALI;
                                    case END_GEN    : return Event.END_GEN_END_ALI;
                                    case END_AMB    : return Event.END_AMB_END_ALI;
                                }
                            break;
        }
        return Event.NONE;
    }

    public void reset( interval_list the_gene, interval_list the_alig )
    {
		gene = the_gene;
        if ( gene != null )
        {
            gene.reset_event_iter();
            if ( !gene.get_event( gene_event ) ) gene_event[ 1 ] = 0;
        }

		alig = the_alig;
        if ( alig != null )
        {
            alig.reset_event_iter();
            if ( !alig.get_event( alig_event ) ) alig_event[ 1 ] = 0;
        }

        state = State.EMPTY;
        gen_ali = amb_ali = ali = 0;
    }

	// return	7 6 5 4 3 2 1 0
	//                        1.... (  1 ) end alignment
	//                      1...... (  2 ) start alignment
	//                    1........ (  4 ) end gene ( not amb. )
	//                  1.......... (  8 ) start gene ( not amb. )
	//                1............ ( 16 ) end gene ( amb. )
	//              1.............. ( 32 ) start gene ( amb. )
	private Event next_event()
	{
		if ( gene_event[ 1 ] == 0 )
			return get_alig_event(); // we have no event from the gene, but eventually from the alignment
		if ( alig_event[ 1 ] == 0 )
			return get_gene_event(); // we have a gene-event but no event from the alignment

		// we have envents from both:
		if ( gene_event[ 0 ] < alig_event[ 0 ] )
			return get_gene_event(); // the gene event comes first...
		if ( alig_event[ 0 ] < gene_event[ 0 ] )
			return get_alig_event();		// the alignment event comes first...

		return combine( get_alig_event(), get_gene_event() ); // from both at the same time!
	}

    private count_result get_result( count_mode mode )
    {
        switch( mode )
        {
            case SIMPLE:   if ( ( gen_ali + amb_ali ) > 0 )
                                return count_result.FEATURE;
                            else
                                return count_result.NO_FEATURE;

            case UNION:    if ( amb_ali > 0 )
                                return count_result.AMBIGUOUS;
                            else if ( gen_ali > 0 )
                                return count_result.FEATURE;
                            else
                                return count_result.NO_FEATURE;

            case STRICT:   if ( ali > 0 )
                                return count_result.NO_FEATURE;
                            else if ( gen_ali > 0 )
                                return count_result.FEATURE;
                            else if ( amb_ali > 0 )
                                return count_result.AMBIGUOUS;
                            else
                                return count_result.NO_FEATURE;

            case NONEMPTY: if ( gen_ali > 0 )
                                return count_result.FEATURE;
                            else if ( amb_ali > 0 )
                                return count_result.AMBIGUOUS;
                            else
                                return count_result.NO_FEATURE;

        }
        return count_result.NO_FEATURE;
    }


    private void enter( final State new_state )
    {
        state = new_state;
        switch( state )
        {
            case ALI        : ali++; break;
            case GEN_ALI    : gen_ali++; break;
            case AMB_ALI    : amb_ali++; break;
        }
    }


    // returns: 0...no_feature, 1...feature, 2..ambiguous
	public count_result walk( interval_list the_gene, interval_list the_alig, count_mode mode )
	{
        reset( the_gene, the_alig );
		Event ev = next_event();
		while ( ev != Event.NONE )
		{
			switch( state )
			{
				case EMPTY : switch( ev )
							{
                                case END_ALI                : break;
                                case START_ALI              : enter( State.ALI ); break;
                                case END_GEN                : break;
                                case END_GEN_END_ALI        : break;
                                case END_GEN_START_ALI      : enter( State.ALI ); break;
                                case START_GEN              : enter( State.GEN ); break;
                                case START_GEN_END_ALI      : enter( State.GEN ); break;
                                case START_GEN_START_ALI    : enter( State.GEN_ALI ); break;
                                case END_AMB                : break;
                                case END_AMB_END_ALI        : break;
                                case END_AMB_START_ALI      : enter( State.ALI ); break;
                                case START_AMB              : enter( State.AMB ); break;
                                case START_AMB_END_ALI      : enter( State.AMB ); break;
                                case START_AMB_START_ALI    : enter( State.AMB_ALI ); break;
							}
							break;

				case GEN :	switch( ev )
							{
                                case END_ALI                : break;
                                case START_ALI              : enter( State.GEN_ALI ); break;
                                case END_GEN                : enter( State.EMPTY ); break;
                                case END_GEN_END_ALI        : enter( State.EMPTY ); break;
                                case END_GEN_START_ALI      : enter( State.ALI ); break;
                                case START_GEN              : break;
                                case START_GEN_END_ALI      : break;
                                case START_GEN_START_ALI    : enter( State.GEN_ALI ); break;
                                case END_AMB                : break;
                                case END_AMB_END_ALI        : break;
                                case END_AMB_START_ALI      : enter( State.GEN_ALI ); break;
                                case START_AMB              : enter( State.AMB ); break;
                                case START_AMB_END_ALI      : enter( State.AMB ); break;
                                case START_AMB_START_ALI    : enter( State.AMB_ALI ); break;
							}
							break;

				case AMB :	switch( ev )
							{
                                case END_ALI                : break;
                                case START_ALI              : enter( State.AMB_ALI ); break;
                                case END_GEN                : break;
                                case END_GEN_END_ALI        : break;
                                case END_GEN_START_ALI      : enter( State.AMB_ALI ); break;
                                case START_GEN              : enter( State.GEN ); break;
                                case START_GEN_END_ALI      : enter( State.GEN ); break;
                                case START_GEN_START_ALI    : enter( State.GEN_ALI ); break;
                                case END_AMB                : enter( State.EMPTY ); break;
                                case END_AMB_END_ALI        : break;
                                case END_AMB_START_ALI      : enter( State.ALI ); break;
                                case START_AMB              : break;
                                case START_AMB_END_ALI      : break;
                                case START_AMB_START_ALI    : enter( State.AMB_ALI ); break;
							}
							break;

				case ALI : switch( ev )
							{
                                case END_ALI                : enter( State.EMPTY ); break;
                                case START_ALI              : break;
                                case END_GEN                : break;
                                case END_GEN_END_ALI        : enter( State.EMPTY ); break;
                                case END_GEN_START_ALI      : break;
                                case START_GEN              : enter( State.GEN_ALI ); break;
                                case START_GEN_END_ALI      : enter( State.GEN ); break;
                                case START_GEN_START_ALI    : enter( State.GEN_ALI ); break;
                                case END_AMB                : break;
                                case END_AMB_END_ALI        : enter( State.EMPTY ); break;
                                case END_AMB_START_ALI      : break;
                                case START_AMB              : enter( State.AMB_ALI ); break;
                                case START_AMB_END_ALI      : enter( State.AMB ); break;
                                case START_AMB_START_ALI    : enter( State.AMB_ALI ); break;
							}
							break;

				case GEN_ALI :	switch( ev )
							{
                                case END_ALI                : enter( State.GEN ); break;
                                case START_ALI              : break;
                                case END_GEN                : enter( State.ALI ); break;
                                case END_GEN_END_ALI        : enter( State.EMPTY ); break;
                                case END_GEN_START_ALI      : enter( State.ALI ); break;
                                case START_GEN              : break;
                                case START_GEN_END_ALI      : enter( State.GEN ); break;
                                case START_GEN_START_ALI    : break;
                                case END_AMB                : break;
                                case END_AMB_END_ALI        : enter( State.GEN ); break;
                                case END_AMB_START_ALI      : break;
                                case START_AMB              : enter( State.AMB_ALI ); break;
                                case START_AMB_END_ALI      : enter( State.AMB ); break;
                                case START_AMB_START_ALI    : enter( State.AMB_ALI ); break;
							}
							break;

				case AMB_ALI :	switch( ev )
							{
                                case END_ALI                : enter( State.AMB ); break;
                                case START_ALI              : break;
                                case END_GEN                : break;
                                case END_GEN_END_ALI        : enter( State.AMB ); break;
                                case END_GEN_START_ALI      : break;
                                case START_GEN              : enter( State.GEN_ALI ); break;
                                case START_GEN_END_ALI      : enter( State.GEN ); break;
                                case START_GEN_START_ALI    : enter( State.GEN_ALI ); break;
                                case END_AMB                : enter( State.ALI ); break;
                                case END_AMB_END_ALI        : enter( State.EMPTY ); break;
                                case END_AMB_START_ALI      : enter( State.ALI ); break;
                                case START_AMB              : break;
                                case START_AMB_END_ALI      : enter( State.AMB ); break;
                                case START_AMB_START_ALI    : break;
							}
							break;

			}
            if ( debug ) System.out.print( "-" + state );
			ev = next_event();
		}
        if ( debug ) System.out.println( " ali=" + ali + " gen+ali=" + gen_ali + " amb+ali=" + amb_ali );

        return get_result( mode );
	}
}
