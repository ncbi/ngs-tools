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

class test_case
{
    interval_list alig;
    interval_list gene;

    String s_hdr, s_ali, s_gen;

    void print()
    {
		System.out.println( s_hdr );
		System.out.println( s_ali );
		System.out.println( s_gen );
        if ( alig != null ) System.out.println( "alig: " + alig.toString() );
		if ( gene != null ) System.out.println( "gene: " + gene.toString() );
    }
}

class test_case_01 extends test_case
{
    test_case_01( String prefix )
    {
		s_hdr = prefix + "case 01: alig before gene";
		s_ali = "alig: XXXXXXXXXX";
		s_gen = "gene:              XXXXXXXXXX";
		alig = new interval_list( 10, 50,25 );
		gene = new interval_list( 10, 100,100 );
    }
}

class test_case_02 extends test_case
{
    test_case_02( String prefix )
    {
		s_hdr = prefix + "case 02: alig after gene";
		s_ali = "alig:               XXXXXXXXXX";
		s_gen = "gene:  XXXXXXXXXX";
		alig = new interval_list( 10, 300,100 );
		gene = new interval_list( 10, 100,100 );
    }
}

class test_case_03 extends test_case
{
    test_case_03( String prefix )
    {
		s_hdr = prefix + "case 03: alig overlaps gene at begin";
		s_ali = "alig:  XXXXXXXXXX";
		s_gen = "gene:       XXXXXXXXXX";
		alig = new interval_list( 10, 50,100 );
		gene = new interval_list( 10, 100,100 );
    }
}

class test_case_04 extends test_case
{
    test_case_04( String prefix )
    {
		s_hdr = prefix + "case 04: alig overlaps gene completely";
		s_ali = "alig:  XXXXXXXXXXXXXXXXXXX";
		s_gen = "gene:       XXXXXXXXXX";
		alig = new interval_list( 10, 50,200 );
		gene = new interval_list( 10, 100,100 );
    }
}

class test_case_05 extends test_case
{
    test_case_05( String prefix )
    {
		s_hdr = prefix + "case 05: alig is overlaped by gene completely";
		s_ali = "alig:         XXXX";
		s_gen = "gene:       XXXXXXXXXX";
		alig = new interval_list( 10, 120,20 );
		gene = new interval_list( 10, 100,100 );
    }
}

class test_case_06 extends test_case
{
    test_case_06( String prefix )
    {
		s_hdr = prefix + "case 06: alig overlaps at the end of the gene";
		s_ali = "alig:             XXXXXXXXX";
		s_gen = "gene:       XXXXXXXXXX";
		alig = new interval_list( 10, 150,100 );
		gene = new interval_list( 10, 100,100 );		
	}
}

class test_case_07 extends test_case
{
	test_case_07( String prefix )
	{
		s_hdr = prefix + "case 07: alig overlaps with an amb. gene at start";
		s_ali = "alig:   XXXXXXXX";
		s_gen = "gene:       AAAAAAAAA";
		alig = new interval_list( 10, 50,100 );
		gene = new interval_list( 10 );		
		gene.add( 100, 100, true );
	}
}

class test_case_08 extends test_case
{
	test_case_08( String prefix )
	{
		s_hdr = prefix + "case 08: alig overlaps amb. gene completely";
		s_ali = "alig:   XXXXXXXXXXXXXXXXXX";
		s_gen = "gene:       AAAAAAAAA";
		alig = new interval_list( 10, 50,200 );
		gene = new interval_list( 10 );		
		gene.add( 100, 100, true );
	}
}

class test_case_09 extends test_case
{
	test_case_09( String prefix )
    {
		s_hdr = prefix + "case 09: alig is overlaped by amb. gene completely";
		s_ali = "alig:         XXXX";
		s_gen = "gene:       AAAAAAAAA";
		alig = new interval_list( 10, 150,20 );
		gene = new interval_list( 10 );		
		gene.add( 100, 100, true );
	}
}

class test_case_10 extends test_case
{
	test_case_10( String prefix )
    {
		s_hdr = prefix + "case 10: alig overlaps with amb. gene at end";
		s_ali = "alig:           XXXXXXXXXX";
		s_gen = "gene:       AAAAAAAAA";
		alig = new interval_list( 10, 150,100 );
		gene = new interval_list( 10 );		
		gene.add( 100, 100, true );
	}
}

class test_case_11 extends test_case
{
	test_case_11( String prefix )
    {
		s_hdr = prefix + "case 11: alig overlaps with non-amb. section of amb. gene";
		s_ali = "alig:  XXXXXXXXX";
		s_gen = "gene:       XXXXAAAAAAAAA";
		alig = new interval_list( 10, 50,100 );
		gene = new interval_list( 10 );		
		gene.add( 100, 50, false );
		gene.add( 150, 50, true );
	}
}

class test_case_12 extends test_case
{
	test_case_12( String prefix )
    {
		s_hdr = prefix + "case 12: alig overlaps with non-amb. and of amb. section of gene";
		s_ali = "alig:  XXXXXXXXXXXXXX";
		s_gen = "gene:       XXXXAAAAAAAAA";
		alig = new interval_list( 10, 50,120 );
		gene = new interval_list( 10 );		
		gene.add( 100, 50, false );
		gene.add( 150, 50, true );
	}
}

class test_case_13 extends test_case
{
	test_case_13( String prefix )
    {
		s_hdr = prefix + "case 13: alig overlaps with non-amb. and of amb. section of gene";
		s_ali = "alig:         XXXXXXX";
		s_gen = "gene:       XXXXAAAAAAAAA";
		alig = new interval_list( 10, 120,50 );
		gene = new interval_list( 10 );		
		gene.add( 100, 50, false );
		gene.add( 150, 50, true );
	}
}

class test_case_14 extends test_case
{
	test_case_14( String prefix )
    {
		s_hdr = prefix + "case 14: alig spanning 2 gene-sections";
		s_ali = "alig:      XXXXXXX";
		s_gen = "gene:    XXXX  XXXXXXX";
		alig = new interval_list( 10, 100,100 );
		gene = new interval_list( 10, 90,20, 190,20 );		
	}
}

class test_case_15 extends test_case
{
	test_case_15( String prefix )
    {
		s_hdr = prefix + "case 15: 2 alig segments matching 2 gene-sections";
		s_ali = "alig:      XX  XX";
		s_gen = "gene:    XXXX  XXXXXXX";
		alig = new interval_list( 10, 100,20, 150,20 );
		gene = new interval_list( 10, 90,30, 150,50 );		
	}
}

class test_case_16 extends test_case
{
	test_case_16( String prefix )
    {
		s_hdr = prefix + "case 16: alig matches start of amb gene-section";
		s_ali = "alig:    XXXX";
		s_gen = "gene:    AAAAAAAAA";
		alig = new interval_list( 10, 100,50 );
		gene = new interval_list( 10 );		
        gene.add( 100, 100, true );
	}
}

class test_case_17 extends test_case
{
	test_case_17( String prefix )
    {
		s_hdr = prefix + "case 17: alig matches end of amb gene-section";
		s_ali = "alig:         XXXX";
		s_gen = "gene:    AAAAAAAAA";
		alig = new interval_list( 10, 180,20 );
		gene = new interval_list( 10 );		
        gene.add( 100, 100, true );
	}
}

class test_case_18 extends test_case
{
	test_case_18( String prefix )
    {
		s_hdr = prefix + "case 18: alig matches end of amb gene-section + 1";
		s_ali = "alig:         XXXXX";
		s_gen = "gene:    AAAAAAAAA";
		alig = new interval_list( 10, 180,21 );
		gene = new interval_list( 10 );		
        gene.add( 100, 100, true );
	}
}

class test_case_19 extends test_case
{
	test_case_19( String prefix )
    {
		s_hdr = prefix + "case 19: alig matches the complete amb gene-section";
		s_ali = "alig:    XXXXX";
		s_gen = "gene:    AAAAA";
		alig = new interval_list( 10, 100,50 );
		gene = new interval_list( 10 );	
        gene.add( 100, 50, true );
	}
}

public class interval_list_test
{

	static boolean print_result( boolean b )
	{
		if ( b )
			System.out.println( "\u001B[32mpassed\u001B[0m\n" );
		else
			System.out.println( "\u001B[31mfailed\u001B[0m\n" );
		return b;
	}
	
	static boolean compare_ivs( interval_list l1, interval_list l2 )
	{
		return print_result( l1.equals( l2 ) );
	}

	/* sorting --------------------------------------------------------------------------------- */
	static boolean test_sort()
	{
		System.out.println( "test sorting intervalls" );
		interval_list l = new interval_list( 10, 220,20, 300,10, 210,20 );
		System.out.println( "before: " + l.toString() );
		l.sort();
		System.out.println( "after : " + l.toString() );		

		return compare_ivs( l, new interval_list( 10, 210,20, 220,20, 300,10 ) );
	}
	
	/* merging --------------------------------------------------------------------------------- */
    static boolean merge_common( interval_list l1, interval_list l2 )
    {
		System.out.println( "before: " + l1.toString() );
		l1.merge_intervals();
		System.out.println( "after : " + l1.toString() );		
		return compare_ivs( l1, l2 );
    }

	static boolean merge_1()
	{
		System.out.println( "test merge: 2 none-overlapping intervalls / do not merge" );
		return merge_common( new interval_list( 10, 100,20, 200,20 ),
                              new interval_list( 10, 100,20, 200,20 ) );
	}

	static boolean merge_2()
	{
		System.out.println( "test merge: 2 overlapping intervalls / do merge" );
		return merge_common( new interval_list( 10, 100,100, 150,100 ),
                              new interval_list( 10, 100,150 ) );
	}
	
	static boolean merge_3()
	{
		System.out.println( "test merge: 2 overlapping intervalls / do merge" );
		return merge_common( new interval_list( 10, 100,100, 50,200 ),
                              new interval_list( 10, 50,200 ) );
	}

	static boolean merge_4()
	{
		System.out.println( "test merge: 2 overlapping intervalls / do merge" );
		return merge_common( new interval_list( 10, 100,100, 150,20 ),
                              new interval_list( 10, 100,100 ) );
	}

	static boolean test_merge()
	{
		boolean res = merge_1();
		if ( res ) res = merge_2();
		if ( res ) res = merge_3();
		if ( res ) res = merge_4();
		return res;
	}
	
	/* detect ambiguity --------------------------------------------------------------------------------- */
	static boolean amb_common( interval_list ovr, interval_list src, interval_list expect )
	{
		System.out.println( "ovr: " + ovr.toString() );
		System.out.println( "src: " + src.toString() );
		src.detect_amb( ovr );
		System.out.println( "res: " + src.toString() );
		return print_result( src.equals( expect ) );
	}
	
	static boolean amb_01()
	{
		System.out.println( "case amb_01: ovr befor src" );
		System.out.println( "ovr: XXXXXXXXXX" );
		System.out.println( "src:              XXXXXXXXXX" );
		System.out.println( "res:              XXXXXXXXXX" );

		interval_list ovr = new interval_list( 10, 50,25 );
		interval_list src = new interval_list( 10, 100,100 );		
		interval_list expect = new interval_list( 10, 100,100 );
		
		return amb_common( ovr, src, expect );
	}

	static boolean amb_02()
	{
		System.out.println( "case amb_02: ovr after src" );
		System.out.println( "ovr:                              XXXXXXXXXX" );
		System.out.println( "src:              XXXXXXXXXX" );
		System.out.println( "res:              XXXXXXXXXX" );

		interval_list ovr = new interval_list( 10, 300,100 );
		interval_list src = new interval_list( 10, 100,100 );		
		interval_list expect = new interval_list( 10, 100,100 );
		
		return amb_common( ovr, src, expect );
	}

	static boolean amb_03()
	{
		System.out.println( "case amb_03: complete overlap between ovr and src, src turns ambiguous" );
		System.out.println( "ovr:        XXXXXXXXXXXXXXXXXXXXXX" );
		System.out.println( "src:              XXXXXXXXXX" );
		System.out.println( "res:              AAAAAAAAAA" );

		interval_list ovr = new interval_list( 10, 50,200 );
		interval_list src = new interval_list( 10, 100,100 );		
		interval_list expect = new interval_list( 10 );
		expect.add( 100, 100, true );

		return amb_common( ovr, src, expect );
	}

	static boolean amb_04()
	{
		System.out.println( "case amb_04: ovr overlapps src partially, ovr starts before src" );
		System.out.println( "ovr:         XXXXXXXXXX" );
		System.out.println( "src:              XXXXXXXXXX" );
		System.out.println( "res:              AAAAAXXXXX" );

		interval_list ovr = new interval_list( 10, 50,100 );
		interval_list src = new interval_list( 10, 100,100 );		
		interval_list expect = new interval_list( 10 );
		expect.add( 100, 50, true );
		expect.add( 150, 50, false );
		
		return amb_common( ovr, src, expect );
	}

	static boolean amb_05()
	{
		System.out.println( "case amb_05: ovr overlapps src partially, ovr starts after src" );
		System.out.println( "ovr:                   XXXXXXXXXX" );
		System.out.println( "src:              XXXXXXXXXX" );
		System.out.println( "res:              XXXXXAAAAA" );

		interval_list ovr = new interval_list( 10, 150,100 );
		interval_list src = new interval_list( 10, 100,100 );		
		interval_list expect = new interval_list( 10 );
		expect.add( 100, 50, false );
		expect.add( 150, 50, true );
		
		return amb_common( ovr, src, expect );
	}

	static boolean amb_06()
	{
		System.out.println( "case amb_06: src complete overlapps ovr" );
		System.out.println( "ovr:                 XXX" );
		System.out.println( "src:              XXXXXXXXXX" );
		System.out.println( "res:              XXXAAAXXXX" );

		interval_list ovr = new interval_list( 10, 125,50 );
		interval_list src = new interval_list( 10, 100,100 );		
		interval_list expect = new interval_list( 10 );
		expect.add( 100, 25, false );
		expect.add( 125, 50, true );
		expect.add( 175, 25, false );
		
		return amb_common( ovr, src, expect );
	}
	
	static boolean amb_07()
	{
		System.out.println( "case amb_07: combination of overlaps" );
		System.out.println( "ovr:      XXXXXXXXXXXXXXXXXX" );
		System.out.println( "src:  XXXXXXXXXX       XXXXXXXXXXXX" );
		System.out.println( "res:  XXXXAAAAAA       AAAAAXXXXXXX" );

		interval_list ovr = new interval_list( 10, 150,200 );
		interval_list src = new interval_list( 10, 100,100, 300,100 );		
		interval_list expect = new interval_list( 10 );
		expect.add( 100, 50, false );
		expect.add( 150, 50, true );
		expect.add( 300, 50, true );
		expect.add( 350, 50, false );
		
		return amb_common( ovr, src, expect );
	}

	static boolean amb_08()
	{
		System.out.println( "case amb_08: combination of overlaps" );
		System.out.println( "ovr:  XXXXXXXXXX       XXXXXXXXXXXX" );
		System.out.println( "src:       XXXXXXXXXXXXXXXXX" );
		System.out.println( "res:       AAAAAXXXXXXXAAAAA" );

		interval_list ovr = new interval_list( 10, 50,75, 175,50 );
		interval_list src = new interval_list( 10, 100,100 );		
		interval_list expect = new interval_list( 10 );
		expect.add( 100, 25, true );
		expect.add( 125, 50, false );
		expect.add( 175, 25, true );
		
		return amb_common( ovr, src, expect );
	}

	static boolean amb_09()
	{
		System.out.println( "case amb_09: combination of overlaps" );
		System.out.println( "ovr:  XXXXX  XXX    XXXXXXXXXXX" );
		System.out.println( "src:    XXXXXXXXXXXXXXXXX   XXXXXXX  XXXXX" );
		System.out.println( "res:    AAAXXAAAXXXXAAAAA   AAAXXXX  XXXXX" );

		interval_list ovr = new interval_list( 10, 50,75, 150,10, 190,60 );
		interval_list src = new interval_list( 10, 100,100, 220,70, 300,50 );
		interval_list expect = new interval_list( 10 );
		expect.add( 100, 25, true );
		expect.add( 125, 25, false );
		expect.add( 150, 10, true );
		expect.add( 160, 30, false );
		expect.add( 190, 10, true );
		expect.add( 220, 30, true );
		expect.add( 250, 40, false );
		expect.add( 300, 50, false );
		
		return amb_common( ovr, src, expect );
	}

	static boolean test_amb()
	{
		boolean res = amb_01();
		if ( res ) res = amb_02();
		if ( res ) res = amb_03();
		if ( res ) res = amb_04();
		if ( res ) res = amb_05();
		if ( res ) res = amb_06();
		if ( res ) res = amb_07();
		if ( res ) res = amb_08();
		if ( res ) res = amb_09();
		return res;
	}

	/* event-iter for interval-list ---------------------------------------------------------------- */	
	static boolean print_event_iter( interval_list l )
	{
		long event[] = new long[ 2 ];
		System.out.println( "l: " + l.toString() );
		l.reset_event_iter();
		while ( l.get_event( event ) )
			System.out.print( "[" + event[ 0 ] + ", " + event[ 1 ] + "] " );
		System.out.println( "." );
		return print_result( true );
	}
	
	static boolean event_iter_01()
	{
		System.out.println( "event_iter 01: read a list of intervals as event-stream" );
		interval_list l = new interval_list( 10, 10,10, 30,10 );
		l.add( 50, 10, true );
		return print_event_iter( l );
	}

	static boolean event_iter_02()
	{
		System.out.println( "event_iter 02: read a list of intervals as event-stream" );
		interval_list l = new interval_list( 10, 10,10, 30,20 );
		l.add( 50, 10, true );
		return print_event_iter( l );
	}

	static boolean event_iter_03()
	{
		System.out.println( "event_iter 03: read a list of intervals as event-stream" );
		interval_list l = new interval_list( 10, 10,10, 30,20 );
		l.add( 51, 10, true );
		return print_event_iter( l );
	}

	static boolean test_event_iter()
	{
		boolean res  = event_iter_01();
		if ( res ) res = event_iter_02();
		if ( res ) res = event_iter_03();
		return res;
	}

	/* ------------------------------------------------------------------------------------------- */		

	static boolean cmp_common( test_case tst_case,
                               count_result exp_union, count_result exp_strict, count_result exp_nonempty )
	{
        tst_case.print();
        interval_list_cmp cmp = new interval_list_cmp( true );

        count_result res_union = cmp.walk( tst_case.gene, tst_case.alig, count_mode.UNION );
        System.out.println( "union...: result=" + res_union + " expected=" + exp_union );

        count_result res_strict = cmp.walk( tst_case.gene, tst_case.alig, count_mode.STRICT );
        System.out.println( "strict..: result=" + res_strict + " expected=" + exp_strict );

        count_result res_nonempty = cmp.walk( tst_case.gene, tst_case.alig, count_mode.NONEMPTY );
        System.out.println( "nonempty: result=" + res_nonempty + " expected=" + exp_nonempty );

		return print_result( res_union == exp_union &&
                              res_strict == exp_strict &&
                              res_nonempty == exp_nonempty );
	}

    static boolean test_cmp()
    {
        String prefix = "cmp.";
		boolean res  = cmp_common( new test_case_01( prefix ),
                            count_result.NO_FEATURE, count_result.NO_FEATURE, count_result.NO_FEATURE );

        if ( res ) res = cmp_common( new test_case_02( prefix ),
                            count_result.NO_FEATURE, count_result.NO_FEATURE, count_result.NO_FEATURE );

        if ( res ) res = cmp_common( new test_case_03( prefix ),
                            count_result.FEATURE, count_result.NO_FEATURE, count_result.FEATURE );

        if ( res ) res = cmp_common( new test_case_04( prefix ),
                            count_result.FEATURE, count_result.NO_FEATURE, count_result.FEATURE );

        if ( res ) res = cmp_common( new test_case_05( prefix ),
                            count_result.FEATURE, count_result.FEATURE, count_result.FEATURE );

        if ( res ) res = cmp_common( new test_case_06( prefix ),
                            count_result.FEATURE, count_result.NO_FEATURE, count_result.FEATURE );

        if ( res ) res = cmp_common( new test_case_07( prefix ),
                            count_result.AMBIGUOUS, count_result.NO_FEATURE, count_result.AMBIGUOUS );

        if ( res ) res = cmp_common( new test_case_08( prefix ),
                            count_result.AMBIGUOUS, count_result.NO_FEATURE, count_result.AMBIGUOUS );

        if ( res ) res = cmp_common( new test_case_09( prefix ),
                            count_result.AMBIGUOUS, count_result.AMBIGUOUS, count_result.AMBIGUOUS );

        if ( res ) res = cmp_common( new test_case_10( prefix ),
                            count_result.AMBIGUOUS, count_result.NO_FEATURE, count_result.AMBIGUOUS );
    
        if ( res ) res = cmp_common( new test_case_11( prefix ),
                            count_result.FEATURE, count_result.NO_FEATURE, count_result.FEATURE );

        if ( res ) res = cmp_common( new test_case_12( prefix ),
                            count_result.AMBIGUOUS, count_result.NO_FEATURE, count_result.FEATURE );

        if ( res ) res = cmp_common( new test_case_13( prefix ),
                            count_result.AMBIGUOUS, count_result.FEATURE, count_result.FEATURE );

        if ( res ) res = cmp_common( new test_case_14( prefix ),
                            count_result.FEATURE, count_result.NO_FEATURE, count_result.FEATURE );

        if ( res ) res = cmp_common( new test_case_15( prefix ),
                            count_result.FEATURE, count_result.FEATURE, count_result.FEATURE );

        if ( res ) res = cmp_common( new test_case_16( prefix ),
                            count_result.AMBIGUOUS, count_result.AMBIGUOUS, count_result.AMBIGUOUS );

        if ( res ) res = cmp_common( new test_case_17( prefix ),
                            count_result.AMBIGUOUS, count_result.AMBIGUOUS, count_result.AMBIGUOUS );

        if ( res ) res = cmp_common( new test_case_18( prefix ),
                            count_result.AMBIGUOUS, count_result.NO_FEATURE, count_result.AMBIGUOUS );

        if ( res ) res = cmp_common( new test_case_19( prefix ),
                            count_result.AMBIGUOUS, count_result.AMBIGUOUS, count_result.AMBIGUOUS );

        return res;
    }

	/* ------------------------------------------------------------------------------------------- */

	public static boolean test()
	{
		boolean res = test_sort();
		if ( res ) res = test_merge();
		if ( res ) res = test_amb();

		if ( res ) res = test_event_iter();
        if ( res ) res = test_cmp();

		System.out.println( "all tests:" );
		return print_result( res );
	}

}