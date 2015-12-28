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
package Bio;

/*
    to combine:
        bio_type ( Invalid / Read-collection / Ref-sequence )
        count ( row-count )
*/
public class BioSpec
{
    private BioAccessionType bio_type;
    private long count;
    private String accession;
            
    public BioAccessionType get_type() { return bio_type; }
    public void set_type( final BioAccessionType value ) { bio_type = value; }
    
    public long get_count() { return count; }
    public void set_count( final long value ) { count = value; }
    
    public String get_accession() { return accession; }
    public void set_accession( final String value ) { accession = value; }
    
    public final void clear()
    {
        this.bio_type = BioAccessionType.INVALID;
        this.count = 0;
        accession = "";
    }
    
    public boolean is_valid()
    {
        boolean res = !( bio_type.equals( BioAccessionType.INVALID ) );
        if ( res ) res = ( count > 0 );
        return res;
    }

    public void copy( final BioSpec other )
    {
        this.bio_type = other.bio_type;
        this.count = other.count;
        this.accession = other.accession;
   }
    
    public BioSpec() { clear(); }
    
    public BioSpec( final BioAccessionType bio_type, final long count, final String accession )
    {
        this.bio_type = bio_type;
        this.count = count;
        this.accession = accession;
    }
}
