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
package HttpWrap;

public class HttpDnldResolved
{
    String accession;
    String obj_id;
    String name;
    long size;
    String mod_date;
    String md5;
    String dnld_ticket;
    String url;
    int result_code;
    String msg;
    
    public long get_size() { return size; }
    public String get_md5() { return md5; }
    public String get_url() { return url; }
    
    @Override public String toString()
    {
        StringBuilder sb = new StringBuilder();
        sb.append( String.format( "accession   : %s\n", accession ) );
        sb.append( String.format( "obj_id      : %s\n", obj_id ) );
        sb.append( String.format( "name        : %s\n", name ) );
        sb.append( String.format( "size        : %d\n", size ) );
        sb.append( String.format( "mod_date    : %s\n", mod_date ) );
        sb.append( String.format( "md5         : %s\n", md5 ) );
        sb.append( String.format( "dnld-ticket : %s\n", dnld_ticket ) );
        sb.append( String.format( "url         : %s\n", url ) );
        sb.append( String.format( "res-code    : %d\n", result_code ) );
        sb.append( String.format( "message     : %s\n", msg ) );
        return sb.toString();
    }
   
}
