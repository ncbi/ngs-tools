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

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;

public class md5_extractor
{
    public static String md5_of_file( String filename )
    {
        String res = "";
        try
        {
            MessageDigest md5 = MessageDigest.getInstance( "MD5" );
            try
            {
                FileInputStream fs = new FileInputStream( filename );
                byte[] buffer = new byte[ 4096 ];
                int nread;
                while ( ( nread = fs.read( buffer ) ) != -1 )
                    md5.update( buffer, 0, nread );
                fs.close();
                byte[] md5bytes = md5.digest();
                StringBuilder sb = new StringBuilder();
                for ( int i = 0; i < md5bytes.length; i++ )
                    sb.append( Integer.toString( ( md5bytes[ i ] & 0xff ) + 0x100, 16 ).substring( 1 ) );
                res = sb.toString();
            }
            catch ( FileNotFoundException ex )
            {
                CLogger.logfmt( "cannot open '%s' : %s" , filename, ex.toString() );
            }
            catch ( IOException ex )
            {
                CLogger.logfmt( "cannot read from '%s' : %s" , filename, ex.toString() );
            }
        }
        catch ( NoSuchAlgorithmException ex )
        {
            CLogger.logfmt( "error in md5_of_file" , ex.toString() );
        }
        return res;
    }
}
