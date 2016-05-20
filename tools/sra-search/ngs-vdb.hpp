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

#ifndef _hpp_ngs_vdb_
#define _hpp_ngs_vdb_

#include <ngs/ReadCollection.hpp>

#include <string>

struct NGS_VDB_ReadCollection;
struct VBlob;

namespace ngs
{
    namespace vdb
    {
        typedef struct VBlob * BlobRef;

        class FragmentBlob
        {
        public:
            const char* Data() const
                throw ();

            uint64_t Size() const
                throw ();

            void GetFragmentInfo ( uint64_t offset, std::string& p_fragId, uint64_t& p_nextFragStart, bool& biological ) const
                throw ( ErrorMsg );

        public:

            // C++ support

            FragmentBlob ( BlobRef ref )
                throw ();

            FragmentBlob & operator = ( const FragmentBlob & obj )
                throw ( ErrorMsg );
            FragmentBlob ( const FragmentBlob & obj )
                throw ( ErrorMsg );

            ~ FragmentBlob ()
                throw ();

        private:

            FragmentBlob & operator = ( BlobRef ref )
                throw ();

        protected:

            BlobRef m_self;
            NGS_VDB_ReadCollection* m_coll; // not owned here

        protected: // temporary for prototyping
            FragmentBlob (struct NGS_VDB_ReadCollection* );
        };

        /*----------------------------------------------------------------------
        * FragmentBlobIterator
        *  iterates across a list of blobs of read fFragments
        */
        class FragmentBlobIterator : public FragmentBlob
        {
        public:

            /* nextBlob
            *  advance to first Blob on initial invocation
            *  advance to next Blob subsequently
            *  returns false if no more Blob are available.
            *  throws exception if more Blob should be available,
            *  but could not be accessed.
            */
            bool nextBlob ()
                throw ( ErrorMsg );

        public:

            // C++ support

            FragmentBlobIterator ( BlobRef ref )
                throw ();

            FragmentBlobIterator & operator = ( const FragmentBlobIterator & obj )
                throw ( ErrorMsg );
            FragmentBlobIterator ( const FragmentBlobIterator & obj )
                throw ( ErrorMsg );

            ~ FragmentBlobIterator ()
                throw ();

        public: // temporary
            FragmentBlobIterator ( struct NGS_VDB_ReadCollection* );

        private:

            FragmentBlob & operator = ( const FragmentBlob & obj ) // copied from ReadIterator; why do we need this here?
                throw ( ErrorMsg );

            FragmentBlobIterator & operator = ( BlobRef ref )
                throw ();
        };

        class VdbReadCollection //: public ngs :: ReadCollection
        {
        public:
            FragmentBlobIterator getFragmentBlobs () const
                throw ( ErrorMsg );

        public:

            // C++ support

            VdbReadCollection & operator = ( ReadCollectionRef ref )
                throw ();
            VdbReadCollection ( ReadCollectionRef ref )
                throw ();

            VdbReadCollection & operator = ( const VdbReadCollection & obj )
                throw ();
            VdbReadCollection ( const VdbReadCollection & obj )
                throw ();

            ~ VdbReadCollection ()
                throw ();

        // temporary for prototyping
        public:
            VdbReadCollection ( const String & spec )
                throw ();

        private:
            struct NGS_VDB_ReadCollection * m_coll;
        };
    }
} // namespace ngs

///////////////////////////////////// Engine

#include <ngs/ncbi/NGS.hpp>

namespace ncbi
{
    class NGS_VDB : public NGS
    {
    public:

        /* openVdbReadCollection
        *  create an object representing a named collection of reads, with support for blob-based search
        *  "spec" may be a path to an object
        *  or may be an id, accession, or URL
        */
        static
        ngs :: vdb :: VdbReadCollection openVdbReadCollection ( const String & spec )
            throw ( ErrorMsg );
    };
}

#endif // _hpp_ngs_vdb_
