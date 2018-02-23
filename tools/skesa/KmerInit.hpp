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

#ifndef _KmerInit_
#define _KmerInit_

#include "LargeInt.hpp"
#include <boost/variant.hpp>

/******************

This is the only place where we manipulate boost::variant directly. The rest of the code MUST use these definitions and boost::visitor

*******************/

using namespace std;
namespace DeBruijn {

#define MaxPrec 16  // kmer up to 512

    // for TKmer
    typedef boost::variant<LargeInt<1>,LargeInt<2>,LargeInt<3>,LargeInt<4>,LargeInt<5>,LargeInt<6>,LargeInt<7>,LargeInt<8>,
                           LargeInt<9>,LargeInt<10>,LargeInt<11>,LargeInt<12>,LargeInt<13>,LargeInt<14>,LargeInt<15>,LargeInt<16>> TLargeIntN;    

    // for TKmerCount
    template<int N> using TLargeIntVec = vector<pair<LargeInt<N>,size_t>>;
    typedef boost::variant<TLargeIntVec<1>,TLargeIntVec<2>,TLargeIntVec<3>,TLargeIntVec<4>,TLargeIntVec<5>,TLargeIntVec<6>,TLargeIntVec<7>,TLargeIntVec<8>,
                           TLargeIntVec<9>,TLargeIntVec<10>,TLargeIntVec<11>,TLargeIntVec<12>,TLargeIntVec<13>,TLargeIntVec<14>,TLargeIntVec<15>,TLargeIntVec<16>> TKmerCountN;
        
    // for TKmerMap
    struct SKmerHash {
        template<typename T> 
        size_t operator() (const T& kmer) const { return kmer.oahash(); }
    };
    template<int N, class V> using TLargeIntMap = unordered_map<LargeInt<N>,V,SKmerHash>;
    template<class V> using TKmerMapN =  boost::variant<TLargeIntMap<1,V>,TLargeIntMap<2,V>,TLargeIntMap<3,V>,TLargeIntMap<4,V>,TLargeIntMap<5,V>,TLargeIntMap<6,V>,TLargeIntMap<7,V>,TLargeIntMap<8,V>,
                                                        TLargeIntMap<9,V>,TLargeIntMap<10,V>,TLargeIntMap<11,V>,TLargeIntMap<12,V>,TLargeIntMap<13,V>,TLargeIntMap<14,V>,TLargeIntMap<15,V>,TLargeIntMap<16,V>>;

    // for CKmerHashCount
    template<int N, class V>
    struct SOneWayList {
        typedef LargeInt<N> large_t;
        typedef pair<large_t,V> element_t;
        typedef forward_list<element_t> list_t;

        element_t m_data;
        list_t m_extra;
    };
    template<int N, class V> using TOneWayListVec = vector<SOneWayList<N,V>>;
    template<class V> using TKmerHashTable = boost::variant<TOneWayListVec<1,V>,TOneWayListVec<2,V>,TOneWayListVec<3,V>,TOneWayListVec<4,V>,TOneWayListVec<5,V>,TOneWayListVec<6,V>,TOneWayListVec<7,V>,TOneWayListVec<8,V>,
                                                            TOneWayListVec<9,V>,TOneWayListVec<10,V>,TOneWayListVec<11,V>,TOneWayListVec<12,V>,TOneWayListVec<13,V>,TOneWayListVec<14,V>,TOneWayListVec<15,V>,TOneWayListVec<16,V>>;

    
    // This variadic template could be used in construsctors of all boost::variants used in this code
    template<typename Variant, template<int, typename...> class BoundedType, typename... Params>
    Variant CreateVariant(int p) {
        switch(p) {
        case 1 :  return BoundedType<1, Params...>();
        case 2 :  return BoundedType<2, Params...>();
        case 3 :  return BoundedType<3, Params...>();
        case 4 :  return BoundedType<4, Params...>();
        case 5 :  return BoundedType<5, Params...>();
        case 6 :  return BoundedType<6, Params...>();
        case 7 :  return BoundedType<7, Params...>();
        case 8 :  return BoundedType<8, Params...>();
        case 9 :  return BoundedType<9, Params...>();
        case 10 : return BoundedType<10, Params...>();
        case 11 : return BoundedType<11, Params...>();
        case 12 : return BoundedType<12, Params...>();
        case 13 : return BoundedType<13, Params...>();
        case 14 : return BoundedType<14, Params...>();
        case 15 : return BoundedType<15, Params...>();
        case 16 : return BoundedType<16, Params...>();
        default :  throw runtime_error("Not supported kmer length");
        }
    }
            
}; // namespace
#endif /* _KmerInit_ */
