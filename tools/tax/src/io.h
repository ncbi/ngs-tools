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

#ifndef IO_H_INCLUDED
#define IO_H_INCLUDED

#include <fstream>
#include <vector>
#include <stdexcept>

struct IO
{
    template <class C>
    static void save_vector_data(std::ofstream &f, const std::vector<C> &v)
    {
        f.write((char*)&v[0], sizeof(C) * v.size());
    }

    template <class C>
    static void load_vector_data(std::ifstream &f, std::vector<C> &v, size_t size)
    {
        v.clear();
        v.resize(size);
        f.read((char*)&v[0], sizeof(C) * size);

	    if (!f)
		    throw std::runtime_error("load_vector_data:: failed to load vector");
    }

    template <class C>
    static void save_vector(std::ofstream &f, const std::vector<C> &v)
    {
        size_t size = v.size();
        write(f, size);
        save_vector_data(f, v);

	    if (!f)
		    throw std::runtime_error("save_vector:: failed to save vector");
    }

    template <class C>
    static void load_vector(std::ifstream &f, std::vector<C> &v)
    {
        size_t size = 0;
        read(f, size);
        load_vector_data(f, v, size);
    }

    template <class X>
    static void write(std::ofstream &f, const X &x)
    {
		f.write((char*)&x, sizeof(x));
	    if (!f)
		    throw std::runtime_error("IO::write failed");
    }

    template <class X>
    static void read(std::ifstream &f, X &x)
    {
		f.read((char*)&x, sizeof(x));
	    if (!f)
		    throw std::runtime_error("IO::read failed");
    }

    template <class C>
    static void load_vector_no_size(std::ifstream &f, std::vector<C> &v, size_t offset, size_t size)
    {
//	    v.clear();
//	    v.resize(size);
	    f.seekg(offset);
	    if (!f)
		    throw std::runtime_error("load_vector_no_size::cannot seek to offset");

//	    f.read((char*)&v[0], v.size() * sizeof(v[0]));
        load_vector_data(f, v, size);
    }

};

#endif