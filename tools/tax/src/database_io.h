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

#ifndef DATABASE_IO_H_INCLUDED
#define DATABASE_IO_H_INCLUDED

#include <string>
#include <vector>
#include <fstream>
#include <stdexcept>

// todo: throw runtime_error

template <class C>
void load_vector(const std::string &filename, std::vector<C> &hash_array, size_t offset = 0)
{
	std::ifstream f(filename, std::ios::binary | std::ios::in);
	if (f.fail() || f.eof())
		throw std::runtime_error(std::string("cannot open ") + filename);

	f.seekg(offset);

	size_t size = 0;
	f.read((char*)&size, sizeof(size));
	if (!size)
		throw std::runtime_error(std::string("cannot read db ") + filename);

	hash_array.resize(size);
	f.read((char*)&hash_array[0], hash_array.size() * sizeof(hash_array[0]));

	if (!f)
		throw std::string("cannot read db completely ") + filename;
}

template <class C>
void load_structure(const std::string &filename, C &c)
{
	std::ifstream f(filename, std::ios::binary | std::ios::in);
	if (f.fail() || f.eof())
		throw std::runtime_error(std::string("cannot open ") + filename);

	f.read((char*)&c, sizeof(c));
}

static size_t filesize(const std::string &filename) 
{
	std::ifstream f(filename, std::ios_base::binary);
	f.seekg(0, std::ios_base::end);
	return f.tellg();
}

template <class C>
static void load_vector_no_size(std::ifstream &f, std::vector<C> &v, size_t offset, size_t size)
{
	v.clear();
	v.resize(size);
	f.seekg(offset);
	if (!f)
		throw std::runtime_error("load_vector_no_size::cannot seek to offset");

	f.read((char*)&v[0], v.size() * sizeof(v[0]));

	if (!f)
		throw std::runtime_error("load_vector_no_size::cannot load vector completely");
}

template <class C>
static size_t size_of_stored_vector(size_t vector_size)
{
	return sizeof(size_t) + vector_size * sizeof(C);
}

#endif
