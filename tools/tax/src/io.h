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
#include <string>
#include "missing_cpp_features.h"

struct IO
{
#if 1
    struct Writer
    {
        std::vector<std::ofstream> out_f;
        std::ostream *stream_f = nullptr;

        Writer(const std::string &_filenames) 
        {
            if (_filenames.empty())
            {
                stream_f = &std::cout;
                return;
            }

            auto filenames = split(_filenames, ',');
            for (auto &filename : filenames)
                out_f.emplace_back(std::ofstream(filename));
            
            if (out_f.size() == 1)
                stream_f = get_out_f(0);
//            f.exceptions( ~std::fstream::goodbit); // todo: think about enabling exceptions here
            check();
        }

        Writer(Writer &writer, int stream_id) : stream_f(writer.get_out_f(stream_id))
        {
        }

        std::ostream &f()
        {
            return *stream_f;
        }

        std::ostream *get_out_f(int stream_id)
        {
            if (stream_id < 0 || stream_id >= out_f.size())
                throw std::runtime_error("Writer:: stream_id < 0 || stream_id >= out_f.size()");

            return &out_f[stream_id];
        }

        void check()
        {
            for (int i = 0; i < out_f.size(); i++)
                if (!out_f[i].good())
                    throw std::runtime_error("failed to write results (no space left on drive?)");

            if (!std::cout.good())
                throw std::runtime_error("failed to write results to std::cout (no space left on drive?)");
        }

        ~Writer()
        {
            for (int i = 0; i < out_f.size(); i++)
                out_f[i].close();

            check();
        }
    };

#else
    struct Writer
    {
        const std::string filename;
        std::ofstream out_f;
        std::ofstream &stream_f;
        int stream_id = -1;

        Writer(const std::string &filename) : filename(filename), out_f(filename), stream_f(out_f)
        {
//            f.exceptions( ~std::fstream::goodbit); // todo: think about enabling exceptions here
            check();
        }

        Writer(const Writer &writer, int stream_id) : filename(writer.filename), stream_f(writer.stream_f), stream_id(stream_id)
        {
        }

        std::ostream &f()
        {
            return filename.empty() ? std::cout : stream_f;
        }

        void check()
        {
            if (!f().good())
                throw std::runtime_error("failed to write results (no space left on drive?)");
        }

        ~Writer()
        {
            if (!filename.empty() && stream_id < 0)
                out_f.close();

            check();
        }
    };
#endif

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
	    f.seekg(offset);
	    if (!f)
		    throw std::runtime_error("load_vector_no_size::cannot seek to offset");

        load_vector_data(f, v, size);
    }

    static size_t filesize(const std::string &filename) 
    {
	    std::ifstream f(filename, std::ios_base::binary);
	    f.seekg(0, std::ios_base::end);
	    return f.tellg();
    }

    static bool file_exists(const std::string &filename)
    {
        std::ifstream f(filename);
        return f.good();
    }

};

#endif