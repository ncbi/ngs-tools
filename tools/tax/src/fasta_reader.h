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

#pragma once

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <assert.h>

#include "reader.h"

class FastaReader final: public Reader 
{
private:
    mutable std::istream *f;
    mutable std::ifstream f_of_filename;

    size_t filesize = 0;
    size_t spot_count = 0, read_count = 0;
    std::string last_desc;
    std::string last_spot_id;
    std::string tmp_line;


    static bool is_description(const std::string &s)
	{
		return !s.empty() && s[0] == '>';
	}

    void read_line(std::string& line) {
        std::getline(*f, line);

        // handling windows line endings
        if (!line.empty() && *line.rbegin() == '\r')
            line.erase(line.size() - 1);
    }

    static bool ends_with(const std::string &s, const std::string &end)
    {
        if (end.size() > s.size()) 
            return false;

        return std::equal(end.rbegin(), end.rend(), s.rbegin());
    }
public:

    static bool is_fasta(const std::string &filename) {
        return ends_with(filename, ".fasta") || ends_with(filename, ".fa") || ends_with(filename, ".fna") || filename == "stdin";
    }
    
	FastaReader(const std::string &filename)
        : f_of_filename(filename, std::ios::binary)
	{
        if (filename != "stdin")
        {
            if (f_of_filename.fail())
    		    throw std::runtime_error(std::string("cannot open the file: ") + filename);

            f = &f_of_filename;
            f->seekg(0, std::ios::end);
            filesize = f->tellg();
            f->seekg(0, std::ios::beg);
        }
        else
        {
            f = &std::cin;
            std::ios_base::sync_with_stdio(false);
            std::cin.tie(nullptr);
        }
        
        std::string line;
        read_line(line);

//		if (line.empty())
//			throw std::runtime_error("fasta file is empty");

//		if (!is_description(line))
//			throw std::runtime_error("this is not a fasta file");

		last_desc = line;
	}

    size_t file_size() const { return filesize; }

    virtual SourceStats stats() const override 
    {
        assert(f->eof());
        return SourceStats(spot_count, read_count);
    }
    
    virtual float progress() const override 
    {
        if (filesize == 0)
            return 0;

        return f->eof() ? 1.0f : float(f->tellg()) / filesize;
    }

    bool read(Fragment* output) override 
    {
        if (f->eof())
            return false;
        if (output) {    
            auto end_pos = last_desc.find_first_of(" /");
            if (end_pos == std::string::npos)
                end_pos = last_desc.size();

#if 0
            auto start_pos = end_pos - 1;
            while (start_pos > 0 && isdigit(last_desc[start_pos]))
                start_pos--;
#else
            int start_pos = 0;
#endif

            output->spotid.assign(last_desc, start_pos + 1, end_pos - 1 - start_pos);
            output->bases.clear();
            output->bases.reserve(300); // todo: tune

            while (!f->eof()) 
            {
                read_line(tmp_line);
                if (is_description(tmp_line)) {
                    last_desc = tmp_line;
                    break;
                } 
                output->bases += tmp_line;
            }

            if (output->bases.empty())
                throw std::runtime_error(std::string("Read is empty: ") + last_desc);
//            output->bases.shrink_to_fit(); todo: tune
            read_count ++;
            if (last_spot_id != output->spotid)
                spot_count ++;
            last_spot_id = output->spotid;

        } else {

            while (!f->eof()) 
            {
                read_line(tmp_line);
                if (is_description(tmp_line)) {
                    last_desc = tmp_line;
                    break;
                } 
            }
        }

        return true;
    }
};

