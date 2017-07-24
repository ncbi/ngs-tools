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
#include <vector>
#include <assert.h>

#include "reader.h"

class FastaReader final: public Reader {
private:
    mutable std::ifstream f;
    size_t fsize;
    size_t spot_idx;
    std::string last_desc;

    static bool is_description(const std::string &s)
	{
		return !s.empty() && s[0] == '>';
	}

    void read_line(std::string& line) {
        std::getline(f, line);
        // handling windows line endings
        if (!line.empty() && *line.rbegin() == '\r') {
            line.erase(line.size() - 1);
        }
    }

    static bool ends_with(const std::string &s, const std::string &end)
    {
        if (end.size() > s.size()) 
            return false;

        return std::equal(end.rbegin(), end.rend(), s.rbegin());
    }
public:

    static bool is_fasta(const std::string &filename) {
        return ends_with(filename, ".fasta") || ends_with(filename, ".fa") || ends_with(filename, ".fna");
    }
    
	FastaReader(const std::string &filename)
        : f(filename, std::ios::binary)
        , spot_idx(0)
	{
        f.seekg(0, std::ios::end);
        fsize = f.tellg();
        f.seekg(0, std::ios::beg);
        
        std::string line;
        read_line(line);
        
		if (line.empty())
			throw std::runtime_error("fasta file is empty");

		if (!is_description(line))
			throw std::runtime_error("this is not a fasta file");

		last_desc = line;
	}

    size_t file_size() const { return fsize; }

    virtual SourceStats stats() const override {
        assert(f.eof());
        return SourceStats(spot_idx);
    }
    
    virtual float progress() const override {
        if (!f.eof()) {
            return float(f.tellg()) / fsize;
        } else {
            return 1;
        }
    }

    bool read(Fragment* output) override {
        if (f.eof()) {
            return false;
        }

        if (output) {
            output->spotid.assign(last_desc, 1, last_desc.size() - 1);
            output->bases.clear();
            output->bases.reserve(10000); // todo: tune
        }

		std::string line;
		while (!f.eof()) {
            read_line(line);

            if (is_description(line)) {
                last_desc = line;
                break;
            } else if (output) {
                output->bases += line;
            }
        }

        if (output) {
            if (output->bases.empty()) {
                throw std::runtime_error("Read is empty");
            }
            output->bases.shrink_to_fit();
        }
        ++spot_idx;
        return true;
    }
};

