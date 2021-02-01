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

#ifndef CONFIG_H_INCLUDED
#define CONFIG_H_INCLUDED

#include <string>
#include <iostream>

struct Config
{
    std::string file_list, dbss_filename;
    int start_window_len = 0, max_window_len = 0, test_windows_count = 0, min_test_windows_count = 0;
    float needed_part_identified = 0, window_upscale_factor = 0;

    Config(int argc, char const *argv[])
    {
        if (argc != 9)
        {
            print_usage();
            exit(1);
        }

        file_list = std::string(argv[1]);
        dbss_filename = std::string(argv[2]);
        start_window_len = std::stoi(std::string(argv[3]));
        max_window_len = std::stoi(std::string(argv[4]));
        test_windows_count = std::stoi(std::string(argv[5]));
        min_test_windows_count = std::stoi(std::string(argv[6]));

        needed_part_identified = std::stof(std::string(argv[7]));
        window_upscale_factor = std::stof(std::string(argv[8]));

        if (window_upscale_factor < 1.0)
            throw std::runtime_error("window_upscale_factor < 1.0");

        if (needed_part_identified < 0 || needed_part_identified > 1.0f)
            throw std::runtime_error("needed_part_identified < 0 || needed_part_identified > 1.0f");

        if (test_windows_count < min_test_windows_count)
            throw std::runtime_error("test_windows_count < min_test_windows_count");

        if (max_window_len < start_window_len)
            throw std::runtime_error("max_window_len < start_window_len");

        if (start_window_len < 0)
            throw std::runtime_error("start_window_len < 0");
    }

    static void print_usage()
    {
        std::cerr << "need <files.list> <dbss file> <start window len> <max window len> <test windows count> <min test windows count> <needed part identified> <window upscale factor>" << std::endl;
    }
};

#endif
