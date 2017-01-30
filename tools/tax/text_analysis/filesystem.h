#ifndef FILESYSTEM_H_INCLUDED
#define FILESYSTEM_H_INCLUDED

#include <list>
#include "kmers.h"
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>

struct FileSystem
{
	static std::list<std::string> get_filenames_by_mask(const std::string &file_mask)
	{
		std::list<std::string> filenames;

		list_dir_files(dir_of(file_mask), [&](const std::string &filename)
			{
//				std::cout << filename << std::endl;
				if (matches_mask(no_dir(filename), no_dir(file_mask)))
					filenames.push_back(filename);
			});

		return std::move(filenames);
	}

	template <class Lambda>
	static void list_dir_files(const std::string &dir, Lambda &&lambda)
	{
		DIR *dp;
		struct dirent *dirp;
		if(!(dp = opendir(dir.c_str()))) 
			throw std::string("error opening dir ") + dir;

		while (dirp = readdir(dp))
			lambda( dir + std::string(dirp->d_name) );

		closedir(dp);
	}

	static std::string no_dir(const std::string &filename)
	{
		auto sep_pos = filename.find_last_of('/'); //std::find(filename.rbegin(), filename.rend(), '/');
		if (sep_pos == std::string::npos)
			return filename;

		return filename.substr(sep_pos + 1);
	}

	static std::string dir_of(const std::string &filename)
	{
		auto sep_pos = filename.find_last_of('/'); //std::find(filename.rbegin(), filename.rend(), '/');
		if (sep_pos == std::string::npos)
			return "./";

		return filename.substr(0, sep_pos + 1);
	}

	static bool matches_mask(const std::string &filename, const std::string &mask)
	{
		auto star_pos = mask.find_first_of('*');
		auto star2_pos = mask.find_last_of('*');
		if (star_pos == std::string::npos)
			throw mask + " is not a mask";

		if (star_pos != star2_pos)
			throw mask + " - mask unsupported. sorry :(";

		auto first_part = mask.substr(0, star_pos);
		auto second_part = mask.substr(star_pos + 1);
		return (first_part.empty() || filename.find(first_part) == 0) && (second_part.empty() || filename.find(second_part) != std::string::npos); // not completely correct
	}
};

#endif