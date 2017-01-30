#ifndef FILE_LIST_LOADER_H_INCLUDED
#define FILE_LIST_LOADER_H_INCLUDED

#include <string>
#include <fstream>
#include <list>

struct FileListLoader
{
	struct File
	{
		size_t filesize;
		std::string filename;
	};

	typedef std::vector<File> Files;
	Files files;

	FileListLoader(const std::string &file_list)
	{
		std::ifstream f(file_list);
		if (f.fail() || f.bad())
			throw std::string("cannot open file list ") + file_list;

		while (!f.eof())
		{
			File f_rec;
			f >> f_rec.filesize;
			char ch;
			f.get(ch);
//			std::getchar(f, ch);
			if (ch != '\t')
				throw "invalid separator";

			std::getline(f, f_rec.filename);
//			f >> ;

			if (f_rec.filename.empty())
				break;

//			std::cout << f_rec.filename << std::endl;
			files.push_back(f_rec);
		}
	}

};

#endif
