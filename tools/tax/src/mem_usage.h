#ifndef MEM_USAGE_H_INCLUDED
#define MEM_USAGE_H_INCLUDED

#include <fstream>
#include <string>
#include <unistd.h>

static unsigned long long int mem_usage()
{
	// http://man7.org/linux/man-pages/man5/proc.5.html
   std::ifstream f("/proc/self/stat", std::ios_base::in);
   std::string x;
   for (int i=0; i<23; i++)
	   f >> x;

   unsigned long long int mem;
   f >> mem;
   return mem * sysconf(_SC_PAGE_SIZE);
}

#endif
