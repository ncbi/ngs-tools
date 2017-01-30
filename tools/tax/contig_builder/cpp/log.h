#ifndef LOG_H_INCLUDED
#define LOG_H_INCLUDED

void print_current_time()
{
	auto t = std::time(nullptr);
	auto timeinfo = std::localtime(&t);
	const int BUFFER_SIZE = 256;
	char buffer[BUFFER_SIZE];
	strftime(buffer, BUFFER_SIZE, "%m/%d/%Y %H:%M:%S", timeinfo);
	std::cerr << "time is " << buffer << std::endl;
}


#endif
