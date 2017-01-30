void print_current_time()
{
	auto t = std::time(nullptr);
	auto timeinfo = std::localtime(&t);
	const int BUFFER_SIZE = 256;
	char buffer[BUFFER_SIZE];
	strftime(buffer, BUFFER_SIZE, "%m/%d/%Y %H:%M:%S", timeinfo);
	cerr << "time is " << buffer << endl;
}

