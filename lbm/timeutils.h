#pragma once

#include <chrono>
#include <stdio.h>

static void timestamp(char *str)
{
	std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();

	// convert from C++11 API to the C struct (milliseconds are dropped)
	std::time_t secs = std::chrono::system_clock::to_time_t(now);
	std::tm* date = std::localtime(&secs);

	// number of milliseconds since the last whole seconds (must be in C++11)
	auto last_sec = std::chrono::system_clock::from_time_t(std::mktime(date));
	int milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(now - last_sec).count();

	sprintf(str, "%04d|%02d|%02d %02d:%02d:%02d.%03d", date->tm_year+1900, date->tm_mon+1, date->tm_mday, date->tm_hour, date->tm_min, date->tm_sec, milliseconds);
}
