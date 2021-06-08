#ifndef _INPUT_H_
#define _INPUT_H_

#include <nlohmann/json.hpp>
using json = nlohmann::json;
#include <matvec.h>

#include <iomanip>	//ONY FOR GCC!!


template <typename T>
bool readValue(const nlohmann::json &j, T &v)
{
	if (j.is_null())
		return false;

	v = j.get<T>();
	return true;
}

// template <typename T, int size>
// bool readVector(const nlohmann::json &j, Eigen::Matrix<T, size, 1, Eigen::DontAlign> &vec)
// {
	// if (j.is_null())
		// return false;

	// std::vector<T> values = j.get<std::vector<T>>();
	// for (unsigned int i = 0; i < values.size(); i++)
		// vec[i] = values[i];
	// return true;
// }

//template <typename T>
bool readVector(const nlohmann::json &j, Vec3_t &vec)
{
	if (j.is_null())
		return false;

	std::vector<double> values = j.get< std::vector<double> >();
	for (unsigned int i = 0; i < values.size(); i++)
		vec[i] = values[i];
	return true;
}

bool readArray(const nlohmann::json &j, std::vector<double> &vec)
{
	if (j.is_null())
		return false;

	std::vector<double> values = j.get< std::vector<double> >();
	vec.resize( values.size() );
	for (unsigned int i = 0; i < values.size(); i++)
		vec[i] = values[i];
	return true;
}

#endif
