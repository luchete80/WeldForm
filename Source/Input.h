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

bool readBoolVector(const nlohmann::json &j, bool vec[])
{
	if (j.is_null())
		return false;

	std::vector<bool> values = j.get< std::vector<bool> >();
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

namespace SPH {

//////////////////////////////////////////////////////////////////////
////////////////////// BASED ON SPLISH SPLASH SCENE LOADER ///////////

class SceneLoader
{
protected:
  nlohmann::json m_jsonData;
  
	public:
		/** \brief Struct for an AABB */
		struct Box
		{
			Vec3_t m_minX;
			Vec3_t m_maxX;
		};

		/** \brief Struct to store a fluid object */
		struct FluidData
		{
			std::string id;
			std::string samplesFile;
			Vec3_t translation;
			Vec3_t rotation;
			Vec3_t scale;
			Vec3_t initialVelocity;
			//unsigned char mode;
			//bool invert;
			//std::array<unsigned int, 3> resolutionSDF;
		};

		/** \brief Struct to store a fluid block */
		struct FluidBlock
		{
			std::string id;
			Box box;
			unsigned char mode;
			//Vector3r initialVelocity;
		};

		struct MaterialData
		{
			std::string id;

		};
    
    //TODO: OUTPUT step CONTROL!
    
		/** \brief Struct to store scene information */
		// struct Scene
		// {
			// std::vector<BoundaryData*> boundaryModels;
			// std::vector<FluidData*> fluidModels;
			// std::vector<FluidBlock*> fluidBlocks;
			// std::vector<EmitterData*> emitters;
			// std::vector<AnimationFieldData*> animatedFields;
			// std::vector<MaterialData*> materials;
			// Real particleRadius;
			// bool sim2D;
			// Real timeStepSize;
			// Vector3r camPosition;
			// Vector3r camLookat;
		// };

}; //Scene Loader
  
}; //SPH

#endif
