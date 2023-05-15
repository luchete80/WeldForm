#pragma once
#include <fstream>
#include <iostream>

#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include <iomanip>	//ONY FOR GCC!!

// namespace nlohmann {
    // template <>
    // struct adl_serializer<SensorData> {
        // // note: the return type is no longer 'void', and the method only takes
        // // one argument
        // static SensorData from_json(const json& j) {
            // SensorData sensor_data;
            // sensor_data.sensor_identify = j.at("SensorIdentify").get<string>();
            // return sensor_data;
        // }

        // // Here's the catch! You must provide a to_json method! Otherwise you
        // // will not be able to convert move_only_type to json, since you fully
        // // specialized adl_serializer on that type
        // static void to_json(json& j, SensorData t) {
            // j = json{ {"SensorIdentify",t.sensor_identify},{"SensorType",t.sensor_type },{"Data",t.data} };
        // }
    // };
    // template <>
    // struct adl_serializer<RfidInfo> {
        // // note: the return type is no longer 'void', and the method only takes
        // // one argument
        // static RfidInfo from_json(const json& j) {
            // RfidInfo rfid_info;
            // rfid_info.identify = j.at("Identify").get<string>();
            // return rfid_info;
        // }

        // // Here's the catch! You must provide a to_json method! Otherwise you
        // // will not be able to convert move_only_type to json, since you fully
        // // specialized adl_serializer on that type
        // static void to_json(json& j, RfidInfo t) {
            // j = json{ { "Identify",t.identify },{ "Position",0 } };
        // }
    // };
// }

template <typename T>
bool readValue(const nlohmann::json &j, T &v)
{
	if (j.is_null())
		return false;

	v = j.get<T>();
	return true;
}


int main(){

	std::ifstream i("input.json");
	json j;
	i >> j;
	
	//std::cout << j{"Configuration"}["pause"] << std::endl;    // returns boolean
	
	nlohmann::json config = j["Configuration"];
	
	double ts;
	//config["timeStepSize"].get<ts>;
	readValue(config["timeStepSize"], /*scene.timeStepSize*/ts);
	std::cout << "Time Step size is: "<<ts<<std::endl;
	// for(auto &array : j["objList"]) {
    // std::cout << array["key1"] << std::endl;    // returns string
    // std::cout << array["key2"] << std::endl;    // returns array
	// }

	// write prettified JSON to another file
	std::ofstream o("pretty.json");
	o << std::setw(4) << j << std::endl;
	
}
