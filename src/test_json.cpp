#include "json/json.h"
#include <iostream>

Json::Value checkin(Json::Value &value){
    if (value.isNull()) {
        std::cerr << "Key not found in inputs" << std::endl;
    }
    return value;
}   

int main(int argc, char** argv) {
    Json::Value root;
    Json::Reader reader;
    Json::StyledWriter writer;

    std::string json_str = "{\"name\":\"test\",\"age\":20,\"nextlevel\": {\"name\":\"test2\",\"age\":21}}";

    if (reader.parse(json_str, root)) {
        std::cout << writer.write(root) << std::endl;
    } else {
        std::cout << "parse error" << std::endl;
    }
    Json::Value root2 = root ;
    root2["nextlevel"]["name"]  = "test3" ;
    std::cout << root2 << std::endl;
    std::cout << root << std::endl;

    // std::cout << get_input({"nextlevel", "name"}, root) << std::endl;
    // std::cout << get_input({"nextlevel2", "name"}, root) << std::endl;
    // std::cout << get_input({"nextlevel", "name2"}, root) << std::endl;

    std::cout << checkin(root["nextlevel"]["name"]).asString() << std::endl;
    std::cout << checkin(root["nextlevel"]["name2"]).asString() << std::endl;





    return 0;
}
