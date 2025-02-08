#include "json/json.h"
#include <iostream>

int main(int argc, char** argv) {
    Json::Value root;
    Json::Reader reader;
    Json::StyledWriter writer;

    std::string json_str = "{\"name\":\"test\",\"age\":20}";
    if (reader.parse(json_str, root)) {
        std::cout << writer.write(root) << std::endl;
    } else {
        std::cout << "parse error" << std::endl;
    }

    return 0;
}