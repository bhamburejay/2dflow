#include <nlohmann/json.hpp>
#include <iostream>


int main() {
    nlohmann::json j;
    j["number"] = 42;
    j["pi"] = 3.14;
    std::cout << j << std::endl;

    try {
        auto value = j.at("key").get<std::string>();
        std::cout << "Key: " << value << std::endl;
    } catch (nlohmann::json::out_of_range& e) {
        std::cerr << "Key not found: " << e.what() << std::endl;
    }
    
    return 0;
}