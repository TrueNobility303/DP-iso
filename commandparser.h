#ifndef SUBGRAPHMATCHING_COMMANDPARSER_H
#define SUBGRAPHMATCHING_COMMANDPARSER_H

#include <string>
#include <algorithm>
#include <vector>

//Parser，将每个选项存储入tokens中
class CommandParser {
private:
    std::vector<std::string> tokens_;

public:
    CommandParser(const int argc, char **argv);
    const std::string getCommandOption(const std::string &option) const;
    bool commandOptionExists(const std::string &option) const;
};

CommandParser::CommandParser(const int argc, char **argv) {
    for (int i = 1; i < argc; ++i)
        tokens_.push_back(std::string(argv[i]));
}

//遍历所有tokens获取输入的命令
const std::string CommandParser::getCommandOption(const std::string &option) const {

    std::vector<std::string>::const_iterator itr;
    itr = find(tokens_.begin(), tokens_.end(), option);
    if (itr != tokens_.end() && ++itr != tokens_.end()) {
        return *itr;
    }
    return "";
}

bool CommandParser::commandOptionExists(const std::string &option) const {
    return find(tokens_.begin(), tokens_.end(), option) != tokens_.end();
}

#endif //SUBGRAPHMATCHING_COMMANDPARSER_H