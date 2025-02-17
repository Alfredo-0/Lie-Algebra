#ifndef PAIRUTILS_H
#define PAIRUTILS_H

#include <vector>
#include <string>
#include <sstream>
#include <cctype>
#include <fstream>
#include <iostream>


struct Pair {
    int left;
    double right;
};

inline std::string trim(const std::string &s) {
    auto start = s.begin();
    while (start != s.end() && std::isspace(*start))
        ++start;
    auto end = s.end();
    do {
        --end;
    } while (std::distance(start, end) > 0 && std::isspace(*end));
    return std::string(start, end + 1);
}

inline std::vector<Pair> parseLine(const std::string &line) {
    std::vector<Pair> pairs;
    std::string cleaned;
    for (char ch : line) {
        if (ch != '(' && ch != ')')
            cleaned.push_back(ch);
    }
    std::istringstream commaStream(cleaned);
    std::string token;
    while (std::getline(commaStream, token, ',')) {
        token = trim(token);
        if (token.empty())
            continue;
        std::istringstream tokenStream(token);
        int left;
        double right;
        tokenStream >> left >> right;
        pairs.push_back({left, right});
    }
    return pairs;
}

struct PairLists {
    std::vector<Pair> list1;
    std::vector<Pair> list2;
};

inline PairLists readPairLists(const std::string &filename) {
    PairLists result;
    std::ifstream input(filename);
    if (!input) {
        std::cerr << "Error opening file " << filename << "\n";
        return result;
    }
    std::string line;
    if (std::getline(input, line))
        result.list1 = parseLine(line);
    if (std::getline(input, line))
        result.list2 = parseLine(line);
    return result;
}

#endif // PAIRUTILS_H
