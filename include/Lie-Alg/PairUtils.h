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

// Updated parseLine: now returns a vector of vectors of Pair.
// Each comma-separated segment in the input is parsed into a vector of pairs.
inline std::vector<std::vector<Pair>> parseLine(const std::string &line) {
    std::vector<std::vector<Pair>> result;
    std::string cleaned;
    // Remove parentheses from the line.
    for (char ch : line) {
        if (ch != '(' && ch != ')')
            cleaned.push_back(ch);
    }
    
    std::istringstream commaStream(cleaned);
    std::string token;
    // Process each comma-separated token.
    while (std::getline(commaStream, token, ',')) {
        token = trim(token);
        if (token.empty())
            continue;
        
        std::istringstream tokenStream(token);
        std::vector<Pair> pairs;
        int left;
        double right;
        // Extract pairs: every two numbers are taken as a pair.
        while (tokenStream >> left >> right) {
            pairs.push_back({left, right});
        }
        result.push_back(pairs);
    }
    
    return result;
}

// Updated PairLists: each list is now a vector of vectors of Pair.
struct PairLists {
    std::vector<std::vector<Pair>> list1;
    std::vector<std::vector<Pair>> list2;
};

// Reads two lines and converts them into PairLists using the new parseLine.
inline PairLists readPairLists(const std::string &line1, const std::string &line2) {
    PairLists result;
    result.list1 = parseLine(line1);
    result.list2 = parseLine(line2);
    return result;
}

#endif // PAIRUTILS_H
