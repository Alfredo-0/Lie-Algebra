#ifndef PAIRUTILS_H
#define PAIRUTILS_H

#include <ginac/ginac.h>

#include <vector>
#include <string>
#include <sstream>
#include <cctype>
#include <fstream>
#include <iostream>

struct Pair {
    int left;
    GiNaC::ex right; //
};

struct PairLists {
    std::vector<std::vector<Pair>> list1;
    std::vector<std::vector<Pair>> list2;
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

inline std::vector<std::vector<Pair>> parseLine(const std::string &line, const GiNaC::ex &x) {
    std::vector<std::vector<Pair>> result;
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
        std::vector<Pair> pairs;
        int left;
        std::string right; //

        while (tokenStream >> left >> right) {
            GiNaC::symtab table;
            table["x"] = x;
            GiNaC::parser reader(table);
            GiNaC::ex right_expr = reader(right);
            pairs.push_back({left, right_expr});
        }
        result.push_back(pairs);
    }
    
    return result;
}

inline PairLists readPairLists(const std::string &line1, const std::string &line2, const GiNaC::ex &x) {
    PairLists result;

    result.list1 = parseLine(line1, x);
    result.list2 = parseLine(line2, x);
    return result;
}

#endif // PAIRUTILS_H
