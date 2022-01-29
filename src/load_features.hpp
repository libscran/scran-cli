#ifndef LOAD_FEATURES_H
#define LOAD_FEATURES_H

#include <type_traits>
#include <vector>
#include <string>
#include <cstdio>
#include "buffin/parse_text_file.hpp"
#include "buffin/parse_gzip_file.hpp"

template<class Function>
struct LineProcessor {
    LineProcessor(Function f) : fun(std::move(f)) {}
    
    template<typename B>
    void add (const B* buffer, size_t n) {
        size_t last = 0;
        size_t i = 0;
        while (i < n) {
            if (buffer[i] == '\n') {
                if (continuing) {
                    collected.push_back(fun(line.c_str(), line.c_str() + line.size()));
                    continuing = false;
                    line = "";
                } else {
                    collected.push_back(fun(buffer + last, buffer + i));
                }
                last = i + 1;
            }
            ++i;
        }

        if (last != n) {
            if (continuing) {
                line += std::string(buffer + last, buffer + n);
            } else {
                continuing = true;
                line = std::string(buffer + last, buffer + n);
            }
        }
    }

    void finish() {
        if (line != "") {
            collected.push_back(fun(line.c_str(), line.c_str() + line.size()));
        }
    }

    bool continuing = false;
    std::string line;

    Function fun;
    std::vector<typename std::invoke_result<Function, const char*, const char*>::type> collected;
};

std::vector<std::pair<std::string, std::string> > load_features(std::string path) {
    auto reader = [](auto start, auto end) -> std::pair<std::string, std::string> {
        auto copy = start;
        while (copy != end && *copy != '\t') {
            ++copy;
        }
        std::string id(start, copy);

        if (copy != end) {
            ++copy; // skip the tab.
        }

        start = copy;
        while (copy != end && *copy != '\t') {
            ++copy;
        }
        std::string sym(start, copy);

        return std::make_pair(id, sym);
    };

    char header[3];
    size_t read;
    {
        buffin::SelfClosingFILE file(path.c_str());
        read = std::fread(header, sizeof(char), 3, file.handle);
    }

    LineProcessor parser(reader);
    if (read == 3 && static_cast<unsigned char>(header[0]) == 0x1f && static_cast<unsigned char>(header[1]) == 0x8b && static_cast<unsigned char>(header[2]) == 0x08) {
        buffin::parse_gzip_file(path.c_str(), parser);
    } else {
        buffin::parse_text_file(path.c_str(), parser);
    }

    parser.finish();
    return parser.collected;
}

#endif
