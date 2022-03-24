#ifndef LOAD_FEATURES_H
#define LOAD_FEATURES_H

#include <type_traits>
#include <vector>
#include <string>
#include <cstdio>
#include "byteme/SomeFileReader.hpp"

std::vector<std::pair<std::string, std::string> > load_features(std::string path) {
    byteme::SomeFileReader reader(path);
    bool remaining = true;

    bool continuing = false;
    std::string line;
    std::vector<std::pair<std::string, std::string> > collected;

    auto fun = [&](auto start, auto end) -> void {
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

        collected.emplace_back(std::move(id), std::move(sym));
    };

    do {
        remaining = reader();
        auto buffer = reinterpret_cast<const char*>(reader.buffer());
        auto n = reader.available();

        size_t last = 0;
        size_t i = 0;
        while (i < n) {
            if (buffer[i] == '\n') {
                if (continuing) {
                    fun(line.c_str(), line.c_str() + line.size());
                    continuing = false;
                    line = "";
                } else {
                    fun(buffer + last, buffer + i);
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
    } while (remaining);

    if (line != "") {
        fun(line.c_str(), line.c_str() + line.size());
    }

    return collected;
}

#endif
