// ─────────────── demo.cpp ───────────────
#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

/* Bring in Poly, rs_decode(), K, etc. */
#include "proj2.cpp"          //  ← replace with  #include "proj2.h"  if you have one
// #include "proj2.cpp"          //  ← replace with  #include "proj2.h"  if you have one

/* helper: did rs_decode signal failure? */
bool is_decode_fail(const Poly& p)
{
    return p.size() == 1 && p[0] == -1;
}

/* ------------------------------------------------------------------------- */
int main(int argc, char* argv[])
{

    std::string inFile  = argc > 1 ? argv[1] : "input.txt";
    std::string outFile = argc > 2 ? argv[2] : "output.txt";

    std::ifstream fin(inFile);
    if (!fin) { std::cerr << "Cannot open " << inFile << '\n'; return 1; }

    std::ofstream fout(outFile);
    if (!fout) { std::cerr << "Cannot open " << outFile << '\n'; return 1; }

    std::string line;
    while (std::getline(fin, line)) {

        /* skip blank or comment lines */
        auto not_space = [](unsigned char c){ return !std::isspace(c); };
        if (std::find_if(line.begin(), line.end(), not_space) == line.end()
            || line[std::find_if(line.begin(), line.end(), not_space) - line.begin()] == '#')
            continue;

        /* ---------- tokenize the line into received[] and erasures[] -------- */
        Poly received;
        std::vector<int> erasures;
        std::istringstream iss(line);
        std::string tok;
        int idx = 0;
        while (iss >> tok) {
            if (tok == "*") {                 // erasure
                received.push_back(0);        // value won't be used
                erasures.push_back(idx);
            } else {                          // regular symbol
                received.push_back(std::stoi(tok) & 0x3F);
            }
            ++idx;
        }
        if (received.empty()) continue;       // nothing useful on this line

        /* ---------- decode -------------------------------------------------- */
        Poly decoded = rs_decode(received, erasures);

        /* ---------- output -------------------------------------------------- */
        if (is_decode_fail(decoded)) {
            fout << -1 << '\n';
        } else {
            for (size_t i = 0; i < decoded.size(); ++i) {
                fout << decoded[i] << (i + 1 == decoded.size() ? '\n' : ' ');
            }
        }
    }

    std::cout << "Finished.  Results written to " << outFile << '\n';
    return 0;
}
