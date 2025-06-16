// ─────────────────── generate_io.cpp ───────────────────
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <set>
#include <string>
#include <vector>

/* bring in your encoder + Poly definition
 * (replace with    #include "proj2.h"    if you have a header) */
#include "ouuo0616.cpp"


/* random helpers ----------------------------------------------------------- */
static std::mt19937 rng(std::random_device{}());
static std::uniform_int_distribution<int> dist_symbol(0, 63);
int rand_sym() { return dist_symbol(rng); }

/* ------------------------------------------------------------------------- */
int main(int argc, char* argv[])
{
    int NUM = 20;                          // default # of test lines
    if (argc > 1) NUM = std::stoi(argv[1]);

    std::ofstream fin("input.txt");
    std::ofstream fans("answer.txt");
    if (!fin || !fans) {
        std::cerr << "cannot open output files\n";
        return 1;
    }

    std::uniform_int_distribution<int> dist_t0(0, 5);   // erasures
    std::uniform_int_distribution<int> dist_t1(0, 10);  // errors
    std::uniform_int_distribution<int> dist_pos(0, N - 1);

    for (int line = 0; line < NUM; ++line) {
        /* 1. make a random message and its codeword ------------------------ */
        Poly msg(K);
        std::generate(msg.begin(), msg.end(), rand_sym);
        Poly cw = rs_encode(msg);                         // your encoder

        /* 2. choose #erasures / #errors ------------------------- */
        bool within_radius = (line % 1 == 0);             // even lines OK
        int t0, t1;
        do {
            t0 = dist_t0(rng);
            t1 = dist_t1(rng);
        } while ((t0 + 2 * t1 <= R) != within_radius);    // enforce class

        /* 3. pick distinct positions for erasures + errors ----------------- */
        std::set<int> pos_set;
        while (static_cast<int>(pos_set.size()) < t0 + t1)
            pos_set.insert(dist_pos(rng));
        auto it = pos_set.begin();
        std::vector<int> erasures, errors;
        for (int i = 0; i < t0; ++i, ++it) erasures.push_back(*it);
        for (int i = 0; i < t1; ++i, ++it) errors.push_back(*it);

        /* 4. build the *garbled* vector ------------------------------------ */
        Poly garbled = cw;
        for (int p : errors) {                // flip symbol at each error pos
            int v;
            do { v = rand_sym(); } while (v == garbled[p]);
            garbled[p] = v;
        }

        /* 5. write INPUT line --------------------------------------------- */
        for (int j = 0; j < N; ++j) {
            if (std::find(erasures.begin(), erasures.end(), j) != erasures.end())
                fin << '*';
            else
                fin << garbled[j];
            fin << (j + 1 == N ? '\n' : ' ');
        }

        /* 6. write ANSWER line -------------------------------------------- */
        if (within_radius) {                  // should decode → print message
            for (int j = 0; j < K; ++j)
                fans << msg[j] << (j + 1 == K ? '\n' : ' ');
        } else {                              // beyond radius
            fans << -1 << '\n';
        }
    }

    std::cout << "Generated input.txt and answer.txt with "
              << NUM << " cases.\n";
    return 0;
}
