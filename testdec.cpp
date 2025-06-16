// testdecod.cpp
// Tests for RS(63,42) encoder/decoder using proj2.cpp

#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <numeric>
#include <random>
#include <cassert>
#include "proj2.cpp"  // includes rs_encode, debug_rs_decode, etc.

using std::vector;
using std::map;
using std::cout;
using std::endl;

vector<int> introduce_errors(vector<int> cw, const vector<int>& pos, const vector<int>& vals) {
    for (size_t i = 0; i < pos.size(); ++i) {
        cw[pos[i]] = gf_add(cw[pos[i]], vals[i]);
    }
    return cw;
}

vector<int> introduce_erasures(vector<int> cw, const vector<int>& erasures) {
    for (int i : erasures) cw[i] = 0;
    return cw;
}

void test_no_errors() {
    std::mt19937 rng(std::random_device{}());
    vector<int> msg(K);
    for (int& x : msg) x = std::uniform_int_distribution<int>(1, 62)(rng);
    auto cw = rs_encode(msg);
    cout << cw << endl;
    auto decoded = debug_rs_decode(cw, {});
    cout << decoded << endl;
    assert(decoded == msg && "No-error decoding failed");
    cout << "✅ No-error decoding passed" << endl;
}

void test_with_errors() {
    std::mt19937 rng(std::random_device{}());
    vector<int> msg(K);
    for (int& x : msg) x = std::uniform_int_distribution<int>(1, 62)(rng);
    auto cw = rs_encode(msg);
    // pick 5 distinct error positions
    vector<int> idx(N);
    std::iota(idx.begin(), idx.end(), 0);
    std::shuffle(idx.begin(), idx.end(), rng);
    vector<int> epos(idx.begin(), idx.begin()+5);
    vector<int> evals(5);
    for (int& v : evals) v = std::uniform_int_distribution<int>(1, 62)(rng);
    auto corrupted = introduce_errors(cw, epos, evals);
    auto decoded = debug_rs_decode(corrupted, {});
    cout << "Error positions: "; for(int p: epos) cout<<p<<" "; cout<<endl;
    cout << "Error values:    "; for(int v: evals) cout<<v<<" "; cout<<endl;
    assert(decoded == msg && "Error-only decoding failed");
    cout << "✅ Error-only decoding passed" << endl;
}

void test_with_erasures() {
    std::mt19937 rng(std::random_device{}());
    vector<int> msg(K);
    for (int& x : msg) x = std::uniform_int_distribution<int>(1, 62)(rng);
    auto cw = rs_encode(msg);
    vector<int> idx(N);
    std::iota(idx.begin(), idx.end(), 0);
    std::shuffle(idx.begin(), idx.end(), rng);
    vector<int> erasures(idx.begin(), idx.begin()+5);
    auto corrupted = introduce_erasures(cw, erasures);
    auto decoded = debug_rs_decode(corrupted, erasures);
    assert(decoded == msg && "Erasure-only decoding failed");
    cout << "✅ Erasure-only decoding passed" << endl;
}

void test_mixed_errors_erasures() {
    std::mt19937 rng(std::random_device{}());
    vector<int> msg(K);
    for (int& x : msg) x = std::uniform_int_distribution<int>(1, 62)(rng);
    auto cw = rs_encode(msg);
    vector<int> idx(N);
    std::iota(idx.begin(), idx.end(), 0);
    std::shuffle(idx.begin(), idx.end(), rng);
    vector<int> erasures(idx.begin(), idx.begin()+2);
    vector<int> rest;
    for (int i = 0; i < N; ++i) if (std::find(erasures.begin(), erasures.end(), i) == erasures.end()) rest.push_back(i);
    std::shuffle(rest.begin(), rest.end(), rng);
    vector<int> epos(rest.begin(), rest.begin()+2);
    vector<int> evals(2);
    for (int& v : evals) v = std::uniform_int_distribution<int>(1, 62)(rng);
    auto corrupted = introduce_errors(cw, epos, evals);
    corrupted = introduce_erasures(corrupted, erasures);
    auto decoded = debug_rs_decode(corrupted, erasures);
    assert(decoded == msg && "Mixed error-erasure decoding failed");
    cout << "✅ Mixed error-erasure decoding passed" << endl;
}

int main() {
    cout << "test start" << endl;
    test_no_errors();
    test_with_errors();
    test_with_erasures();
    test_mixed_errors_erasures();
    return 0;
}
