// proj2.cpp
// Reed-Solomon Encoder & Decoder (RS(63,42)) over GF(2^6)
// Converted from Python to C++

#include <iostream>
#include <vector>
#include <stdexcept>
#include <map>

using namespace std;

// -------------------- GF(2^6) Constants --------------------
const int FIELD_SIZE = 64;
const int PRIMITIVE_POLY = 0x43;  // x^6 + x + 1
const int N = 63;
const int K = 42;
const int R = N - K;

const int EXP_TABLE[N] = {
    1, 2, 4, 8, 16, 32, 3, 6, 12, 24, 48, 35, 5, 10, 20,
    40, 19, 38, 15, 30, 60, 59, 53, 41, 17, 34, 7, 14, 28, 56,
    51, 37, 9, 18, 36, 11, 22, 44, 27, 54, 47, 29, 58, 55, 45,
    25, 50, 39, 13, 26, 52, 43, 21, 42, 23, 46, 31, 62, 63, 61,
    57, 49, 33
};

int LOG_TABLE[FIELD_SIZE];

void init_log_table() {
    fill(LOG_TABLE, LOG_TABLE + FIELD_SIZE, -1);
    for (int i = 0; i < N; ++i) {
        LOG_TABLE[EXP_TABLE[i]] = i;
    }
}

// -------------------- GF(2^6) Arithmetic --------------------
int gf_add(int a, int b) { return a ^ b; }
int gf_sub(int a, int b) { return a ^ b; }

int gf_mul(int a, int b) {
    if (a == 0 || b == 0) return 0;
    return EXP_TABLE[(LOG_TABLE[a] + LOG_TABLE[b]) % 63];
}

int gf_div(int a, int b) {
    if (b == 0) throw runtime_error("GF division by zero");
    if (a == 0) return 0;
    return EXP_TABLE[(LOG_TABLE[a] - LOG_TABLE[b] + 63) % 63];
}

int gf_inv(int a) {
    if (a == 0) throw runtime_error("GF inverse of zero");
    return EXP_TABLE[(63 - LOG_TABLE[a]) % 63];
}

int gf_pow(int a, int n) {
    if (a == 0) return 0;
    return EXP_TABLE[(LOG_TABLE[a] * n) % 63];
}

// -------------------- Polynomial Operations --------------------
using Poly = vector<int>;

Poly poly_add(const Poly& f, const Poly& g) {
    Poly result(max(f.size(), g.size()), 0);
    for (size_t i = 0; i < f.size(); ++i) result[i] = f[i];
    for (size_t i = 0; i < g.size(); ++i) result[i] ^= g[i];
    return result;
}

Poly poly_mul(const Poly& f, const Poly& g) {
    Poly result(f.size() + g.size() - 1, 0);
    for (size_t i = 0; i < f.size(); ++i)
        for (size_t j = 0; j < g.size(); ++j)
            result[i + j] ^= gf_mul(f[i], g[j]);
    return result;
}

Poly poly_shift(const Poly& p, int n) {
    Poly result(n, 0);
    result.insert(result.end(), p.begin(), p.end());
    return result;
}

int poly_deg(const Poly& p) {
    for (int i = p.size() - 1; i >= 0; --i)
        if (p[i] != 0) return i;
    return -1;
}

Poly poly_trim(const Poly& p) {
    int deg = poly_deg(p);
    if (deg < 0) return {0};
    return Poly(p.begin(), p.begin() + deg + 1);
}

Poly poly_scale(const Poly& p, int scalar) {
    Poly result;
    for (auto c : p) result.push_back(gf_mul(c, scalar));
    return result;
}

pair<Poly, Poly> poly_divmod(Poly f, Poly g) {
    f = poly_trim(f); g = poly_trim(g);
    Poly q(max((int)f.size() - (int)g.size() + 1, 0), 0);
    while (poly_deg(f) >= poly_deg(g)) {
        int shift = poly_deg(f) - poly_deg(g);
        int coef = gf_div(f[poly_deg(f)], g[poly_deg(g)]);
        Poly g_scaled = poly_scale(g, coef);
        Poly g_shifted(shift, 0);
        g_shifted.insert(g_shifted.end(), g_scaled.begin(), g_scaled.end());
        q[shift] = coef;
        f = poly_add(f, g_shifted);
        f = poly_trim(f);
    }
    return {q, f};
}



Poly poly_derivative(const Poly& p) {
    if (p.size() < 2) return {0};
    Poly deriv(p.size() - 1, 0);
    for (size_t i = 1; i < p.size(); i += 2) deriv[i - 1] = p[i];
    while (deriv.size() > 1 && deriv.back() == 0) deriv.pop_back();
    return deriv;
}

int poly_eval(const Poly& poly, int x) {
    int result = 0, power = 1;
    for (int c : poly) {
        result = gf_add(result, gf_mul(c, power));
        power = gf_mul(power, x);
    }
    return result;
}

Poly build_erasure_locator(const vector<int>& erasures) {
    Poly sigma0 = {1};
    for (int i : erasures) {
        Poly term = {gf_sub(0, EXP_TABLE[i]), 1};
        sigma0 = poly_mul(sigma0, term);
    }
    return sigma0;
}

vector<int> erase_positions(const vector<int>& received, const vector<int>& erasures) {
    vector<int> erased = received;
    for (int i : erasures) erased[i] = 0;
    return erased;
}

Poly compute_syndrome(const vector<int>& received) {
    Poly S;
    for (int j = 1; j <= R; ++j) S.push_back(poly_eval(received, EXP_TABLE[j % 63]));
    return S;
}

Poly modified_syndrome(const Poly& S, const Poly& sigma0) {
    Poly S0 = poly_mul(sigma0, S);
    S0.resize(R);
    return S0;
}

pair<Poly, Poly> extended_euclidean(Poly a, Poly b, int mu, int nu) {
    a = poly_trim(a);
    b = poly_trim(b);
    Poly r_prev = a, r_curr = b;
    Poly v_prev = {0}, v_curr = {1};
    while (poly_deg(r_curr) > nu || poly_deg(v_curr) > mu) {
        auto [q, r_next] = poly_divmod(r_prev, r_curr);
        r_prev = r_curr; r_curr = r_next;
        Poly v_next = poly_trim(poly_add(v_prev, poly_mul(q, v_curr)));
        v_prev = v_curr; v_curr = v_next;
    }
    return {v_curr, r_curr};
}

vector<int> find_error_positions(const Poly& sigma) {
    vector<int> positions;
    for (int i = 0; i < N; ++i) {
        int xinv = EXP_TABLE[(63 - i) % 63];
        if (poly_eval(sigma, xinv) == 0) positions.push_back(i);
    }
    return positions;
}

map<int, int> evaluate_error_magnitudes(const Poly& sigma, const Poly& omega, const vector<int>& positions) {
    map<int, int> error_vector;
    Poly sigma_prime = poly_derivative(sigma);
    for (int i : positions) {
        int xinv = EXP_TABLE[(63 - i) % 63];
        int num = poly_eval(omega, xinv);
        int denom = poly_eval(sigma_prime, xinv);
        if (denom == 0) throw runtime_error("Forney: divide by zero");
        error_vector[i] = gf_sub(0, gf_div(num, denom));
    }
    return error_vector;
}

vector<int> apply_error_correction(const vector<int>& received, const map<int, int>& error_vector) {
    vector<int> corrected = received;
    for (auto& [i, mag] : error_vector)
        corrected[i] = gf_sub(corrected[i], mag);
    return corrected;
}
vector<int> rs_decode(const vector<int>& received, const vector<int>& erasures) {
    Poly sigma0 = build_erasure_locator(erasures);
    vector<int> R_erased = erase_positions(received, erasures);
    Poly S = compute_syndrome(R_erased);
    Poly S0 = modified_syndrome(S, sigma0);

    int t0 = erasures.size();
    int mu = (R - t0) / 2;
    int nu = (R + t0 - 1) / 2;

    auto [sigma1, omega] = extended_euclidean(vector<int>(R + 1, 0), S0, mu, nu);
    if (sigma1[0] == 0) return {-1};
    int sigma1_0_inv = gf_inv(sigma1[0]);
    for (int& c : sigma1) c = gf_mul(c, sigma1_0_inv);
    for (int& c : omega) c = gf_mul(c, sigma1_0_inv);

    Poly sigma = poly_mul(sigma0, sigma1);
    vector<int> error_pos = find_error_positions(sigma);
    if (error_pos.size() != poly_deg(sigma)) return {-1};
    map<int, int> error_mag = evaluate_error_magnitudes(sigma, omega, error_pos);
    if (error_mag.count(-1)) return {-1};
    vector<int> corrected = apply_error_correction(received, error_mag);

    if (any_of(compute_syndrome(corrected).begin(), compute_syndrome(corrected).end(), [](int x) { return x != 0; }))
        return {-1};

    return vector<int>(corrected.end() - K, corrected.end());
}