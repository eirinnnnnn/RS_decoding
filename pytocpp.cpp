// proj2.cpp
// Converted from gf.py, poly.py, decoder.py
#include <vector>
#include <map>
#include <iostream>
#include <stdexcept>
#include <algorithm>

using std::vector;
using std::map;
using std::cout;
using std::endl;

//------------------------------------------------------------------------------
// Field parameters (from encoder.py / proj2 spec)
//------------------------------------------------------------------------------

const int N = 63;
const int K = 42;
const int R = 21;

//------------------------------------------------------------------------------
// GF(2^6) arithmetic (gf.py)
//------------------------------------------------------------------------------

static const int EXP_TABLE[63] = {
    1,  2, 4, 8, 16, 
    32, 3, 6, 12, 24, 
    48, 35, 5, 10, 20, 
    40, 19, 38, 15, 30, 
    60, 59, 53, 41, 17, 
    34, 7, 14, 28, 56, 
    51, 37, 9, 18, 36, 
    11, 22, 44, 27, 54, 
    47, 29, 58, 55, 45, 
    25, 50, 39, 13, 26, 
    52, 43, 21, 42, 23, 
    46, 31, 62, 63, 61, 
    57, 49, 33
};

static int LOG_TABLE[64];
static const int GENERATOR_POLY_COEFFS[22] = {
    58, 62, 59, 7, 35, 58, 63, 47, 51, 6, 33,
    43, 44, 27, 7, 53, 39, 62, 52, 41, 44, 1
};
static const int G_EVAL_LIST[63] = {
    34, 
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,16,13,43,41,48,15,11,52,33,
    1,22,20,28,17,13,19,14,56,25,
    45,10,46,45,8,12,14,44,14,17,
    26,31,22,59,29,52,31,57,48,45,
    51,13
};

inline int gf_add(int a, int b) { return a ^ b; }
inline int gf_sub(int a, int b) { return a ^ b; }
int gf_mul(int a, int b) {
    if (!a || !b) return 0;
    return EXP_TABLE[(LOG_TABLE[a] + LOG_TABLE[b]) % 63];
}
int gf_div(int a, int b) {
    if (b == 0) throw std::runtime_error("division by zero in GF(2^6)");
    if (a == 0) return 0;
    int idx = (LOG_TABLE[a] - LOG_TABLE[b]) % 63;
    if (idx < 0) idx += 63;
    return EXP_TABLE[idx];
}
int gf_inv(int a) {
    if (a == 0) throw std::runtime_error("inverse of 0 in GF(2^6)");
    return EXP_TABLE[(63 - LOG_TABLE[a]) % 63];
}
int gf_pow(int a, int n) {
    if (a == 0) return 0;
    return EXP_TABLE[(LOG_TABLE[a] * n) % 63];
}
int poly_eval(const vector<int>& poly, int x) {
    int result = 0;
    int power = 1;
    for (int coeff : poly) {
        result = gf_add(gf_mul(coeff, power), result);
        power = gf_mul(power, x);
    }
    return result;
}
void init_gf_tables() {
    for (int i = 0; i < 63; ++i) {
        LOG_TABLE[EXP_TABLE[i]] = i;
    }
}
void verify_generator_polynomial() {
    bool passed = true;
    for (int i = 0; i < 63; ++i) {
        int alpha_i = EXP_TABLE[i];
        int expected = G_EVAL_LIST[i];
        int actual = poly_eval(
            vector<int>(GENERATOR_POLY_COEFFS, GENERATOR_POLY_COEFFS+22),
            alpha_i
        );
        if (actual != expected) {
            cout << "Mismatch at alpha^" << i
                 << ": expected " << expected
                 << ", got " << actual << "\n";
            passed = false;
        }
    }
    if (passed) cout << "All g(alpha^i) values match the official table.\n";
}

//------------------------------------------------------------------------------
// Polynomial operations (poly.py)
//------------------------------------------------------------------------------

vector<int> poly_trim(const vector<int>& p) {
    int i = (int)p.size() - 1;
    while (i >= 0 && p[i] == 0) --i;
    if (i < 0) return {0};
    return vector<int>(p.begin(), p.begin() + i + 1);
}

int poly_deg(const vector<int>& p) {
    for (int i = (int)p.size() - 1; i >= 0; --i) {
        if (p[i] != 0) return i;
    }
    return -1;
}

vector<int> poly_add(const vector<int>& f, const vector<int>& g) {
    size_t m = std::max(f.size(), g.size());
    vector<int> F(f), G(g);
    F.resize(m);
    G.resize(m);
    vector<int> out(m);
    for (size_t i = 0; i < m; ++i)
        out[i] = gf_add(F[i], G[i]);
    return out;
}

vector<int> poly_mul(const vector<int>& f, const vector<int>& g) {
    vector<int> res(f.size() + g.size() - 1, 0);
    for (size_t i = 0; i < f.size(); ++i)
        for (size_t j = 0; j < g.size(); ++j)
            res[i+j] = gf_add(res[i+j], gf_mul(f[i], g[j]));
    return res;
}

vector<int> poly_scale(const vector<int>& p, int scalar) {
    vector<int> out(p.size());
    for (size_t i = 0; i < p.size(); ++i)
        out[i] = gf_mul(p[i], scalar);
    return out;
}

vector<int> poly_shift(const vector<int>& p, int n) {
    vector<int> out(n, 0);
    out.insert(out.end(), p.begin(), p.end());
    return out;
}

std::pair<vector<int>,vector<int>>
poly_divmod(vector<int> f, vector<int> g) {
    f = poly_trim(f);
    g = poly_trim(g);
    int df = poly_deg(f), dg = poly_deg(g);
    if (dg < 0) throw std::runtime_error("division by zero");
    vector<int> quot(df - dg + 1, 0), rem(f);
    while (poly_deg(rem) >= dg) {
        int shift = poly_deg(rem) - dg;
        int lc = gf_div(rem.back(), g.back());
        vector<int> sg = poly_scale(g, lc);
        vector<int> aligned(shift, 0);
        aligned.insert(aligned.end(), sg.begin(), sg.end());
        quot[shift] = lc;
        rem = poly_trim(poly_add(rem, aligned));
    }
    return {quot, rem};
}

vector<int>
extended_euclidean(const vector<int>& a, const vector<int>& b, int mu, int nu) {
    vector<int> r_prev = poly_trim(a), r_curr = poly_trim(b);
    vector<int> u_prev = {1}, u_curr = {0};
    vector<int> v_prev = {0}, v_curr = {1};
    while (poly_deg(r_curr) > nu || poly_deg(v_curr) > mu) {
        auto [q, r_next] = poly_divmod(r_prev, r_curr);
        r_prev = r_curr; r_curr = r_next;
        vector<int> u_next = poly_trim(poly_add(u_prev, poly_mul(q, u_curr)));
        u_prev = u_curr; u_curr = u_next;
        vector<int> v_next = poly_trim(poly_add(v_prev, poly_mul(q, v_curr)));
        v_prev = v_curr; v_curr = v_next;
    }
    return v_curr;  // Note: Python returns (v_curr, r_curr) but decode_key_equation
}

//------------------------------------------------------------------------------
// Reed–Solomon encoder
//------------------------------------------------------------------------------

vector<int> rs_encode(const vector<int>& message) {
    if ((int)message.size() != K)
        throw std::runtime_error("Expected message length K");
    auto m_shifted = poly_shift(message, R);
    vector<int> gen(GENERATOR_POLY_COEFFS, GENERATOR_POLY_COEFFS+R);
    auto [quot, remainder] = poly_divmod(m_shifted, gen);
    auto parity = poly_trim(remainder);
    auto codeword = poly_trim(poly_add(parity, m_shifted));
    return codeword;
}

//------------------------------------------------------------------------------
// Reed–Solomon decoder (decoder.py)
//------------------------------------------------------------------------------

vector<int> build_erasure_locator(const vector<int>& erasures) {
    vector<int> sigma0 = {1};
    for (int i : erasures) {
        // term = (1 - α^i x) → [gf_sub(0, α^i), 1]
        int coeff = gf_sub(0, EXP_TABLE[i]);
        sigma0 = poly_mul(sigma0, { coeff, 1 });
    }
    return sigma0;
}

vector<int> erase_positions(vector<int> received, const vector<int>& erasures) {
    for (int i : erasures) received[i] = 0;
    return received;
}

vector<int> compute_syndrome(const vector<int>& rec) {
    vector<int> S(R);
    for (int j = 1; j <= R; ++j)
        S[j-1] = poly_eval(rec, EXP_TABLE[j]);
    return S;
}

vector<int> modified_syndrome(const vector<int>& S, const vector<int>& sigma0) {
    vector<int> S0 = poly_mul(sigma0, S);
    S0.resize(R);
    return S0;
}

std::pair<vector<int>,vector<int>>
solve_key_equation(const vector<int>& S0, int e0) {
    vector<int> r_poly(R+1, 0);
    r_poly[R] = 1;
    int mu = (R - e0) / 2;
    int nu = (R + e0 - 1) / 2;
    auto v_curr = extended_euclidean(r_poly, S0, mu, nu);
    vector<int> omega = poly_trim(v_curr); // but we need r_curr too; we re-run EU:
    // Let's re-run to get both:
    vector<int> r_prev = poly_trim(r_poly), r_curr = poly_trim(S0);
    vector<int> u_prev = {1}, u_curr = {0};
    vector<int> v_prev = {0}, v_curr2 = {1};
    while (poly_deg(r_curr) > nu || poly_deg(v_curr2) > mu) {
        auto [q, r_next] = poly_divmod(r_prev, r_curr);
        r_prev = r_curr; r_curr = r_next;
        vector<int> u_next = poly_trim(poly_add(u_prev, poly_mul(q, u_curr)));
        u_prev = u_curr; u_curr = u_next;
        vector<int> v_next = poly_trim(poly_add(v_prev, poly_mul(q, v_curr2)));
        v_prev = v_curr2; v_curr2 = v_next;
    }
    vector<int> sigma1 = v_curr2;
    omega = r_curr;
    // Normalize by sigma1[0]
    int s0 = sigma1[0];
    if (s0 == 0) throw std::runtime_error("sigma1(0) = 0, cannot normalize");
    for (auto &c : sigma1) c = gf_div(c, s0);
    for (auto &c : omega) c = gf_div(c, s0);
    return {sigma1, omega};
}

vector<int> combine_locators(const vector<int>& a, const vector<int>& b) {
    return poly_mul(a, b);
}

vector<int> find_error_positions(const vector<int>& sigma) {
    vector<int> pos;
    for (int i = 0; i < N; ++i) {
        int x_inv = EXP_TABLE[(63 - i) % 63];
        if (poly_eval(sigma, x_inv) == 0) pos.push_back(i);
    }
    return pos;
}

vector<int> poly_derivative(vector<int> p) {
    if (p.size() < 2) return {0};
    vector<int> d(p.size()-1, 0);
    for (size_t i = 1; i < p.size(); ++i) {
        if (i & 1) d[i-1] = p[i];
    }
    while (d.size()>1 && d.back()==0) d.pop_back();
    return d;
}

map<int,int> evaluate_error_magnitudes(
    const vector<int>& sigma,
    const vector<int>& omega,
    const vector<int>& positions
) {
    map<int,int> ev;
    vector<int> sigma_p = poly_derivative(sigma);
    for (int i : positions) {
        int x_inv = EXP_TABLE[(63 - i) % 63];
        int num = poly_eval(omega, x_inv);
        int den = poly_eval(sigma_p, x_inv);
        if (den == 0) throw std::runtime_error("σ'(...)=0 during Forney eval");
        ev[i] = gf_sub(0, gf_div(num, den));
    }
    return ev;
}

vector<int> apply_error_correction(
    vector<int> rec,
    const map<int,int>& ev
) {
    for (auto [i, mag] : ev) {
        rec[i] = gf_sub(rec[i], mag);
    }
    return rec;
}

vector<int> rs_decode(const vector<int>& received, const vector<int>& erasures) {
    auto sigma0 = build_erasure_locator(erasures);
    auto R_erased = erase_positions(received, erasures);
    auto S = compute_syndrome(R_erased);
    auto S0 = modified_syndrome(S, sigma0);
    auto [sigma1, omega] = solve_key_equation(S0, erasures.size());
    auto sigma = combine_locators(sigma0, sigma1);
    auto pos = find_error_positions(sigma);
    auto ev  = evaluate_error_magnitudes(sigma, omega, pos);
    auto corr= apply_error_correction(received, ev);
    return vector<int>(corr.end() - K, corr.end());
}

vector<int> debug_rs_decode(const vector<int>& received, const vector<int>& erasures) {
    cout << "\n🧪 BEGIN RS DECODING DEBUG 🧪\n";
    cout << "🔹 Received codeword (LSB-first): "; 
    for (int x: received) cout << x << " "; cout << "\n";
    auto sigma0 = build_erasure_locator(erasures);
    auto R_erased = erase_positions(received, erasures);
    auto S = compute_syndrome(R_erased);
    cout << "🔸 Syndrome S(x): "; 
    for (int s: S) cout << s << " "; cout << "\n";
    auto S0 = modified_syndrome(S, sigma0);
    auto [sigma1, omega] = solve_key_equation(S0, erasures.size());
    cout << "🔸 σ₁(x) from EE: "; for (int c: sigma1) cout<<c<<" "; cout<<"\n";
    cout << "🔸 ω(x) from EE:  "; for (int c: omega)    cout<<c<<" "; cout<<"\n";
    auto sigma = combine_locators(sigma0, sigma1);
    auto pos   = find_error_positions(sigma);
    cout << "🔸 Error positions found: "; for (int i: pos) cout<<i<<" "; cout<<"\n";
    auto ev    = evaluate_error_magnitudes(sigma, omega, pos);
    cout << "🔸 Error magnitudes: ";
    for (auto [i,m] : ev) cout<<"("<<i<<":"<<m<<") ";
    cout << "\n";
    auto corr = apply_error_correction(received, ev);
    cout << "✅ Corrected codeword: ";
    for (int x: corr) cout<<x<<" "; cout<<"\n";
    cout << "🧪 END RS DECODING DEBUG 🧪\n\n";
    return vector<int>(corr.end() - K, corr.end());
}

//------------------------------------------------------------------------------
// main()
//------------------------------------------------------------------------------

int main() {
    init_gf_tables();
    verify_generator_polynomial();
    // additional tests / usage here...
    return 0;
}
