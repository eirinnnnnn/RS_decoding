/*  proj2.cpp  ―  (63, 42) Reed-Solomon encoder/decoder
 *  ───────────────────────────────────────────────────
 *
 */

#include <bits/stdc++.h>
using namespace std;

/* ──────────────────────────────────────────────
   GF(2⁶) constants and tables
   ────────────────────────────────────────────── */
static const int N = 63;          // code length
static const int K = 42;          // message length
static const int R = N - K;       // # parity symbols (21)

/* EXP_TABLE[i] == α^i  (mod α⁶+α+1)  for i=0..62 */
static const int EXP_TABLE[63] = {
    1,  2,  4,  8, 16,
   32,  3,  6, 12, 24,
   48, 35,  5, 10, 20,
   40, 19, 38, 15, 30,
   60, 59, 53, 41, 17,
   34,  7, 14, 28, 56,
   51, 37,  9, 18, 36,
   11, 22, 44, 27, 54,
   47, 29, 58, 55, 45,
   25, 50, 39, 13, 26,
   52, 43, 21, 42, 23,
   46, 31, 62, 63, 61,
   57, 49, 33
};


static const int GENERATOR_POLY_COEFFS[R+1] = {
    58, 62, 59,  7, 35, 58, 63, 47, 51,  6, 33,
    43, 44, 27,  7, 53, 39, 62, 52, 41, 44,  1
};

static int LOG_TABLE[64];
struct _InitTables {
    _InitTables() {
        fill(begin(LOG_TABLE), end(LOG_TABLE), -1);
        for (int i = 0; i < 63; ++i) LOG_TABLE[EXP_TABLE[i]] = i;
    }
} _initTables;


/* ──────────────────────────────────────────────
   GF(2⁶) arithmetic (⊕ is XOR, characteristic-2)
   ────────────────────────────────────────────── */
inline int gf_add(int a, int b) { return a ^ b; }
inline int gf_sub(int a, int b) { return a ^ b; }

inline int gf_mul(int a, int b) {
    if (!a || !b) return 0;
    return EXP_TABLE[(LOG_TABLE[a] + LOG_TABLE[b]) % 63];
}
inline int gf_div(int a, int b) {
    if (!b) throw runtime_error("GF div by 0");
    if (!a) return 0;
    return EXP_TABLE[(LOG_TABLE[a] - LOG_TABLE[b] + 63) % 63];
}
inline int gf_inv(int a) {
    if (!a) throw runtime_error("GF inv(0)");
    return EXP_TABLE[(63 - LOG_TABLE[a]) % 63];
}
inline int gf_pow(int a, int n) {
    if (!a) return 0;
    return EXP_TABLE[(LOG_TABLE[a] * n) % 63];
}

/* ──────────────────────────────────────────────
   Polynomial helpers  (vectors are LSB-first)
   ────────────────────────────────────────────── */
using Poly = vector<int>;

inline int poly_deg(const Poly& p) {
    for (int i = (int)p.size()-1; i >= 0; --i)
        if (p[i]) return i;
    return -1;
}
inline void poly_trim(Poly& p) {
    while (p.size() > 1 && p.back() == 0) p.pop_back();
}
inline Poly poly_add(const Poly& f, const Poly& g) {
    Poly res(max(f.size(), g.size()), 0);
    for (size_t i = 0; i < res.size(); ++i) {
        int a = (i < f.size()? f[i] : 0);
        int b = (i < g.size()? g[i] : 0);
        res[i] = gf_add(a, b);
    }
    poly_trim(res); return res;
}
inline Poly poly_scale(const Poly& p, int s) {
    Poly res(p.size(), 0);
    for (size_t i = 0; i < p.size(); ++i) res[i] = gf_mul(p[i], s);
    return res;
}
inline Poly poly_shift(const Poly& p, int n) {
    Poly res(n, 0); res.insert(res.end(), p.begin(), p.end());
    return res;
}
inline Poly poly_mul(const Poly& f, const Poly& g) {
    Poly res(f.size() + g.size() - 1, 0);
    for (size_t i = 0; i < f.size(); ++i)
        for (size_t j = 0; j < g.size(); ++j)
            res[i+j] = gf_add(res[i+j], gf_mul(f[i], g[j]));
    poly_trim(res); return res;
}

/* Polynomial division (LSB-first)  :  f = q·g + r */
pair<Poly,Poly> poly_divmod(Poly f, Poly g) {
    poly_trim(f); poly_trim(g);
    int deg_g = poly_deg(g);
    if (deg_g < 0) throw runtime_error("poly div by zero");
    Poly q(max(1, poly_deg(f)-deg_g+1), 0);
    Poly r = f;
    while (poly_deg(r) >= deg_g) {
        int shift = poly_deg(r) - deg_g;
        int lc = gf_div(r.back(), g.back());
        Poly scaled = poly_scale(g, lc);
        Poly aligned(shift, 0); aligned.insert(aligned.end(), scaled.begin(), scaled.end());
        if ((size_t)shift >= q.size()) q.resize(shift+1,0);
        q[shift] = lc;
        r = poly_add(r, aligned);   // subtraction == addition in GF(2)
    }
    poly_trim(q); poly_trim(r); return {q, r};
}

/* Polynomial evaluation */
inline int poly_eval(const Poly& coeff, int x) {
    int res = 0, power = 1;
    for (int c : coeff) {
        res = gf_add(res, gf_mul(c, power));
        power = gf_mul(power, x);
    }
    return res;
}

/* Formal derivative in char 2 */
Poly poly_derivative(const Poly& p) {
    if (p.size() < 2) return {0};
    Poly d(p.size()-1, 0);
    for (size_t i = 1; i < p.size(); ++i)
        if (i & 1) d[i-1] = p[i];
    poly_trim(d); return d;
}

/* ──────────────────────────────────────────────
   Extended Euclidean algorithm (for key equation)
   stopping when  deg(r) ≤ ν  &&  deg(v) ≤ μ
   ────────────────────────────────────────────── */
pair<Poly,Poly> extended_euclidean(
        Poly a, Poly b, int mu, int nu)
{
    Poly r_prev = a, r_curr = b;
    Poly u_prev = {1}, u_curr = {0};
    Poly v_prev = {0}, v_curr = {1};

    while (poly_deg(r_curr) > nu || poly_deg(v_curr) > mu) {
        auto [q, r_next] = poly_divmod(r_prev, r_curr);
        r_prev = r_curr;  r_curr = r_next;

        Poly u_next = poly_add(u_prev, poly_mul(q, u_curr));
        Poly v_next = poly_add(v_prev, poly_mul(q, v_curr));

        u_prev = u_curr;  u_curr = u_next;
        v_prev = v_curr;  v_curr = v_next;
    }
    poly_trim(v_curr); poly_trim(r_curr);
    return {v_curr, r_curr};      // (σ₁, ω)
}

/* ──────────────────────────────────────────────
   RS(63,42)  ENCODER
   ────────────────────────────────────────────── */
Poly rs_encode(const Poly& message) {
    if ((int)message.size() != K)
        throw runtime_error("message length != 42");
    Poly m_shift = poly_shift(message, R);                // I(x)·x^r
    auto [quot, rem] = poly_divmod(m_shift, Poly(GENERATOR_POLY_COEFFS, GENERATOR_POLY_COEFFS+R+1));
    Poly codeword = poly_add(rem, m_shift);               // r(x) + I(x)x^r
    poly_trim(codeword);
    return codeword;                                      // [parity, message]
}

/* ──────────────────────────────────────────────
   Helper routines for decoding
   ────────────────────────────────────────────── */
Poly build_erasure_locator(const vector<int>& erasures) {
    Poly sigma0 = {1};
    for (int i : erasures) {
        Poly term = { gf_sub(0, EXP_TABLE[i]), 1 };       // (1 - α^i x)
        sigma0 = poly_mul(sigma0, term);
    }
    return sigma0;
}
Poly erase_positions(const Poly& received, const vector<int>& erasures){
    Poly out = received;
    for (int i : erasures) out[i] = 0;
    return out;
}
Poly compute_syndrome(const Poly& received){
    Poly S;
    for (int j = 1; j <= R; ++j)
        S.push_back( poly_eval(received, EXP_TABLE[j]) );
    return S;
}
Poly modified_syndrome(const Poly& S, const Poly& sigma0){
    Poly S0 = poly_mul(sigma0, S);
    S0.resize(R, 0);      // mod x^r  (keep <R terms)
    poly_trim(S0); return S0;
}

/* solves σ₁(x)S₀(x) = ω(x)  (mod x^r) */
pair<Poly,Poly> solve_key_equation(const Poly& S0, int e0){
    Poly r_poly(R+1, 0); r_poly.back() = 1;   // x^r
    int mu = (R - e0) / 2;
    int nu = (R + e0 - 1) / 2;
    auto [sigma1, omega] = extended_euclidean(r_poly, S0, mu, nu);

    /* normalize by σ₁(0) */
    int sigma1_0 = sigma1.empty()? 0 : sigma1[0];
    if (!sigma1_0) throw runtime_error("σ₁(0)=0, cannot normalize");
    for (int &c : sigma1) c = gf_div(c, sigma1_0);
    for (int &c : omega ) c = gf_div(c, sigma1_0);
    return {sigma1, omega};
}
Poly combine_locators(const Poly& a, const Poly& b){
    return poly_mul(a,b);
}
vector<int> find_error_positions(const Poly& sigma){
    vector<int> pos;
    for (int i = 0; i < N; ++i){
        int x_inv = EXP_TABLE[(63 - i) % 63];
        if (poly_eval(sigma, x_inv) == 0) pos.push_back(i);
    }
    return pos;
}
unordered_map<int,int> evaluate_error_magnitudes(
        const Poly& sigma, const Poly& omega,
        const vector<int>& positions)
{
    unordered_map<int,int> error_vec;
    Poly sigma_prime = poly_derivative(sigma);
    for (int i : positions){
        int x_inv = EXP_TABLE[(63 - i) % 63];
        int num = poly_eval(omega, x_inv);
        int den = poly_eval(sigma_prime, x_inv);
        if (!den) throw runtime_error("σ'(x_inv)=0 in Forney");
        error_vec[i] = gf_sub(0, gf_div(num, den));
    }
    return error_vec;
}
Poly apply_error_correction(Poly received,
        const unordered_map<int,int>& error_vec)
{
    for (auto [i,mag] : error_vec) received[i] = gf_sub(received[i], mag);
    return received;
}

/* ──────────────────────────────────────────────
   RS(63,42)  DECODER  (errors + erasures)
   Returns last K symbols (the decoded message)
   ────────────────────────────────────────────── */
// Poly rs_decode(const Poly& received, const vector<int>& erasures){
//     Poly sigma0   = build_erasure_locator(erasures);
//     Poly R_erased = erase_positions(received, erasures);
//     Poly S        = compute_syndrome(R_erased);
//     Poly S0       = modified_syndrome(S, sigma0);

//     auto [sigma1, omega] = solve_key_equation(S0, (int)erasures.size());
//     Poly sigma  = poly_mul(sigma0, sigma1);
//     auto errPos = find_error_positions(sigma);
//     auto errMag = evaluate_error_magnitudes(sigma, omega, errPos);
//     Poly corrected = apply_error_correction(received, errMag);
//     return Poly(corrected.end()-K, corrected.end());  // message (LSB-first)
// }
/* ─── Add near the top, after other constants ──────────────────────────────── */
static const Poly DECODE_FAIL = {-1};     // sentinel “decode failure”
/* ──────────────────────────────────────────────────────────────────────────── */

Poly rs_decode(const Poly& received, const vector<int>& erasures)
{
    try {
        /* 0. basic counts */
        const int t0 = static_cast<int>(erasures.size());

        /* 1. erasure preprocessing + syndromes */
        Poly sigma0   = build_erasure_locator(erasures);
        Poly R_erased = erase_positions(received, erasures);
        Poly S        = compute_syndrome(R_erased);
        Poly S0       = modified_syndrome(S, sigma0);

        /* 2. key-equation (σ₁, ω) */
        auto [sigma1, omega] = solve_key_equation(S0, t0);

        /* early budget test  t0 + 2·deg(σ₁) ≤ R ? */
        int t1_est = poly_deg(sigma1);
        // if (t0 + 2 * t1_est > R)          // outside radius
        //     return DECODE_FAIL;

        /* 3. locate errors, check root count */
        Poly sigma  = combine_locators(sigma0, sigma1);
        auto errPos = find_error_positions(sigma);

        // if (static_cast<int>(errPos.size()) != poly_deg(sigma))
        //     return DECODE_FAIL;           // locator/root mismatch ⇒ failure

        /* 4. magnitudes + correction */
        auto errMag    = evaluate_error_magnitudes(sigma, omega, errPos);
        Poly corrected = apply_error_correction(received, errMag);

        /* 5. late definitive test – all R syndromes must be zero */
        Poly S_chk = compute_syndrome(corrected);          // 21 syndromes
        for (int j = 0; j < R; ++j)
        if (S_chk[j]) {
            std::cerr << "FAIL  S" << std::setw(2) << j+1
                    << " = " << S_chk[j]
                    << "   (after correcting " << errPos.size()
                    << " errors + " << erasures.size() << " erasures)\n";
            break;                                     // one line per word is enough
        }


        /* 6. success – return decoded message (last K symbols) */
        return Poly(corrected.end() - K, corrected.end());
    }
    catch (const std::exception&) {
        /* any arithmetic / division-by-zero etc. ⇒ treat as decode failure */
        return DECODE_FAIL;
    }
}
/* ──────────────────────────────────────────────────────────────────────────── */



/*  DEMO-DRIVER  **************************************************************
 *  Each line in the input stream is one garbled codeword:
 *     45 63 * 0 2 … 15
 *  An asterisk means “erasure”; every other token is an integer 0-63.
 *  Blank lines (and lines that start with '#') are ignored.
 *  The function prints one line of output per input line:
 *    decoded symbols separated by spaces   – if decoding succeeds
 *    -1                                     – if t0 + 2t1 > 21  (or any error)
 ******************************************************************************/
void process_demo_input(std::istream& in, std::ostream& out = std::cout)
{
    std::string line;
    while (std::getline(in, line)) {
        // /* Skip empty / comment lines ---------------------------------------- */
        // std::string trimmed;
        // std::remove_copy_if(line.begin(), line.end(),
        //                     std::back_inserter(trimmed),
        //                     [](unsigned char c){ return std::isspace(c); });
        // if (trimmed.empty() || trimmed[0] == '#')  continue;

        /* Tokenise one codeword -------------------------------------------- */
        Poly received;
        std::vector<int> erasures;
        std::istringstream iss(line);
        std::string tok;
        int idx = 0;
        while (iss >> tok) {
            if (tok == "*") {                     // erasure
                received.push_back(0);            // placeholder value; ignored
                erasures.push_back(idx);
            } else {                              // regular symbol 0-63
                int val = std::stoi(tok);
                received.push_back(val);   
            }
            ++idx;
        }

        if (received.empty())                     // nothing useful on the line
            continue;

        /* Run decoder ------------------------------------------------------- */
        Poly decoded = rs_decode(received, erasures);

        /* Print result ------------------------------------------------------ */
        if (decoded == DECODE_FAIL) {
            out << -1 << '\n';
        } else {
            for (size_t i = 0; i < decoded.size(); ++i) {
                out << decoded[i] << (i + 1 == decoded.size() ? '\n' : ' ');
            }
        }
    }
}


/* ──────────────────────────────────────────────
   Simple main() with quick sanity tests.
   Comment out or adapt to assignment I/O as needed.
   ────────────────────────────────────────────── */
#ifdef RS_MAIN_DEMO
// int main() {
//     /* message 0,1,2,…,41 */
//     Poly msg(K); iota(msg.begin(), msg.end(), 0);

//     Poly cw = rs_encode(msg);
//     cout << "Encoded codeword ("<<cw.size()<<"):\n";
//     for (int x: cw) cout<<x<<' '; cout<<"\n\n";

//     /* test: 2 errors + 3 erasures  (t₀+2t₁ = 8 ≤ 21) */
//     Poly rx = cw;
//     rx[32]  = gf_add(rx[32], 12);        // error
//     rx[60] = gf_add(rx[60], 7);        // error
//     vector<int> ers = {2, 3, 45};      // erasures
//     for (int i: ers) rx[i] = 0;        // mark with 0 (erasures handled separately)
//     cout << "Rx codeword ("<<rx.size()<<"):\n";
//     for (int x: rx) cout<<x<<' '; cout<<"\n\n";

//     /* decode */
//     Poly dec_msg = rs_decode(rx, ers);

//     cout << "Decoded message:\n";
//     for (int x: dec_msg) cout<<x<<' '; cout<<"\n";
//     return 0;
// }
#endif
