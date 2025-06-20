#define DBG(x) do { if (debug) { x; } } while (0)
static bool debug = true;      // flip to false when you’re satisfied

#include <iostream>
#include <vector>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <map>
#include <algorithm>
#include <cmath> 


using namespace std;
const int N = 63;
const int K = 42;
const int R = N - K;

struct RSDecodeError : public std::exception {
    std::string message;
    RSDecodeError(const std::string& msg) : message(msg) {}
    const char* what() const noexcept override { return message.c_str(); }
};

// EXP_TABLE[i] = α^i
const int EXP_TABLE[63] = {
    1, 2, 4, 8, 16, 32, 3, 6, 12, 24,
    48, 35, 5, 10, 20, 40, 19, 38, 15, 30,
    60, 59, 53, 41, 17, 34, 7, 14, 28, 56,
    51, 37, 9, 18, 36, 11, 22, 44, 27, 54,
    47, 29, 58, 55, 45, 25, 50, 39, 13, 26,
    52, 43, 21, 42, 23, 46, 31, 62, 63, 61,
    57, 49, 33
};

int LOG_TABLE[64];

const vector<int> GENERATOR_POLY_COEFFS = {
    58, 62, 59, 7, 35, 58, 63, 47, 51, 6, 33,
    43, 44, 27, 7, 53, 39, 62, 52, 41, 44, 1
};

const int G_EVAL_LIST[63] = {
    34, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 16, 13, 43, 41, 48, 15, 11, 52, 33,
    1, 22, 20, 28, 17, 13, 19, 14, 56, 25,
    45, 10, 46, 45, 8, 12, 14, 44, 14, 17,
    26, 31, 22, 59, 29, 52, 31, 57, 48, 45,
    51, 13
};

// Initialization of LOG_TABLE
void init_log_table() {
    for (int i = 0; i < 64; ++i) LOG_TABLE[i] = -1;
    for (int i = 0; i < 63; ++i)
        LOG_TABLE[EXP_TABLE[i]] = i;
}

// Field operations
inline int gf_add(int a, int b) {
    return a ^ b;
}

inline int gf_sub(int a, int b) {
    return a ^ b;
}

int gf_mul(int a, int b) {
    if (a == 0 || b == 0) return 0;
    return EXP_TABLE[(LOG_TABLE[a] + LOG_TABLE[b]) % 63];
}

int gf_div(int a, int b) {
    if (b == 0) throw runtime_error("Division by zero in GF(2^6)");
    if (a == 0) return 0;
    return EXP_TABLE[(LOG_TABLE[a] - LOG_TABLE[b] + 63) % 63];
}

int gf_inv(int a) {
    if (a == 0) throw runtime_error("Inverse of 0 in GF(2^6)");
    return EXP_TABLE[(63 - LOG_TABLE[a]) % 63];
}

int gf_pow(int a, int n) {
    if (a == 0) return 0;
    return EXP_TABLE[(LOG_TABLE[a] * n) % 63];
}

// Polynomial evaluation at x
int poly_eval(const vector<int>& poly_coeffs, int x) {
    int result = 0, power = 1;
    for (int coeff : poly_coeffs) {
        result = gf_add(gf_mul(coeff, power), result);
        power = gf_mul(power, x);
    }
    return result;
}

// Verification routine
void verify_generator_polynomial() {
    bool passed = true;
    for (int i = 0; i < 63; ++i) {
        int alpha_i = EXP_TABLE[i];
        int expected = G_EVAL_LIST[i];
        int actual = poly_eval(GENERATOR_POLY_COEFFS, alpha_i);
        if (actual != expected) {
            cout << "Mismatch at alpha^" << i << ": expected " << expected << ", got " << actual << endl;
            passed = false;
        }
    }
    if (passed)
        cout << "All g(alpha^i) values match the official table." << endl;
}

// polynomial helper functions

vector<int> poly_add(const vector<int>& f, const vector<int>& g) {
    size_t max_len = max(f.size(), g.size());
    vector<int> f_pad = f, g_pad = g;
    f_pad.resize(max_len, 0);
    g_pad.resize(max_len, 0);
    vector<int> result(max_len);
    for (size_t i = 0; i < max_len; ++i)
        result[i] = gf_add(f_pad[i], g_pad[i]);
    return result;
}

vector<int> poly_mul(const vector<int>& f, const vector<int>& g) {
    vector<int> result(f.size() + g.size() - 1, 0);
    for (size_t i = 0; i < f.size(); ++i)
        for (size_t j = 0; j < g.size(); ++j)
            result[i + j] = gf_add(result[i + j], gf_mul(f[i], g[j]));
    return result;
}

vector<int> poly_scale(const vector<int>& p, int scalar) {
    vector<int> result;
    for (int c : p)
        result.push_back(gf_mul(c, scalar));
    return result;
}

vector<int> poly_shift(const vector<int>& p, int n) {
    vector<int> result(n, 0);
    result.insert(result.end(), p.begin(), p.end());
    return result;
}

int poly_deg(const vector<int>& p) {
    for (int i = static_cast<int>(p.size()) - 1; i >= 0; --i)
        if (p[i] != 0) return i;
    return -1;
}

vector<int> poly_trim(const vector<int>& p) {
    int i = static_cast<int>(p.size()) - 1;
    while (i >= 0 && p[i] == 0) --i;
    return vector<int>(p.begin(), p.begin() + i + 1);
}

vector<int> poly_make_monic(const vector<int>& p) {
    vector<int> p_trimmed = poly_trim(p);
    if (p_trimmed.empty()) return {};
    int inv = gf_inv(p_trimmed.back());
    vector<int> result;
    for (int c : p_trimmed)
        result.push_back(gf_mul(c, inv));
    return result;
}

pair<vector<int>, vector<int>> poly_divmod(vector<int> f, vector<int> g) {
    f = poly_trim(f);
    g = poly_trim(g);
    int deg_f = poly_deg(f);
    int deg_g = poly_deg(g);

    if (deg_g < 0)
        throw runtime_error("Division by zero polynomial");

    vector<int> quotient(deg_f - deg_g + 1, 0);
    vector<int> remainder = f;

    while (poly_deg(remainder) >= deg_g) {
        int shift = poly_deg(remainder) - deg_g;
        int lead_coeff = gf_div(remainder.back(), g.back());

        vector<int> scaled_g = poly_scale(g, lead_coeff);
        vector<int> aligned_g(shift, 0);
        aligned_g.insert(aligned_g.end(), scaled_g.begin(), scaled_g.end());

        quotient[shift] = lead_coeff;
        remainder = poly_trim(poly_add(remainder, aligned_g));
    }

    return {quotient, remainder};
}

vector<int> poly_gcd(vector<int> a, vector<int> b) {
    a = poly_trim(a);
    b = poly_trim(b);
    if (b.empty()) return poly_make_monic(a);
    while (!b.empty()) {
        auto [q, r] = poly_divmod(a, b);
        a = b;
        b = poly_trim(r);
    }
    return poly_make_monic(a);
}

pair<vector<int>, vector<int>> extended_euclidean(vector<int> a, vector<int> b, int mu, int nu) {
    a = poly_trim(a);
    b = poly_trim(b);

    vector<int> r_prev = a, r_curr = b;
    vector<int> u_prev = {1}, u_curr = {0};
    vector<int> v_prev = {0}, v_curr = {1};

    while (poly_deg(r_curr) > nu || poly_deg(v_curr) > mu) {
        auto [q, r_next] = poly_divmod(r_prev, r_curr);
        r_prev = r_curr;
        r_curr = r_next;

        vector<int> u_next = poly_trim(poly_add(u_prev, poly_mul(q, u_curr)));
        u_prev = u_curr;
        u_curr = u_next;

        vector<int> v_next = poly_trim(poly_add(v_prev, poly_mul(q, v_curr)));
        v_prev = v_curr;
        v_curr = v_next;
    }

    return {v_curr, r_curr}; // (σ1, ω)
}

vector<int> build_erasure_locator(const vector<int>& erasures) {
    vector<int> sigma0 = {1};
    for (int i : erasures) {
        vector<int> term = {1, gf_sub(0, EXP_TABLE[i])};
        // vector<int> term = {gf_sub(0, EXP_TABLE[i]), 1};
        sigma0 = poly_mul(sigma0, term);
    }
    return sigma0;
}

vector<int> erase_positions(vector<int> received, const vector<int>& erasures) {
    for (int i : erasures)
        received[i] = 0;
    return received;
}

vector<int> compute_syndrome(const vector<int>& received) {
    vector<int> S;
    for (int j = 1; j <= R; ++j)
        S.push_back(poly_eval(received, EXP_TABLE[j]));
    return S;
}

vector<int> modified_syndrome(const vector<int>& syndrome, const vector<int>& sigma0) {
    vector<int> S0 = poly_mul(sigma0, syndrome);
    S0.resize(R);
    return S0;
}

pair<vector<int>, vector<int>> solve_key_equation(const vector<int>& s0, int e0) {
    vector<int> r_poly(R + 1, 0);
    r_poly[R] = 1;
    // int mu = (R - e0) / 2;
    // int nu = (R + e0 - 1) / 2;

    int mu = static_cast<int>(floor((R - e0) / 2.0));
    int nu = static_cast<int>(ceil((R + e0) / 2.0)) - 1;

    auto [sigma1, omega] = extended_euclidean(r_poly, s0, mu, nu);

    // Optional normalization
    int sigma1_0 = sigma1[0];
    if (sigma1_0 == 0)
        throw RSDecodeError("sigma1(0) = 0, cannot normalize");

    for (int& c : sigma1) c = gf_div(c, sigma1_0);
    for (int& c : omega) c = gf_div(c, sigma1_0);

    return {sigma1, omega};
}
bool solve_linear_system(vector<vector<int>> A, vector<int> b, vector<int>& sol) {
    int m = A.size();
    int n = A[0].size();
    sol.assign(n, 0);

    for (int col = 0, row = 0; col < n && row < m; ++col) {
        int sel = -1;
        for (int i = row; i < m; ++i) {
            if (A[i][col] != 0) {
                sel = i;
                break;
            }
        }
        if (sel == -1) continue;

        swap(A[sel], A[row]);
        swap(b[sel], b[row]);

        int inv = gf_inv(A[row][col]);
        for (int j = col; j < n; ++j)
            A[row][j] = gf_mul(A[row][j], inv);
        b[row] = gf_mul(b[row], inv);

        for (int i = 0; i < m; ++i) {
            if (i != row && A[i][col] != 0) {
                int f = A[i][col];
                for (int j = col; j < n; ++j)
                    A[i][j] = gf_sub(A[i][j], gf_mul(f, A[row][j]));
                b[i] = gf_sub(b[i], gf_mul(f, b[row]));
            }
        }
        ++row;
    }

    // Back-substitution: now A should be reduced to row-echelon form
    for (int i = 0, j = 0; i < m && j < n; ++j) {
        if (A[i][j] == 1) {
            sol[j] = b[i];
            ++i;
        }
    }

    return true; 
}
pair<vector<int>, vector<int>> solve_key_wb(const vector<int>& s0, int e0) {
    int mu = (R - e0) / 2;
    int nu = (R + e0 + 1) / 2 - 1;

    int n_eqs = R;
    int n_unknowns = (mu + 1) + (nu + 1); // coefficients of sigma1 and omega

    vector<vector<int>> A(n_eqs, vector<int>(n_unknowns, 0));
    vector<int> b(n_eqs, 0);

    for (int j = 0; j < R; ++j) {
        int x = EXP_TABLE[j + 1];
        int Sj = poly_eval(s0, x);  // S0(α^j)
        b[j] = Sj;

        int power = 1;
        for (int i = 0; i <= mu; ++i) {
            A[j][i] = gf_mul(Sj, power); // S0 * x^i
            power = gf_mul(power, x);
        }

        power = 1;
        for (int i = 0; i <= nu; ++i) {
            A[j][mu + 1 + i] = gf_sub(0, power); // -x^i for ω(x)
            power = gf_mul(power, x);
        }
    }

    // Solve the system A * [sigma1_coeffs | omega_coeffs]^T = b
    vector<int> sol;
    bool ok = solve_linear_system(A, b, sol); // You’ll need to implement this

    if (!ok) throw RSDecodeError("WB solve failed");

    vector<int> sigma1(sol.begin(), sol.begin() + mu + 1);
    vector<int> omega(sol.begin() + mu + 1, sol.end());

    if (sigma1[0] == 0)
        throw RSDecodeError("WB failure: sigma1(0) = 0");

    int inv = gf_inv(sigma1[0]);
    for (int& c : sigma1) c = gf_mul(c, inv);
    for (int& c : omega) c = gf_mul(c, inv);

    return {sigma1, omega};
}


vector<int> combine_locators(const vector<int>& sigma0, const vector<int>& sigma1) {
    return poly_mul(sigma0, sigma1);
}

vector<int> find_error_positions(const vector<int>& sigma) {
    vector<int> positions;
    for (int i = 0; i < N; ++i) {
        int x_inv = EXP_TABLE[(63 - i) % 63];
        if (poly_eval(sigma, x_inv) == 0)
            positions.push_back(i);
    }
    return positions;
}

vector<int> poly_derivative(const vector<int>& p) {
    if (p.size() < 2) return {0};
    vector<int> deriv(p.size() - 1, 0);
    for (size_t i = 1; i < p.size(); ++i)
        if (i & 1) deriv[i - 1] = p[i];
    while (deriv.size() > 1 && deriv.back() == 0)
        deriv.pop_back();
    return deriv;
}

map<int, int> evaluate_error_magnitudes(const vector<int>& sigma, const vector<int>& omega, const vector<int>& positions) {
    map<int, int> error_vector;
    vector<int> sigma_prime = poly_derivative(sigma);
    for (int i : positions) {
        int x_inv = EXP_TABLE[(63 - i) % 63];
        int num = poly_eval(omega, x_inv);
        int denom = poly_eval(sigma_prime, x_inv);
        if (denom == 0)
            throw RSDecodeError("σ'(α^-" + to_string(i) + ") = 0 during Forney eval");
        error_vector[i] = gf_sub(0, gf_div(num, denom));
    }
    return error_vector;
}

vector<int> apply_error_correction(vector<int> received, const map<int, int>& error_vector) {
    for (const auto& [i, mag] : error_vector)
        received[i] = gf_sub(received[i], mag);
    return received;
}
auto print_poly = [](const string& name, const vector<int>& p) {
    cout << "  " << name << " = [";
    for (size_t i = 0; i < p.size(); ++i) {
        cout << p[i];
        if (i < p.size() - 1) cout << " ";
    }
    cout << "] (deg=" << poly_deg(p) << ")\n";
};



vector<int> rs_decode(const vector<int>& received, const vector<int>& erasures) {
    vector<int> sigma0 = build_erasure_locator(erasures);
    vector<int> R_erased = erase_positions(received, erasures);
    vector<int> S = compute_syndrome(R_erased);
    vector<int> S0 = modified_syndrome(S, sigma0);
    if (poly_deg(S0)==-1) return vector<int>(received.end()-K, received.end());

    // auto [sigma1, omega] = solve_key_wb(S0, erasures.size());
    auto [sigma1, omega] = solve_key_equation(S0, erasures.size());



    int t0 = erasures.size();
    cout << "  t0=" << t0 
     << "  deg(σ1)=" << poly_deg(sigma1)
     << "  deg(ω)=" << poly_deg(omega)
     << "  deg(S)=" << poly_deg(S)
     << "  deg(S0)=" << poly_deg(S0)
     << "  deg(σ0)=" << poly_deg(sigma0) << endl;
    print_poly("σ1", sigma1);
    print_poly("ω", omega);
    print_poly("S", S);
    print_poly("S0", S0);
    print_poly("σ0", sigma0);

    if (poly_deg(omega) >= erasures.size() + poly_deg(sigma1)) {
        //  << "deg(sigma1): " << poly_deg(sigma1) << ", deg(omega): " 
        // << poly_deg(omega) << ", deg(sigma0) = " << poly_deg(sigma0) << endl;

        throw RSDecodeError("Condition (A) violated: deg(ω)≥ t₀ + deg(σ₁)");
    }
    int t1_budget = (R - t0) / 2;
    if (poly_deg(sigma1) > t1_budget) {
        throw RSDecodeError("Decoding failure: t0=" + to_string(t0)
            + ", 2·deg(σ₁)=" + to_string(2 * poly_deg(sigma1))
            + ", R=" + to_string(R));
    }
    if (t0 > R) {
        throw RSDecodeError("Radius exceeded: t₀ = " +
            to_string(t0) + " > R = " + to_string(R));
    }

    int t1_est = poly_deg(sigma1);
    if (t0 + 2 * t1_est > R) {
        throw RSDecodeError("Radius exceeded: t₀ + 2·t₁ = " +
            to_string(t0 + 2 * t1_est) + " > R = " + to_string(R));
    }
    // if (t0 == R && t1_est == 0 && poly_deg(omega) == R - 1)
    //     throw RSDecodeError("σ1 trivial while t0 == R : uncorrectable");

    
    vector<int> sigma = combine_locators(sigma0, sigma1);
    vector<int> error_pos = find_error_positions(sigma);
    if ((int)error_pos.size() != poly_deg(sigma)) {
        throw RSDecodeError("Locator degree ≠ number of roots found (beyond radius)");
    }

    auto error_mag = evaluate_error_magnitudes(sigma, omega, error_pos);
    vector<int> corrected = apply_error_correction(received, error_mag);

    if (any_of(compute_syndrome(corrected).begin(), compute_syndrome(corrected).end(), [](int s){ return s != 0; })) {
        throw RSDecodeError("Syndrome non-zero after correction (beyond radius)");
    }

    return vector<int>(corrected.end() - K, corrected.end());  // last K symbols
}

vector<vector<int>> interpolate_gs_poly(const vector<int>& xs, const vector<int>& ys, int dA, int dB) {
    int m = xs.size();
    int num_unknowns = dA + dB + 2;

    vector<vector<int>> A(m, vector<int>(num_unknowns, 0));
    vector<int> b(m, 0);

    for (int i = 0; i < m; ++i) {
        int x = xs[i];
        int y = ys[i];

        int xp = 1;
        for (int j = 0; j <= dA; ++j) {
            A[i][j] = xp;
            xp = gf_mul(xp, x);
        }

        int yp = y;
        xp = 1;
        for (int j = 0; j <= dB; ++j) {
            A[i][dA + 1 + j] = gf_mul(yp, xp);
            xp = gf_mul(xp, x);
        }

        b[i] = 0;
    }

    vector<int> sol;
    bool ok = solve_linear_system(A, b, sol);
    if (!ok) throw RSDecodeError("GS interpolation failed");

    // Split solution into A(x) and B(x)
    vector<int> Acoeffs(sol.begin(), sol.begin() + dA + 1);
    vector<int> Bcoeffs(sol.begin() + dA + 1, sol.end());
    return {Acoeffs, Bcoeffs}; // Q(x, y) = A(x) + B(x) * y
}
vector<int> gs_decode(const vector<int>& received) {
    vector<int> xs, ys;
    for (int i = 0; i < N; ++i) {
        xs.push_back(EXP_TABLE[i + 1]);
        ys.push_back(received[i]);
    }

    vector<vector<int>> Q = interpolate_gs_poly(xs, ys, K - 1, K - 2);
    vector<int> A = Q[0];
    vector<int> B = Q[1];

    if (poly_deg(B) < 0) throw RSDecodeError("Degenerate Q(x,y) in GS decode");

    vector<int> f = poly_divmod(poly_trim(poly_scale(A, gf_sub(0, 1))), poly_trim(B)).first;

    if (poly_deg(f) >= K) throw RSDecodeError("Decoded polynomial exceeds degree");

    return f;
}

int main() {
    init_log_table();
    verify_generator_polynomial();

    ifstream infile("input.txt");
    ofstream outfile("output.txt");

    if (!infile.is_open()) {
        cerr << "❌ Cannot open input.txt\n";
        return 1;
    }
    if (!outfile.is_open()) {
        cerr << "❌ Cannot open output.txt\n";
        return 1;
    }

    string line;
    int line_num = 0;
    while (getline(infile, line)) {
        ++line_num;
        vector<int> received;
        vector<int> erasures;
        istringstream iss(line);
        string token;

        int pos = 0;
        while (iss >> token) {
            if (token == "*") {
                received.push_back(0); // use 0 for erased symbol
                erasures.push_back(pos);
            } else {
                received.push_back(stoi(token));
            }
            ++pos;
        }

        try {
            // vector<int> decoded = gs_decode(received);
            vector<int> decoded = rs_decode(received, erasures);
            // for (int v : decoded) outfile << v << " ";
            for (int i = 0; i < K; ++i) {
                outfile << decoded[i]; if (i < K - 1) outfile << " ";
            }

            outfile << "\n";
        } catch (const RSDecodeError& e) {
            cerr << "❌ Line " << line_num << ": " << e.what() << "\n";
            for (int i = 0; i < K; ++i) outfile << "* ";
            outfile << "\n";
        } catch (const exception& e) {
            cerr << "❌ Line " << line_num << ": Unexpected error: " << e.what() << "\n";
            for (int i = 0; i < K; ++i) outfile << "* ";
            outfile << "\n";
        }
    }

    infile.close();
    outfile.close();

    cout << "✅ Decoding complete. Results saved to output.txt\n";
    return 0;
}
