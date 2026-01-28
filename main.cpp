#include <cmath>
#include <cctype>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

struct Particle {
    double m;
    double x, y, z;
    double vx, vy, vz;
    double fx, fy, fz;
};

static void usage(const char* prog) {
    std::cerr << "Usage: " << prog
              << " <N_or_inputfile> <dt> <steps> <dump_every> [output.tsv]\n"
              << "  If first arg is a number => random init with N particles.\n"
              << "  If first arg is a file   => load first state from TSV (solar.tsv format).\n";
}

static bool is_number(const std::string& s) {
    if (s.empty()) return false;
    for (char c : s) {
        if (!std::isdigit((unsigned char)c)) return false;
    }
    return true;
}


static void dump_state(std::ostream& os, const std::vector<Particle>& p) {
    os << p.size();
    os << std::setprecision(17);
    for (const auto& a : p) {
        os << '\t' << a.m
           << '\t' << a.x << '\t' << a.y << '\t' << a.z
           << '\t' << a.vx << '\t' << a.vy << '\t' << a.vz
           << '\t' << a.fx << '\t' << a.fy << '\t' << a.fz;
    }
    os << "\n";
}

static std::vector<std::string> split_tsv(const std::string& line) {
    std::vector<std::string> out;
    std::stringstream ss(line);
    std::string item;
    while (std::getline(ss, item, '\t')) out.push_back(item);
    return out;
}


static std::vector<Particle> load_first_state_tsv(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Could not open input file: " + path);

    std::string line;
    if (!std::getline(in, line)) throw std::runtime_error("Empty input file: " + path);

    auto tok = split_tsv(line);
    if (tok.empty()) throw std::runtime_error("Bad TSV line in: " + path);

    int nb = std::stoi(tok.at(0));
    if (nb <= 0) throw std::runtime_error("nbpart must be > 0");

    if ((int)tok.size() < 1 + nb * 10) {
        throw std::runtime_error("Not enough columns in TSV for nbpart=" + std::to_string(nb));
    }

    std::vector<Particle> p(nb);
    int idx = 1;
    for (int i = 0; i < nb; i++) {
        p[i].m  = std::stod(tok.at(idx + 0));
        p[i].x  = std::stod(tok.at(idx + 1));
        p[i].y  = std::stod(tok.at(idx + 2));
        p[i].z  = std::stod(tok.at(idx + 3));
        p[i].vx = std::stod(tok.at(idx + 4));
        p[i].vy = std::stod(tok.at(idx + 5));
        p[i].vz = std::stod(tok.at(idx + 6));
        p[i].fx = std::stod(tok.at(idx + 7));
        p[i].fy = std::stod(tok.at(idx + 8));
        p[i].fz = std::stod(tok.at(idx + 9));
        idx += 10;
    }
    return p;
}

static std::vector<Particle> init_random(int n) {
    std::mt19937_64 rng(1234567);

    std::uniform_real_distribution<double> mass(1.0, 10.0);
    std::uniform_real_distribution<double> pos(-1.0, 1.0);
    std::uniform_real_distribution<double> vel(-0.01, 0.01);

    std::vector<Particle> p(n);
    for (int i = 0; i < n; i++) {
        p[i].m = mass(rng);
        p[i].x = pos(rng); p[i].y = pos(rng); p[i].z = pos(rng);
        p[i].vx = vel(rng); p[i].vy = vel(rng); p[i].vz = vel(rng);
        p[i].fx = p[i].fy = p[i].fz = 0.0;
    }
    return p;
}

static void reset_forces(std::vector<Particle>& p) {
    for (auto& a : p) a.fx = a.fy = a.fz = 0.0;
}


// F = G*m1*m2/r^2 
// F_vec = G*m1*m2 * (r_vec) / |r|^3
static void compute_forces(std::vector<Particle>& p) {
    const double G = 6.674e-11;
    const double eps2 = 1e-9; 

    reset_forces(p);

    const int n = (int)p.size();
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double dx = p[j].x - p[i].x;
            double dy = p[j].y - p[i].y;
            double dz = p[j].z - p[i].z;

            double r2 = dx*dx + dy*dy + dz*dz + eps2;
            double inv_r = 1.0 / std::sqrt(r2);
            double inv_r3 = inv_r * inv_r * inv_r;

            double f = G * p[i].m * p[j].m * inv_r3;
            double fx = f * dx;
            double fy = f * dy;
            double fz = f * dz;

            // i pulled toward j, j pulled toward i
            p[i].fx += fx; p[i].fy += fy; p[i].fz += fz;
            p[j].fx -= fx; p[j].fy -= fy; p[j].fz -= fz;
        }
    }
}


static void integrate(std::vector<Particle>& p, double dt) {
    for (auto& a : p) {
        double ax = a.fx / a.m;
        double ay = a.fy / a.m;
        double az = a.fz / a.m;

        a.vx += ax * dt;
        a.vy += ay * dt;
        a.vz += az * dt;

        a.x += a.vx * dt;
        a.y += a.vy * dt;
        a.z += a.vz * dt;
    }
}

int main(int argc, char** argv) {
    if (argc < 5 || argc > 6) {
        usage(argv[0]);
        return 1;
    }

    std::string src = argv[1];
    double dt = 0.0;
    long long steps = 0;
    long long dump_every = 0;

    try {
        dt = std::stod(argv[2]);
        steps = std::stoll(argv[3]);
        dump_every = std::stoll(argv[4]);
    } catch (...) {
        usage(argv[0]);
        return 1;
    }

    if (dt <= 0 || steps < 0 || dump_every <= 0) {
        std::cerr << "Error: dt must be > 0, steps >= 0, dump_every > 0\n";
        return 1;
    }

    std::vector<Particle> p;
    try {
        if (is_number(src)) {
            int n = std::stoi(src);
            if (n <= 0) { std::cerr << "N must be > 0\n"; return 1; }
            p = init_random(n);
        } else {
            p = load_first_state_tsv(src);
        }
    } catch (const std::exception& e) {
        std::cerr << "Init error: " << e.what() << "\n";
        return 1;
    }

    std::ofstream ofs;
    std::ostream* osp = &std::cout;
    if (argc == 6) {
        ofs.open(argv[5]);
        if (!ofs) {
            std::cerr << "Could not open output file: " << argv[5] << "\n";
            return 1;
        }
        osp = &ofs;
    }

    // Time 
    auto t0 = std::chrono::high_resolution_clock::now();

    for (long long t = 0; t < steps; t++) {
        compute_forces(p);
        integrate(p, dt);

        if (t % dump_every == 0) {
            dump_state(*osp, p);
        }
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> sec = t1 - t0;

    // Print timing to stderr
    std::cerr << "SIM_TIME_SECONDS " << sec.count() << "\n";

    return 0;
}
