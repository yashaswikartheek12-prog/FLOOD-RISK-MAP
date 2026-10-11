// Flood risk analysis (simplified C++ version)
// Same core idea as the CruzHacks project: elevation + rainfall + land use -> flood risk map.
// Build: g++ -O2 -std=c++17 flood_risk.cpp -o flood_risk
// Run (synthetic demo): ./flood_risk
// Run (your data):      ./flood_risk data.csv WIDTH HEIGHT
//   data.csv rows (row-major, one per cell): elevation_m,rainfall_mm,land_class
//   land_class: 0=Water 1=Building 2=Vegetation 3=Crop 4=Bare

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

enum class Land { Water = 0, Building, Vegetation, Crop, Bare };

// How much rain becomes surface runoff for each land type
double runoffCoeff(Land l) {
    switch (l) {
        case Land::Water:      return 1.00;
        case Land::Building:   return 0.90;
        case Land::Bare:       return 0.60;
        case Land::Crop:       return 0.40;
        case Land::Vegetation: return 0.25;
    }
    return 0.5;
}

struct Cell {
    double elev = 0, rain = 0;
    Land land = Land::Bare;
    double flow = 0;   // accumulated runoff reaching this cell
    double risk = 0;   // final score 0..1
    bool sink = false; // no lower neighbour -> water pools here
};

class FloodModel {
public:
    int W, H;
    vector<Cell> g;

    FloodModel(int w, int h) : W(w), H(h), g(w * h) {}
    Cell& at(int x, int y) { return g[y * W + x]; }

    void run() {
        accumulateFlow();
        scoreRisk();
    }

private:
    // D8 flow accumulation: each cell sends its water to its steepest lower neighbour
    void accumulateFlow() {
        for (auto& c : g) { c.flow = c.rain * runoffCoeff(c.land); c.sink = false; }

        vector<int> order(g.size());
        iota(order.begin(), order.end(), 0);
        sort(order.begin(), order.end(), [&](int a, int b) { return g[a].elev > g[b].elev; });

        for (int idx : order) {
            int x = idx % W, y = idx / W;
            double bestSlope = 0;
            int best = -1;
            for (int dy = -1; dy <= 1; dy++)
                for (int dx = -1; dx <= 1; dx++) {
                    if (!dx && !dy) continue;
                    int nx = x + dx, ny = y + dy;
                    if (nx < 0 || ny < 0 || nx >= W || ny >= H) continue;
                    double dist = (dx && dy) ? 1.4142 : 1.0;
                    double slope = (g[idx].elev - at(nx, ny).elev) / dist;
                    if (slope > bestSlope) { bestSlope = slope; best = ny * W + nx; }
                }
            if (best >= 0) g[best].flow += g[idx].flow;
            else g[idx].sink = true;
        }
    }

    // Risk = low-lying + lots of water flowing in + surface that doesn't absorb
    void scoreRisk() {
        double minE = 1e18, maxE = -1e18, maxF = 0;
        for (auto& c : g) {
            minE = min(minE, c.elev); maxE = max(maxE, c.elev);
            maxF = max(maxF, c.flow);
        }
        for (auto& c : g) {
            double lowness  = 1.0 - (c.elev - minE) / (maxE - minE + 1e-9);
            double flowNorm = log1p(c.flow) / log1p(maxF + 1e-9);
            double surface  = runoffCoeff(c.land);
            double r = 0.35 * lowness + 0.45 * flowNorm + 0.20 * surface;
            if (c.sink) r += 0.10;  // depressions collect water
            c.risk = min(1.0, r);
        }
    }
};

// ---------- I/O helpers ----------
bool loadCSV(const string& path, FloodModel& m) {
    ifstream f(path);
    if (!f) return false;
    string line;
    size_t i = 0;
    while (getline(f, line) && i < m.g.size()) {
        if (line.empty() || !(isdigit(line[0]) || line[0] == '-' || line[0] == '.')) continue; // skip header
        stringstream ss(line);
        string a, b, c;
        getline(ss, a, ','); getline(ss, b, ','); getline(ss, c, ',');
        m.g[i].elev = stod(a);
        m.g[i].rain = stod(b);
        int cls = stoi(c);
        m.g[i].land = static_cast<Land>(clamp(cls, 0, 4));
        i++;
    }
    return i == m.g.size();
}

void makeSynthetic(FloodModel& m) {
    mt19937 rng(42);
    normal_distribution<double> noise(0, 0.6);
    uniform_real_distribution<double> u(0, 1);
    for (int y = 0; y < m.H; y++)
        for (int x = 0; x < m.W; x++) {
            Cell& c = m.at(x, y);
            // a valley running through the middle with hills on the sides
            double valley = fabs(x - m.W / 2.0) * 0.8;
            c.elev = 10 + valley + 2.0 * sin(y / 5.0) + noise(rng);
            c.rain = 120 + 30 * u(rng);  // mm in a storm event
            if (c.elev < 11)          c.land = Land::Water;
            else if (c.elev < 16)     c.land = (u(rng) < 0.5) ? Land::Building : Land::Crop;
            else if (c.elev < 22)     c.land = (u(rng) < 0.6) ? Land::Crop : Land::Bare;
            else                      c.land = Land::Vegetation;
        }
}

void printMap(const FloodModel& m) {
    cout << "Legend: ~ water   . low risk   o medium   # HIGH\n\n";
    for (int y = 0; y < m.H; y++) {
        for (int x = 0; x < m.W; x++) {
            const Cell& c = m.g[y * m.W + x];
            char ch = c.land == Land::Water ? '~'
                    : c.risk >= 0.66 ? '#'
                    : c.risk >= 0.40 ? 'o' : '.';
            cout << ch;
        }
        cout << '\n';
    }
}

void writeCSV(const FloodModel& m, const string& path) {
    ofstream f(path);
    f << "x,y,elevation,land_class,flow,risk,level\n";
    for (int y = 0; y < m.H; y++)
        for (int x = 0; x < m.W; x++) {
            const Cell& c = m.g[y * m.W + x];
            const char* lvl = c.risk >= 0.66 ? "HIGH" : c.risk >= 0.40 ? "MEDIUM" : "LOW";
            f << x << ',' << y << ',' << c.elev << ',' << static_cast<int>(c.land)
              << ',' << c.flow << ',' << c.risk << ',' << lvl << '\n';
        }
}

int main(int argc, char** argv) {
    int W = 60, H = 25;
    FloodModel* m;
    if (argc >= 4) {
        W = stoi(argv[2]); H = stoi(argv[3]);
        m = new FloodModel(W, H);
        if (!loadCSV(argv[1], *m)) { cerr << "Bad CSV or size mismatch\n"; return 1; }
    } else {
        m = new FloodModel(W, H);
        makeSynthetic(*m);
        cout << "(no input given, using synthetic terrain)\n\n";
    }

    m->run();
    printMap(*m);
    writeCSV(*m, "risk_output.csv");

    int hi = 0, med = 0, lo = 0;
    for (auto& c : m->g) { if (c.risk >= 0.66) hi++; else if (c.risk >= 0.40) med++; else lo++; }
    cout << "\nHIGH: " << hi << "  MEDIUM: " << med << "  LOW: " << lo
         << "   (full grid saved to risk_output.csv)\n";
    delete m;
    return 0;
}
