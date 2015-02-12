#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <sstream>
#include <getopt.h>
using namespace std;

int N = 0;
int M = 0;
double E1 = -1;
double E2 = -1;
double E3 = -1;
double W = -1;
double P = -1;
bool isNormal = false;
double lambda = -1;
double mu = -1;
double sigma = -1;
int maxRejected = 0;

double normal_pmf(double x, double m, double s) {
    static const double inv_sqrt_2pi = 0.3989422804014327;
    double a = (x - m) / s;
    return inv_sqrt_2pi / s * exp(-0.5f * a * a);
}

double poisson_pmf(const double k, const double l) {
  return exp(k * log(l) - lgamma(k + 1.0) - l);
}

string getHash(const int o, const int oi, const int c, const int ci, const int w, const int r) {
  stringstream ss;
  ss << o << "|" << oi << "|" << c << "|" << ci << "|" << w << "|" << r;
  return ss.str();
}

struct State {
  State() {}
  State(int o, int oi, int c, int ci, int w, int r) : o(o), oi(oi), c(c), ci(ci), w(w), r(r) {
    hash = getHash(o, oi, c, ci, w, r);
    u = 0;
  }
  bool operator<(const State& rhs) const {
    return hash < rhs.hash;
  }
  bool operator==(const State& rhs) const {
    return o == rhs.o &&
           oi == rhs.oi &&
           c == rhs.c &&
           ci == rhs.ci &&
           w == rhs.w &&
           r == rhs.r &&
           (abs(u-rhs.u) < 10e-8);
  }
  double getReward() const {
    return -1 * o*E3 + oi*E1 + ci*E2 + w*W * r*P;
  }

  int o; // open
  int oi; // opening
  int c; // closed
  int ci; // closing
  int w; // waiting
  int r; // rejected
  double u; // utility
  string hash;
};

typedef unordered_map<string, State> StateMap_t;
StateMap_t stateMap;

void populateStates() {
  for (int numOpen = 0; numOpen <= M; ++numOpen) {
    for (int numOpening = 0; numOpening <= M; ++numOpening) {
      for (int numClosed = 0; numClosed <= M; ++numClosed) {
        for (int numClosing = 0; numClosing <= M; ++numClosing) {
          for (int numWaiting = 0; numWaiting <= N; ++numWaiting) {
            for (int numRejected = 0; numRejected <= maxRejected; ++numRejected) {
              if (numOpen + numOpening + numClosed + numClosing > M) continue;
              State s(numOpen, numOpening, numClosed, numClosing, numWaiting, numRejected);
              stateMap[s.hash] = s; 
            }
          }
        }
      }
    }
  }
}

double getBestExpectedUtility(const State& s) {
  double bestUtility = -std::numeric_limits<double>::max();
  for (int closed = 0; closed <= M; ++closed) {
    for (int opened = 0; opened <= M; ++opened) {
      if (closed < s.o || opened > s.c) continue;
      int numAvailRooms = s.o + s.oi - s.w;
      int numLeftoverPatients = 0;
      if (numAvailRooms < 0) {
        numLeftoverPatients = abs(numAvailRooms);
        numAvailRooms = 0;
      }
      const int numAvailSeats = N - numLeftoverPatients;
      for (int i = 0; i <= numAvailSeats+maxRejected; ++i) {
        double p = 0;
        if (isNormal) {
          p = normal_pmf(i, mu, sigma);
        } else {
          p = poisson_pmf(i, lambda);
        }
        if (p < 0.001) continue;

        const int r = (i > numAvailSeats ? i-numAvailSeats : i);
        const int w = (r > 0 ? numAvailSeats : i);
        State s1(numAvailRooms, opened, s.c+s.ci, closed, w, r);
        if (stateMap.find(s1.hash) == stateMap.end()) {
          continue;
        }
        bestUtility = max(bestUtility, p*s1.getReward());
      }
    }
  }
  return bestUtility;
}

void valueIteration() {
  int nIter = 0;
  while (true) {
    StateMap_t tmp;
    int i = 0;
    for (const auto& p : stateMap) {
      State s = p.second; // copy
      s.u = s.getReward() + getBestExpectedUtility(s);
      tmp[s.hash] = s;
      ++i;
      if (i % 10000 == 0) 
        cout << "i: " << i << endl;
    }
    int j = 0;
    for (const auto& p : stateMap) {
      j++;
      if (j < 10) break;
      cout << j << ": " << p.second.u << endl;
    }
    if (stateMap == tmp) break;
    stateMap = tmp; // copy
    nIter++;
    cout << nIter << endl;
  }
}

int main(int argc, char* const argv[]) {
  const struct option longopts[] = {
    { "help", 0, 0, 'h' },
    { "N", 1, 0, 'N' },
    { "M", 1, 0, 'M' },
    { "E1", 1, 0, 'E' },
    { "E2", 1, 0, 'F' },
    { "E3", 1, 0, 'G' },
    { "W", 1, 0, 'W' },
    { "P", 1, 0, 'P' },
    { "dist", 1, 0, 'D' },
    { "lambda", 1, 0, 'l' },
    { "mu", 1, 0, 'm' },
    { "sigma", 1, 0, 's' },
    { 0, 0, 0, 0 }
  };
  char c = '\0';
  int index = 0;
  while ((c = getopt_long(argc, argv, "hN:M:E:F:G:W:P:D:l:m:s:", longopts, &index)) != -1) {
    switch (c) {
      case 'h':
        cout << "Usage: " << argv[0] << endl;
        for (int i = 0; longopts[i].name; ++i) {
          const string name = longopts[i].name;
          const string val = name == "dist" ? "<normal|poisson>" :
                             name == "help" ? "" : "<value>";
          cout << "\t--" << longopts[i].name << " " << val << endl;
        }
        return 0;
      case 'N':
        N = atoi(optarg);
        break;
      case 'M':
        M = atoi(optarg);
        break;
      case 'E':
        E1 = atof(optarg);
        break;
      case 'F':
        E2 = atof(optarg);
        break;
      case 'G':
        E3 = atof(optarg);
        break;
      case 'W':
        W = atof(optarg);
        break;
      case 'P':
        P = atof(optarg);
        break;
      case 'D':
        isNormal = strcmp("normal", optarg) == 0;
        if (!isNormal && strcmp("poisson", optarg) != 0) {
          cerr << "Error: distribution must be either 'normal' or 'poisson'" << endl;
          return -1;
        }
        break;
      case 'l':
        lambda = atof(optarg);
        break;
      case 'm':
        mu = atof(optarg);
        break;
      case 's':
        sigma = atof(optarg);
        break;
    }
  }
  cout << "N: " << N << endl
       << "M: " << M << endl
       << "E1: " << E1 << endl
       << "E2: " << E2 << endl
       << "E3: " << E3 << endl
       << "W: " << W << endl
       << "P: " << P << endl
       << "dist: " << (isNormal ? "normal" : "poisson") << endl
       << "lambda: " << lambda << endl
       << "mu: " << mu << endl
       << "sigma: " << sigma << endl;

  if (isNormal) {
    for (int i = mu; ; ++i) {
      double p = normal_pmf(i, mu, sigma);
      if (p < 0.001) break;
      maxRejected = i;
    }
  } else {
    for (int i = lambda; ; ++i) {
      double p = poisson_pmf(i, lambda);
      if (p < 0.001) break;
      maxRejected = i;
    }
  }

  cout << "maxRejected: " << maxRejected << endl;
  populateStates();
  cout << "numStates: " << stateMap.size() << endl;
  valueIteration();
  return 0;
}
