#include <cassert>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <map>
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
double discountFactor = 0.9;
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
  ss << o << " " << oi << " " << c << " " << ci << " " << w << " " << r;
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
    double reward = o*E3 + oi*(E1+E3) + ci*E2 + w*W * r*P;
    return -reward;
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

typedef map<string, State> StateMap_t;
StateMap_t stateMap;

void populateStates() {
  for (int numOpen = 0; numOpen <= M; ++numOpen) {
    for (int numOpening = 0; numOpening <= M; ++numOpening) {
      for (int numClosed = 0; numClosed <= M; ++numClosed) {
        for (int numClosing = 0; numClosing <= M; ++numClosing) {
          for (int numWaiting = 0; numWaiting <= N; ++numWaiting) {
            for (int numRejected = 0; numRejected <= maxRejected; ++numRejected) {
              if (numOpen + numOpening + numClosed + numClosing != M) continue;
              if (numWaiting != N && numRejected > 0) continue; // only reject if no seats left
#if 0
              double p = 0;
              if (isNormal) {
                p = normal_pmf(numWaiting+numRejected, mu, sigma);
              } else {
                p = poisson_pmf(numWaiting+numRejected, lambda);
              }
              if (p < 0.001) continue;
#endif
              State s(numOpen, numOpening, numClosed, numClosing, numWaiting, numRejected);
              stateMap[s.hash] = s; 
            }
          }
        }
      }
    }
  }
}

double getBestExpectedUtility(const State& s, string& bestPolicy) {
  double bestUtility = -std::numeric_limits<double>::max();
  for (int closing = 0; closing <= s.o + s.oi; ++closing) {
    for (int opening = 0; opening <= s.c + s.ci; ++opening) {
      int numAvailRooms = s.o + s.oi - s.w;
      int numLeftoverPatients = 0;
      if (numAvailRooms < 0) {
        numLeftoverPatients = abs(numAvailRooms);
        numAvailRooms = 0;
      }
      const int numAvailSeats = N - numLeftoverPatients;
      double currUtility = 0.0;
      for (int i = 0; i <= numAvailRooms+numAvailSeats+maxRejected; ++i) {
        double p = 0;
        if (isNormal) {
          p = normal_pmf(i, mu, sigma);
        } else {
          p = poisson_pmf(i, lambda);
        }
        if (p < 0.001) continue;

        int w = 0;
        int r = 0;
        if (numAvailRooms < i) {
          w = i - numAvailRooms;
          if (numAvailSeats < w) {
            r = w - numAvailSeats;
            w = numAvailSeats;
          } else {
            // have enough seats
          }
        } else {
          // have enough rooms
        }
        w += numLeftoverPatients;
        const string hash = getHash(s.o+s.oi-closing, opening, s.c+s.ci-opening, closing, w, r);
        const auto& it = stateMap.find(hash);
        if (it == stateMap.end()) {
          continue;
        }
        currUtility += it->second.u * p;
      }

      if (currUtility > bestUtility) {
        bestUtility = currUtility;
        stringstream ss;
        ss << "OPEN=" << opening << " CLOSE=" << closing;
        bestPolicy = ss.str();
      }
    }
  }
  return bestUtility;
}

void policyIteration() {

}

void valueIteration() {
  int nIter = 0;
  while (true) {
    StateMap_t tmp;
    for (const auto& p : stateMap) {
      State s = p.second; // copy
      string str = "";
      s.u = s.getReward() + discountFactor * getBestExpectedUtility(s, str);
      tmp[s.hash] = s;
    }
#if 0
    int j = 0;
    for (const auto& p : tmp) {
      j++;
      if (j > 100) break;
      cout << p.second.hash << ": " << p.second.u << endl;
    }
#endif
    if (stateMap == tmp) break;
    stateMap = tmp; // copy
    nIter++;
    cout << "iter: " << nIter << endl;
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
    { "iter", 1, 0, 'I' },
    { "lambda", 1, 0, 'l' },
    { "mu", 1, 0, 'm' },
    { "sigma", 1, 0, 's' },
    { "gamma", 1, 0, 'f' },
    { 0, 0, 0, 0 }
  };
  bool isValue = true;
  char c = '\0';
  int index = 0;
  while ((c = getopt_long(argc, argv, "hN:M:E:F:G:W:P:D:l:m:s:", longopts, &index)) != -1) {
    switch (c) {
      case 'h':
        cout << "Usage: " << argv[0] << endl;
        for (int i = 0; longopts[i].name; ++i) {
          const string name = longopts[i].name;
          const string val = name == "dist" ? "<normal|poisson>" :
                             name == "iter" ? "<value|policy>" :
                             name == "help" ? "" : "<value>";
          cout << "\t--" << longopts[i].name << " " << val << endl;
        }
        return 0;
      case 'N':
        N = atoi(optarg);
        assert(N > 0);
        break;
      case 'M':
        M = atoi(optarg);
        assert(M > 0);
        break;
      case 'E':
        E1 = atof(optarg);
        assert(E1 >= 0);
        break;
      case 'F':
        E2 = atof(optarg);
        assert(E2 >= 0);
        break;
      case 'G':
        E3 = atof(optarg);
        assert(E3 >= 0);
        break;
      case 'W':
        W = atof(optarg);
        assert(W >= 0);
        break;
      case 'P':
        P = atof(optarg);
        assert(P >= 0);
        break;
      case 'D':
        isNormal = strcmp("normal", optarg) == 0;
        if (!isNormal && strcmp("poisson", optarg) != 0) {
          cerr << "Error: distribution must be either 'normal' or 'poisson'" << endl;
          return -1;
        }
        break;
      case 'I':
        isValue = strcmp("value", optarg) == 0;
        if (!isValue && strcmp("policy", optarg) != 0) {
          cerr << "Error: iter must be either 'value' or 'policy'" << endl;
          return -1;
        }
        break;
      case 'l':
        lambda = atof(optarg);
        assert(lambda >= 0);
        break;
      case 'm':
        mu = atof(optarg);
        assert(mu >= 0);
        break;
      case 's':
        sigma = atof(optarg);
        assert(sigma >= 0);
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
  if (isValue) {
    valueIteration();
  } else {
    policyIteration();
  }

  cout << "Converged! Enter a state to view the optimal policy." << endl
       << "eg. <open_rooms> <opening_rooms> <closed_rooms> <closing_rooms> <num_waiting> <num_rejected>" << endl;
  string str;
  while (getline(cin, str)) {
    if (str.empty()) continue;
    stringstream ss(str);
    int o, oi, c, ci, w, r;
    ss >> o >> oi >> c >> ci >> w >> r;
    State s(o, oi, c, ci, w, r);
    if (stateMap.find(s.hash) == stateMap.end()) {
      cout << "Error! Invalid state: " << str << endl;
    } else {
      string bestPolicy;
      double d = getBestExpectedUtility(s, bestPolicy);
      cout << s.hash << "utility: " << d << " policy: " << bestPolicy << endl;
    }
  }
  return 0;
}
