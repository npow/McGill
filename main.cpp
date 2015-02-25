#include <cassert>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <sstream>
#include <vector>
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
int maxPeople = 0;

double normal_pmf(double x, double m, double s) {
  static const double inv_sqrt_2pi = 0.3989422804014327;
  double a = (x - m) / s;
  return inv_sqrt_2pi / s * exp(-0.5f * a * a);
}

double poisson_pmf(const double k, const double l) {
  return exp(k * log(l) - lgamma(k + 1.0) - l);
}

string getHash(const int o, const int p) {
  stringstream ss;
  ss << o << " " << p;
  return ss.str();
}

struct State {
  State() {}
  State(int o, int p) : o(o), p(p) {
    hash = getHash(o, p);
    u = 0;
  }
  bool operator<(const State& rhs) const {
    return hash < rhs.hash;
  }
  bool operator==(const State& rhs) const {
    return o == rhs.o &&
           p == rhs.p &&
           (abs(u-rhs.u) < 0.01*(1-discountFactor)/discountFactor);
  }
  double getReward() const {
    double reward = open()*E3 + waiting()*W + rejected()*P;
    return -reward;
  }
  int open() const {
    return o;
  }
  int closed() const {
    return M-o;
  }
  int waiting() const {
    int x = p > o ? (p-o > N ? N : p-o) : 0;
    assert(x >= 0);
    return x;
  }
  int rejected() const {
    int x = p > N+o ? p-N-o : 0;
    assert(x >= 0);
    return x;
  }

  int o; // open rooms
  int p; // people that arrived
  double u; // utility
  string hash;
};

typedef unordered_map<string, State> StateMap_t;
StateMap_t stateMap;

void populateStates() {
  for (int o = 0; o <= M; ++o) {
    for (int p = 0; p <= maxPeople; ++p) {
      State s(o, p);
      stateMap[s.hash] = s;
    }
  }
}

double getExpectedUtility(const State& s, const int o) {
  const int opening = o - s.o > 0 ? o - s.o : 0;
  const int closing = o - s.o < 0 ? s.o - o : 0;
  double currUtility = 0;
  for (int p = 0; p <= maxPeople; ++p) {
    double prob = 0;
    if (isNormal) {
      prob = normal_pmf(p, mu, sigma);
    } else {
      prob = poisson_pmf(p, lambda);
    }
    if (prob < 0.001) continue;

    const string hash = getHash(o, p);
    const auto& it = stateMap.find(hash);
    assert(it != stateMap.end());
    if (it == stateMap.end()) {
      continue;
    }
    const double pp = prob * (0 -(E1+E3)*opening - E2*closing + it->second.u);
    assert(pp <= 0);
    currUtility += pp;
  }
  return currUtility;
}

double getBestExpectedUtility(const State& s, int& bestPolicy) {
  double bestUtility = -std::numeric_limits<double>::max();
  for (int o = 0; o <= M; ++o) {
    double currUtility = getExpectedUtility(s, o);
    assert(currUtility > -std::numeric_limits<double>::max());

    if (currUtility > bestUtility) {
      bestUtility = currUtility;
      bestPolicy = o;
    }
  }
  return bestUtility;
}

void policyIteration() {
  int nIter = 1;
  unordered_map<string, int> currPolicy;
  for (const auto& p : stateMap) {
    currPolicy[p.second.hash] = M; // initial policy
  }
  while (true) {
    StateMap_t tmp;
    for (const auto& p : stateMap) {
      State s = p.second; // copy
      assert(currPolicy.find(s.hash) != currPolicy.end());
      const int o = currPolicy[s.hash];
      s.u = s.getReward() + discountFactor * getExpectedUtility(s, o);
      assert(s.u <= 0);
      tmp[s.hash] = s;
    }
    stateMap = tmp; // copy

    bool unchanged = true;
    for (const auto& p : stateMap) {
      const double bestVal = getExpectedUtility(p.second, currPolicy[p.second.hash]);
      int o = 0;
      if (getBestExpectedUtility(p.second, o) > bestVal) {
        currPolicy[p.second.hash] = o;
        unchanged = false;
      }
    }
    if (unchanged) {
      break;
    }

    nIter++;
    cout << '\r' << "iter: " << nIter << flush;
  }
}

void valueIteration() {
  int nIter = 1;
  cout << " ";
  while (true) {
    stringstream ss;
    StateMap_t tmp;
    for (const auto& p : stateMap) {
      State s = p.second; // copy
      int policy = 0;
      s.u = s.getReward() + discountFactor * getBestExpectedUtility(s, policy);
      assert(s.u <= 0);
      tmp[s.hash] = s;
    }
    if (stateMap == tmp) break;
    stateMap = tmp; // copy
    nIter++;
    cout << '\r' << "iter: " << nIter << flush;
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
  while ((c = getopt_long(argc, argv, "hN:M:E:F:G:W:P:D:f:l:m:s:", longopts, &index)) != -1) {
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
      case 'f':
        discountFactor = atof(optarg);
        assert(discountFactor < 1 && discountFactor >= 0);
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
       << "iter: " << (isValue ? "value" : "policy") << endl
       << "dist: " << (isNormal ? "normal" : "poisson") << endl
       << "lambda: " << lambda << endl
       << "mu: " << mu << endl
       << "sigma: " << sigma << endl;

  if (isNormal) {
    for (int i = mu; ; ++i) {
      double prob = normal_pmf(i, mu, sigma);
      if (prob < 0.001) break;
      maxPeople = i;
    }
  } else {
    for (int i = lambda; ; ++i) {
      double prob = poisson_pmf(i, lambda);
      if (prob < 0.001) break;
      maxPeople = i;
    }
  }

  cout << "maxPeople: " << maxPeople << endl;
  populateStates();
  cout << "numStates: " << stateMap.size() << endl;
  if (isValue) {
    valueIteration();
  } else {
    policyIteration();
  }

  return 0;
  cout << endl << "Converged! Enter a state to view the optimal policy." << endl
       << "eg. <open_rooms> <num_arriving_people>" << endl;
  string str;
  while (getline(cin, str)) {
    if (str.empty()) continue;
    stringstream ss(str);
    int o, p;
    ss >> o >> p;
    State s(o, p);
    if (stateMap.find(s.hash) == stateMap.end()) {
      cout << "Error! Invalid state: " << str << endl;
    } else {
      int o = 0;
      double d = getBestExpectedUtility(s, o);
      const int opening = o - s.o > 0 ? o - s.o : 0;
      const int closing = o - s.o < 0 ? s.o - o : 0;
      stringstream ss;
      ss << "OPEN=" << opening << " CLOSE=" << closing;

      cout << "state: " << s.hash << ", utility: " << d << ", policy: " << ss.str() << endl;
    }
  }
  return 0;
}
