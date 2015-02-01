#include <cstring>
#include <iostream>
#include <getopt.h>
using namespace std;

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
  char c = '\0';
  int index = 0;
  while ((c = getopt_long(argc, argv, "hN:M:E:F:G:W:P:D:l:m:s:", longopts, &index)) != -1) {
    switch (c) {
      case 'h':
        break;
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

  return 0;
}
