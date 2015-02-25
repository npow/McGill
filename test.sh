#!/bin/bash
set -x
for N in `jot 10`; do
  for M in `jot 10`; do
    for E1 in `jot 10`; do
      for E2 in `jot 10`; do
        for E3 in `jot 10`; do
          for P in `jot 10`; do
            for W in `jot 10`; do
              for lambda in `jot 10`; do
                for sigma in `jot 3`; do
                  for mu in `jot 10`; do 
                    ./main --N $N --M $M --E1 $E1 --E2 $E2 --E3 $E3 --W $W --P $P --dist poisson --mu $mu --sigma $sigma --lambda $lambda --iter policy > `mktemp -t npow`
                  done
                done
              done
            done
          done
        done
      done
    done
  done
done
