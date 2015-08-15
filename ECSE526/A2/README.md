# ECSE526_A2

##### To build, run `make`. This produces an executable `main`.
```
make
```

##### Usage
```
./main -h                                                                                                                                                                                 [0:05:46]
Usage: ./main
	--help
	--N <value>
	--M <value>
	--E1 <value>
	--E2 <value>
	--E3 <value>
	--W <value>
	--P <value>
	--dist <normal|poisson>
	--iter <value|policy>
	--lambda <value>
	--mu <value>
	--sigma <value>
	--gamma <value>
```

##### Example run
```
./main --M 50 --N 49 --E1 10 --E2 10 --E3 10 --W 1 --P 100 --dist poisson --lambda 2 --iter value --gamma 0.1
N: 49
M: 50
E1: 10
E2: 10
E3: 10
W: 1
P: 100
iter: value
dist: poisson
lambda: 2
mu: -1
sigma: -1
maxPeople: 7
numStates: 408
iter: 4
Converged! Enter a state to view the optimal policy.
eg. <open_rooms> <num_arriving_people>
1 0
state: 1 0, utility: -12.1987, policy: OPEN=0 CLOSE=1
0 0
state: 0 0, utility: -2.20967, policy: OPEN=0 CLOSE=0
```
