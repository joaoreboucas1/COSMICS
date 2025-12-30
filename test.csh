#! /bin/csh -f
# Script for testing cosmics.
#
cd test_results
#
# Testing linger_con with cdm + Lambda model.
#
../bin/linger_con >& test.out <<!
.0355 .3145 0.65 0.
65 2.726 0.24 0
1
.01 100 0.
!
time >>& test.out
\mv -f linger.dat linger.con
\mv -f lingerg.dat lingerg.con
#
# Testing linger_syn with isocurvature seed + cdm model.
#
../bin/linger_syn >>& test.out <<!
.05 .95 0. 0.
50 2.726 0.24 0
0
0.001 10 41 0.
4
!
time >>& test.out
\mv -f linger.dat linger.syn
\mv -f lingerg.dat lingerg.syn
#
# Testing deltat with linger_con results.
#
../bin/deltat >>& test.out <<!
60 0.96
2
lingerg.con
!
time >>& test.out
#
# Testing grafic with linger_syn results.
#
../bin/grafic >>& test.out <<!
1
linger.syn
1.6
-0.8
.001 1
1.0 0.1 0.05
34159265
1
!
time >>& test.out
\mv -f p3m.dat p3m.a
#
# Testing grafic with fit to cdm spectrum.
#
../bin/grafic >>& test.out <<!
2
1.0 0.0 50.0
1.0
20.0
0,0
0.5 0.1 0.05
314159265
1
!
time >>& test.out
\mv -f p3m.dat p3m.b
