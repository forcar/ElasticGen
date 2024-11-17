# ElasticGen
**Java port of Fortran (e,e'p) elastic event generator package [elastgen](https://github.com/forcar/elastgen).** 

Original author: **Ralph Minehart** (University of Virginia)

Modifications by K. Joo and L.C. Smith (UVA)

## Installation
```
git clone https://github.com/forcar/ElasticGen
cd ElasticGen
./build.sh
alias elas <path to ElasticGen>/bin/elas
```
## Usage
Invoke [demo](https://github.com/forcar/ElasticGen/blob/5eec9a2ef1166d2c40f0f95220de8bb854f5cee3/src/main/java/org/clas/lib/MoTsai.java#L463) which reproduces Table I in [Mo-Tsai](https://github.com/forcar/elastgen/blob/master/pdf/RevModPhys.41.205.pdf) page 209 
```
ls-imac.lan 52: elas demo
     Ebeam      Angle    Eelec      -q2         Z0         Z1         Z2
    17.314       35.1    3.975  25.0311    -0.2297    -0.0467    -0.0411
    15.999       19.7   8.0074  14.9965    -0.2516    -0.0262    -0.0307
    14.649       18.8    7.992  12.4921    -0.2522    -0.0236     -0.028
    13.329       17.6   8.0055   9.9897    -0.2524    -0.0206    -0.0247
    11.999     16.082   7.9969   7.5101    -0.2518    -0.0173    -0.0207
    10.723       14.0   8.0054   5.0997    -0.2499    -0.0137     -0.016
     6.032     17.186   4.6867   2.5245    -0.2403    -0.0122    -0.0106
     2.201     38.601   1.4552   1.3996    -0.2257    -0.0161    -0.0083
     2.206     15.999   2.0219   0.3455    -0.2139    -0.0056    -0.0022
```
Generate [table](https://github.com/forcar/ElasticGen/blob/b17b74fa3a60b603e2b1ed5198c4e87b602c8096/src/main/java/org/clas/lib/MoTsai.java#L463) as follows: `elas table <Ebeam> <theta_min> <thata_max> <theta_bin_width>` 
```
ls-imac.lan 53: elas table 7.645 5 30 2
   Ebeam Angle Eelec   -q2 xsraw(nb) xsrad(nb)  rc_int  rc_ext      rc     bcc
   7.645     5 7.415 0.431  5875.699  4518.491  -0.227  -0.036   0.769   0.775
   7.645     7 7.207 0.821   635.191   482.683  -0.239  -0.036    0.76   0.835
   7.645     9 6.948 1.308    97.789    73.628  -0.248  -0.036   0.753   0.877
   7.645    11  6.65 1.868    19.358    14.466  -0.256  -0.036   0.747   0.902
   7.645    13 6.324 2.478      4.68     3.475  -0.262  -0.036   0.743   0.918
   7.645    15 5.984 3.117     1.339     0.989  -0.268  -0.035   0.739   0.943
   7.645    17 5.638 3.767     0.442     0.325  -0.273  -0.035   0.735   0.952
   7.645    19 5.295 4.411     0.166     0.121  -0.277  -0.035   0.732   0.961
   7.645    21  4.96 5.038     0.069      0.05   -0.28  -0.035    0.73   0.966
   7.645    23  4.64  5.64     0.032     0.023  -0.283  -0.034   0.728   0.978
   7.645    25 4.335 6.211     0.016     0.011  -0.286  -0.034   0.726   0.976
   7.645    27 4.049 6.748     0.008     0.006  -0.288  -0.034   0.725   0.984
   7.645    29 3.782  7.25     0.005     0.003   -0.29  -0.033   0.724   0.985
```
## Notes

This Java package is currently not setup as an event generator, but contains a single class [MoTsai.java](https://github.com/forcar/ElasticGen/blob/main/src/main/java/org/clas/lib/MoTsai.java) which includes calculations from the Mo-Tsai paper
[RMP, Vol. 41, 205 (1969)](https://github.com/forcar/elastgen/blob/master/pdf/RevModPhys.41.205.pdf).  The code can be used to calculate radiative corrections (including external straggling in the target).  

The method [radcor](https://github.com/forcar/ElasticGen/blob/b17b74fa3a60b603e2b1ed5198c4e87b602c8096/src/main/java/org/clas/lib/MoTsai.java#L216) codes the radiative correction equations (II.6) and (II.9):

<img width="1088" alt="MoTsai 1" src="https://github.com/user-attachments/assets/60f3293a-d647-41d7-a805-518f687e5994">


The method [radtail](https://github.com/forcar/ElasticGen/blob/b17b74fa3a60b603e2b1ed5198c4e87b602c8096/src/main/java/org/clas/lib/MoTsai.java#L348) calculates internal bremmstrahlung (radiative tails) by evaluating the integrand of equation (B.4), which is sampled in the Fortran event generator [elastgen](https://github.com/forcar/elastgen) 

<img width="1028" alt="Screenshot 2024-11-16 at 2 27 15 PM" src="https://github.com/user-attachments/assets/9c295729-eda5-4d17-a783-4860f5152054">
