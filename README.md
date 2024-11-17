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
Invoke [demo](https://github.com/forcar/ElasticGen/blob/8c2ecbeadb4cd9a7fb3304a4d09f1e3660a5e7b7/src/main/java/org/clas/lib/MoTsai.java#L451) which reproduces Table I in [Mo-Tsai](https://github.com/forcar/elastgen/blob/master/pdf/RevModPhys.41.205.pdf) page 209 
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
Generate [table](https://github.com/forcar/ElasticGen/blob/8c2ecbeadb4cd9a7fb3304a4d09f1e3660a5e7b7/src/main/java/org/clas/lib/MoTsai.java#L470) as follows: `elas table <Ebeam> <theta_min> <theta_max> <theta_bin_width> <wcut>`
```
ls-imac.lan 104: elas table 7.546 5 30 2 1.05
   Ebeam Angle Eelec   -q2 xsraw(nb) xsrad(nb)  rc_int  rc_ext      rc     bcc
   7.546     5 7.322  0.42   6.201E3   4.925E3  -0.198  -0.033   0.794   0.763
   7.546     7 7.119 0.801   6.789E2   5.338E2  -0.208  -0.032   0.786   0.837
   7.546     9 6.866 1.276   1.056E2   8.238E1  -0.216  -0.032    0.78   0.879
   7.546    11 6.575 1.823   2.108E1   1.634E1  -0.223  -0.032   0.775   0.912
   7.546    13 6.256  2.42   5.127E0   3.954E0  -0.228  -0.032   0.771   0.919
   7.546    15 5.923 3.046   1.473E0   1.131E0  -0.233  -0.032   0.768   0.941
   7.546    17 5.584 3.682  4.881E-1  3.733E-1  -0.237  -0.032   0.765   0.951
   7.546    19 5.247 4.314  1.830E-1  1.396E-1   -0.24  -0.031   0.762   0.956
   7.546    21 4.919  4.93  7.638E-2  5.809E-2  -0.243  -0.031   0.761   0.974
   7.546    23 4.603 5.523  3.498E-2  2.655E-2  -0.245  -0.031   0.759   0.973
   7.546    25 4.303 6.085  1.736E-2  1.316E-2  -0.247   -0.03   0.758    0.98
   7.546    27 4.021 6.615  9.246E-3  6.996E-3  -0.249   -0.03   0.757   0.983
   7.546    29 3.757  7.11  5.234E-3  3.956E-3   -0.25   -0.03   0.756   0.986
```
## Notes

This Java package is not yet configured as an event generator, but contains a single class [MoTsai.java](https://github.com/forcar/ElasticGen/blob/main/src/main/java/org/clas/lib/MoTsai.java) which includes calculations from the Mo-Tsai paper
[RMP, Vol. 41, 205 (1969)](https://github.com/forcar/elastgen/blob/master/pdf/RevModPhys.41.205.pdf).  The code can be used to calculate radiative corrections (including external straggling in the target).  

The method [radcor](https://github.com/forcar/ElasticGen/blob/4775773439641bcd4d87f13549be66366a58db73/src/main/java/org/clas/lib/MoTsai.java#L219) codes the radiative correction equations (II.6) and (II.9):

<img width="1088" alt="MoTsai 1" src="https://github.com/user-attachments/assets/60f3293a-d647-41d7-a805-518f687e5994">


The method [radtail](https://github.com/forcar/ElasticGen/blob/4775773439641bcd4d87f13549be66366a58db73/src/main/java/org/clas/lib/MoTsai.java#L348) calculates internal bremmstrahlung (radiative tails) by evaluating the integrand of equation (B.4), which is sampled in the Fortran event generator [elastgen](https://github.com/forcar/elastgen) :

<img width="1028" alt="Screenshot 2024-11-16 at 2 27 15 PM" src="https://github.com/user-attachments/assets/9c295729-eda5-4d17-a783-4860f5152054">
