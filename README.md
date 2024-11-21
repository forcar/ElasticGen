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
Generate [table](https://github.com/forcar/ElasticGen/blob/8c2ecbeadb4cd9a7fb3304a4d09f1e3660a5e7b7/src/main/java/org/clas/lib/MoTsai.java#L470) as follows: `elas table <Ebeam> <theta_min> <theta_max> <theta_bin_width> <wcut> <ff>`
```
ls-imac.lan 94: elas table 7.546 5 30 2 1.05 1 
   ebeam theta eprime   -q2   tau xsraw(nb) xsrad(nb)  rc_int  rc_ext      rc     bcc
   7.546     5  7.322  0.42 0.119   6.553E3   5.205E3  -0.198  -0.033   0.794   0.757
   7.546     7  7.119 0.801 0.227   6.832E2   5.371E2  -0.208  -0.032   0.786    0.83
   7.546     9  6.866 1.276 0.362   1.010E2   7.878E1  -0.216  -0.032    0.78    0.87
   7.546    11  6.575 1.823 0.518   1.952E1   1.513E1  -0.223  -0.032   0.775   0.903
   7.546    13  6.256  2.42 0.687   4.697E0   3.622E0  -0.228  -0.032   0.771   0.918
   7.546    15  5.923 3.046 0.865   1.356E0   1.041E0  -0.233  -0.032   0.768   0.937
   7.546    17  5.584 3.682 1.046  4.564E-1  3.491E-1  -0.237  -0.032   0.765   0.952
   7.546    19  5.247 4.314 1.225  1.748E-1  1.333E-1   -0.24  -0.031   0.763   0.964
   7.546    21  4.919  4.93   1.4  7.469E-2  5.680E-2  -0.243  -0.031   0.761   0.969
   7.546    23  4.603 5.523 1.568  3.506E-2  2.661E-2  -0.245  -0.031   0.759   0.975
   7.546    25  4.303 6.085 1.728  1.783E-2  1.351E-2  -0.247   -0.03   0.758   0.982
   7.546    27  4.021 6.615 1.878  9.721E-3  7.356E-3  -0.249   -0.03   0.757   0.985
   7.546    29  3.757  7.11 2.019  5.626E-3  4.253E-3   -0.25   -0.03   0.756   0.989
```
Generate [ff](https://github.com/forcar/ElasticGen/blob/f6241b00b4625e0e54573b1bf45d143a71cc3771/src/main/java/org/clas/lib/MoTsai.java#L562) comparison as follows: `elas ff <Ebeam> <theta_min> <theta_max> <theta_bin_width> <wcut>` 
```
ls-imac.lan 96: elas ff 7.546 5 30 2 1.05 
   ebeam theta   -q2   tau dipole(nb) bosted(nb)  brash(nb)     ye(nb)  ratio1  ratio2  ratio3
   7.546     5  0.42 0.119    6.553E3    6.207E3    6.074E3    6.373E3    1.06    1.08    1.03
   7.546     7 0.801 0.227    6.832E2    6.796E2    6.683E2    6.980E2    1.01    1.02    0.98
   7.546     9 1.276 0.362    1.010E2    1.057E2    1.038E2    1.077E2    0.96    0.97    0.94
   7.546    11 1.823 0.518    1.952E1    2.110E1    2.062E1    2.130E1    0.93    0.95    0.92
   7.546    13  2.42 0.687    4.697E0    5.133E0    4.998E0    5.155E0    0.92    0.94    0.91
   7.546    15 3.046 0.865    1.356E0    1.474E0    1.433E0    1.481E0    0.92    0.95    0.92
   7.546    17 3.682 1.046   4.564E-1   4.886E-1   4.748E-1   4.923E-1    0.93    0.96    0.93
   7.546    19 4.314 1.225   1.748E-1   1.832E-1   1.782E-1   1.856E-1    0.95    0.98    0.94
   7.546    21  4.93   1.4   7.469E-2   7.645E-2   7.454E-2   7.789E-2    0.98       1    0.96
   7.546    23 5.523 1.568   3.506E-2   3.501E-2   3.424E-2   3.587E-2       1    1.02    0.98
   7.546    25 6.085 1.728   1.783E-2   1.738E-2   1.706E-2   1.790E-2    1.03    1.05       1
   7.546    27 6.615 1.878   9.721E-3   9.255E-3   9.116E-3   9.581E-3    1.05    1.07    1.01
   7.546    29  7.11 2.019   5.626E-3   5.239E-3   5.180E-3   5.446E-3    1.07    1.09    1.03
```
## Notes

This Java package is not yet configured as an event generator, but contains a single class [MoTsai.java](https://github.com/forcar/ElasticGen/blob/main/src/main/java/org/clas/lib/MoTsai.java) which includes calculations from the Mo-Tsai paper
[RMP, Vol. 41, 205 (1969)](https://github.com/forcar/elastgen/blob/master/pdf/RevModPhys.41.205.pdf).  The code can be used to calculate radiative corrections (including external straggling in the target).  

The method [radcor](https://github.com/forcar/ElasticGen/blob/8c2ecbeadb4cd9a7fb3304a4d09f1e3660a5e7b7/src/main/java/org/clas/lib/MoTsai.java#L219) codes the radiative correction equations (II.6) and (II.9):

<img width="1088" alt="MoTsai 1" src="https://github.com/user-attachments/assets/60f3293a-d647-41d7-a805-518f687e5994">


The method [radtail](https://github.com/forcar/ElasticGen/blob/8c2ecbeadb4cd9a7fb3304a4d09f1e3660a5e7b7/src/main/java/org/clas/lib/MoTsai.java#L348) calculates internal bremmstrahlung (radiative tails) by evaluating the integrand of equation (B.4), which is sampled in the Fortran event generator [elastgen](https://github.com/forcar/elastgen) :

<img width="1028" alt="Screenshot 2024-11-16 at 2 27 15 PM" src="https://github.com/user-attachments/assets/9c295729-eda5-4d17-a783-4860f5152054">
