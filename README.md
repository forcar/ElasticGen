# ElasticGen
Original author: Ralph Minehart (University of Virginia)
Modifications by K. Joo and L.C. Smith (UVA)

Java port of Fortran (e,e'p) elastic event generator package [elastgen](https://github.com/forcar/elastgen).  

This Java package is currently not setup as an event generator, but contains a single class [MoTsai.java](https://github.com/forcar/ElasticGen/blob/main/src/main/java/org/clas/lib/MoTsai.java) which includes calculations based on the Mo-Tsai paper
[RMP, Vol. 41, 205 (1969)](https://github.com/forcar/elastgen/blob/master/pdf/RevModPhys.41.205.pdf).  The code can be used to calculate radiative corrections (including external straggling in the target).  The method [radcor](https://github.com/forcar/ElasticGen/blob/b17b74fa3a60b603e2b1ed5198c4e87b602c8096/src/main/java/org/clas/lib/MoTsai.java#L216) codes the radiative correction equations (II.6) and (II.9):

<img width="1088" alt="MoTsai 1" src="https://github.com/user-attachments/assets/60f3293a-d647-41d7-a805-518f687e5994">


and the method [radtail](https://github.com/forcar/ElasticGen/blob/b17b74fa3a60b603e2b1ed5198c4e87b602c8096/src/main/java/org/clas/lib/MoTsai.java#L348) codes internal bremmstrahlung (radiative tails) calculated using equation (B.4):

<img width="1028" alt="Screenshot 2024-11-16 at 2 27 15 PM" src="https://github.com/user-attachments/assets/9c295729-eda5-4d17-a783-4860f5152054">
