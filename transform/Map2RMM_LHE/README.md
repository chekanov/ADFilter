# Map2RMM

This is a library for mapping collision data to the rapidity-mass matrices (RMM) for machine learning.
This example uses input Monte Carlo data in the form of ProMC files from the [HepSim Monte Carlo repository](http://atlaswww.hep.anl.gov/hepsim/). The example builds anti-KT jets, and fill ROOT histograms.

The output RMM  is "4t3n" since we use 4 types of objects (jet, e, mu, gamma) up to multiplicity 3.
You can redefine the type of the RMM using these values inside example.cc: 

```
   const int maxNumber=3; // max number for each object (MET is not counted)
   const int maxTypes=4;  // max numbers of types (met not counted)
```

To increase the multiplicity for each object, simply increase the value of maxNumber.
If you need to go beyond 4 types in this example, you would need to extend the function in map2rmm/src/map2rmm.cxx.


## Installation 

This installation instruction works for any Linux computer with gcc4 and above.


 1. Install ProMC (http://atlaswww.hep.anl.gov/asc/promc/), ROOT and FastJet 
 2. Check the installation. The variables: 

```
   echo $PROMC
   echo $ROOTSYS
   echo $FASTJET
```
  they should return the installation paths. 

 3. Go to "map2rmm/" and compile the library as "make"
    This create 2 libraries :

```
      lib/libmap2rmm.so
      lib/libmap2rmm_static.a
```

 4. Go to the upper level and compile the example.cc  "make". The compilation will link the above library

 5. Download ProMC files from HepSim and put them to the "data" directory. Use hs-tools as: 
  
``` 
   hs-get tev100_higgs_ttbar_mg5 data
```
   See the HepSim documentation. 

 6. Process all files inside the directory "data" using the command "./example".

 7. Look at the output root file "output.root" with histograms.
    The RMM data are stored as a tree "inputNN" with "proj" branch. We store only non-zero values, row and column indicies. 

## Output

The output with RMM matrices is "output.root". 
Typically,  you need "inputNN" tree which contains non-zero values of RMM (with indices). We keep only non-zero values since the RMM matrix is sparse. You need to add zeros if you want to get the real matrix.
The script "prepare.py" gives an example of how to unpack the matrix stored inside this root file.
  
## Reference

This algorithm is described in: S.Chekanov "Imaging particle collision data for event classification using machine learning",  Nucl Inst. and Meth. in Phys. Research (NIMA), A931 (2019) p92 (https://arxiv.org/abs/1805.11650).

Different use cases for classification of experimental data using this algorithm were discussed
in the article: 
"Machine learning using rapidity-mass matrices for event classification problems in HEP", 
(S.V.Chekanov),  ANL-HEP-147750, https://arxiv.org/abs/1810.06669
 

S.Chekanov (ANL) 

