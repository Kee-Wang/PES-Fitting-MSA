# In order to use MSA:

1. Download contents of the git repository.
2. Read readme.txt.

# About MSA

MSA is a software that does a fit of electronic energies using fitting bases that are invariant with respect to permutations of like atoms. The detailed theory is discussed in the Xie and Bowman paper (see reference at the end of file). The software also provides the potential gradient, if desired.




# Credits
MSA codes: Zhen Xie

Gradient extension: Chen Qu

Python script: Qingfeng Kee Wang



# More Info

For more details and video turorial: https://scholarblogs.emory.edu/bowman/msa/

# Contact
kee.wang@emory.edu

# Reference
* Xie, Z., Bowman, J.M. Permutationally Invariant Polynomial Basis for Molecular Energy Surface Fitting via Monomial Symmetrization. J. Chem. Theory Comput. 2010, 6, 26-34.
*  R. Conte, C. Qu, and J. M. Bowman, Permutationally Invariant Fitting of Many-Body, Non-covalent Interactions with Application to Three-Body Methane–Water–Water, R. Conte, C. Qu, and J. M. Bowman, J. Chem. Theory Comput. 2015, 11, 1631-1638.


# Version log

## ver 1.5, Dec. 19, 2016
* Compatibility improved so it works for both Python 2 and Python 3. -- Kee 

## ver 1.4, Dec. 19, 2016
* Added feature so that people can enter their own desired a0 value. Use `sed` to change a0 in `gradient.f90` directly.  -- by Chen Qu

## ver 1.3, Dec. 16, 2016
* a0 value in derivative.pl was fixed from 2d0 to 2d5. -- Thanks Bryan (University of Colorado Boulder) for pointing that out.


## ver 1.2, Nov. 29, 2016

* Updated README.text info so that video link was pointed to group page. -- by Chen Qu
* Minor correction. Corrected energies in sample file -- by Chen Qu

## ver 1.1, Nov. 21, 2016

* Added feature so that people can enter their own weights. -- by Chen Qu

## ver 1.0

* First usable MSA package was distributed. -- by Kee

