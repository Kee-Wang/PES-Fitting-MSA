This is the MSA software that fits the potential energy surface using
fitting bases that are invariant with respect to permutations of like
atoms. The detailed theory is described in the paper
Xie, Z.; Bowman, J. M. J. Chem. Theory Comput. 2010, 6, 26-34

Requirements:
To use this software, you need a
(A) a Fortran 90 and a C++ compiler;
(B) "dgelss" subroutine from LAPACK
(C) Perl and Pythond

We've tested on Linux operating system, using "ifort" and "g++" as the
Fortran and C++ compilers, and Intel MKL library for the "dgelss".
Other compilers should also work, but you may need to use different
flags. LAPACK is available at www.netlib.org, and is also included in
the optimized libraries such as ACML and MKL.

After downloading and unzipping the "MSA", the following folders and
files should be in the package:

(A) A folder called "src", which contains the source codes. Inside it
    there are
    (a) A folder called "emsa", which contains the C++ codes that
        generate the monomials and polynomials.
    (b) a test program "getpot.f90" that uses the fitted potential to
        calculate the potential of a given geometry.
    (c) A Makefile that is used to compile the Fortran codes.
    (d) A test data file "test.xyz" that contains an arbitrary geometry
        of the H2-H2O. This test data is used by the "getpot" program
        and the expected output is written in the file "expected.out".
    (e) Two Perl scripts.

(B) A example data file "points.dat" that contains 44623 geometries
    of H2-H2O and the corresponding interaction energies. This is
    the database used to fit the potential.

(C) A Python script that should be executed to fit the potential

(D) A script called "reset" to clean up everything and reset the whole
    MSA directory to its original condition. (Be careful because all
    the output and results will be removed)
    
(E) This README.txt file

Here are the procedure and a few explanations for using the software.
The procedure is also described in our tutorial video.
https://www.youtube.com/watch?v=vGATfUzEeAE
We recommend that you watch the video first and get familiar with this
software.

(0) First take a look at the "points.dat" file.
    It contains geometries and corresponding energies. The first line
    is the number of atoms, and second line is the energy in Hartree.
    All the following lines are the Cartesian coordinates in Angstrom.
    The order of atoms must agree with the molecular formula. In our
    example, it's A2B2C. The energy is the interaction between H2 and
    H2O, but our fitting software can be used for both interaction and
    full potential.

(1) Modify the two Makefile's, based on your compiler and LAPACK.

(2) Run the Python script:
    ./msa.py
    (or python msa.py)

    (A) Input the polynomial order you would like to use for the fitting.
        In our example we use 4. The number of coefficients increases
        rapidly when the polynomial order becomes larger. So you may want
        to start with low polynomial orders.
    (B) Input the molecular formula (or the permutational symmetry group)
        Our example is H2-H2O, and we use 2 2 1 here. The full symmetry is
        4 1, but we don't expect any H exchange between H2 and H2O, so
        2 2 1 is also a reasonable choice.
    (C) Input the name of the data file, which is "points.dat".

(3) The program tells you the number of coefficients and the number of
    points in the data file. And it asks you if you would like to continue.
    If the number of coefficients is too small (which leads to large
    fitting error), or too large (which may cause over-fitting), you can
    type "n" to terminate and then pick another polynomial order. If you
    would like continue, just enter y, and the program will continue.

(4) The program fits the potential energy surface and when it finishes,
    the root-mean-square fitting error is printed on the screen.
    The coefficients of the fit is written in "coeff.dat", and three
    Fortran code files "pes_shell.f90", "basis.f90", and "gradient.f90"
    are also generated.

(5) The test program "getpot.x" is compiled, and if you would like to run
    the test, use the command
    ./getpot.x test.xyz
    The results is written in test.out, and if you use polynomial order 4
    and symmetry group 2 2 1 as we did in the video, the results in the
    "test.out" should be the same as those in "expected.out". Of course,
    small numerical error is allowed.

(6) If you would like to use the fit in your own program, pes_shell.f90,
    basis.f90, gradient.f90, and coeff.dat are necessary. Copy these four
    files to the folder that contains your own program, and in your own
    Fortran code, insert "use pes_shell", and "call pes_init()", (as we
    do in the "getpot.f90" example), and you can calculate the potential
    of any configuration using the "f" function, and the gradient using the
    "g" function.

Credits:
MSA codes: Zhen Xie
gradient extension: Chen Qu
Python script: Qingfeng Kee Wang
