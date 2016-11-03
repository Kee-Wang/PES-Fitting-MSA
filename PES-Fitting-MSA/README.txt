This is the MSA software that fits the potential energy surface using
fitting bases that are invariant with respect to permutations of like
atoms. The detailed theory is described in the paper
Xie, Z.; Bowman, J. M. J. Chem. Theory Comput. 2010, 6, 26-34

Requirements:
To use this software, you need a Fortran 90 compiler, a C++ compiler, 
and the "dgelss" subroutine from LAPACK. Also Perl is required.
We've tested on Linux operating system, using "ifort" and "g++" as the
Fortran and C++ compilers, and Intel MKL library for the "dgelss".
Other compilers should also work, but you may need to use different
flags. LAPACK is available at
www.netlib.org,
and is also included in the optimized libraries such as ACML and MKL.

Here is a list of folders or files in the package:
    (A) A folder called "emsa", which contains the C++ codes that
        generate the monomials and polynomials.
    (B) Two Fortran code templates "fit.f90" and "pes_shell.f90". You
        need to modify these two files later.
    (C) One example Fortran program "getpot.f90" that uses the fitted
        potential energy surface to calculate the potential of a given
        geometry.
    (D) A Makefile that is used to compile the Fortran codes.
    (E) A example data file "points.dat" that contains 44623 geometries
        of H2-H2O and the corresponding interaction energies. This is
        the database used to fit the potential.
    (F) Test data "test.xyz" that contains an arbitrary geometry of
        the H2-H2O. This test data is used by the "getpot" program
        and the expected output is written in the file "expected.out".
    (G) This "README" file.

Here are the procedure and a few explanations for using the software.
The procedure is also described in our tutorial video.
www.youtube.com/watch?v=WzChdPWc7tQ
We recommend that you watch the video first and get familiar with this
software.

(1) Go to "emsa" folder, and compile the source code by using the
    Makefile in this folder. A executable "msa" should be generated.
    This Makefile uses g++ compiler for C++ and works on Linux machines.
    If you are using other compilers or operating systems, you may need
    to modify this Makefile.

(2) Move the executable to the main folder "PES-Fitting-MSA-master" and
    run it as:
    ./msa [max order] [molecule type]
    for example: ./msa 4 2 2 1
    (A) 4 is the polynomial order. Typically we use higher orders.
        But for this demonstration, we just use 4 because it's fast.
        The number of coefficients increases rapidly with the polynomial
        order. So you may want to start with low polynomial orders.
        Also, we always want the number of coefficients to be
        sufficiently smaller than the number of ab initio energies
        to avoid overfitting.
    (B) 2 2 1 is the formula A2B2C of the system. This fitting basis
        works for all the systems with this general formula.
        The example system is H2-H2O. The full symmetry is 4 1, because
        there are 4 hydrogen atoms. But since we don't expect any H
        exchange between H2 and H2O, 2 2 1 is also a reasonable choice.

(4) Run postemsa.pl and derivative.pl:
    ./postemsa.pl 4 2 2 1
    The Perl script reads the output of msa (MOL_2_2_1_4.*) and creates
    the Fortran file called basis.f90, which defines the monomials
    and polynomials as the fitting bases

    ./derivative.pl 4 2 2 1
    It also creates a Fortran file gradient.f90, which defines the
    first derivative of the monomials and polynomials with respect to
    Cartesian coordinates.

(5) Modify the fit.f90 code properly:
    (A) change XXX to the number of atoms
    (B) change YYY to the number of data points
    (C) change ZZZ to the number of coefficients
        (basis.f90 tells the number of coeff.)

(6) Compile the fitting code by typing:
    make fit.x
    Here another Makefile to compile the Fortran code is provided,
    and depending on your operating system and compiler, you may need
    to modify it. This Makefile uses "ifort" as the Fortran compiler.
    The fit.f90 uses "dgelss" from the linear algebra package (LAPACK)
    for the standard linear least-squares problem.

(7) The points.dat file contains geometries and corresponding energies.
    First line is the number of atoms, and second line is the energy in
    Hartree. All the following lines are the Cartesian coordinates in
    Angstrom. The order of the atoms must agree with A2B2C
    In our example, the energy is the interaction between H2 and H2O,
    which is the difference between the dimer energy and the monomer
    energies. But our fitting software can be used for both interaction
    and full potential.

(8) Run the fitting code: ./fit.x
    All the coefficients will be in "coeff.dat"

(9) Modify the pes_shell.f90:
    pes_shell.f90 defines the subroutine to use the PES. It is also a
    template.
    (A) replace XXX with the number of coefficients (size of c in
        basis.f90)
    (B) replace YYY with the number of polynomials (size of p in
        basis.f90, should be the same as the number of coefficients)
    (C) replace ZZZ with the number of monomials (size of m in
        basis.f90)

(10) Compile the getpot using the Makefile:
     make getpot.x

(11) Run the test:
     ./getpot.x test.xyz
     The energy and first derivative of energy (gradient) will be
     written in the file test.out. It should be the same as the
     expected.out (of course small numerical difference like 10e-7 or
     10e-8 is expected)
