
*****************************************************
*
*		MATAMG - MATLAB Algebraic Multigrid toolbox
*
*	Version v.0.1.4
*	Copyright(C) 2009-2011, Minho Park
*	
*****************************************************

			CONTENTS
===============================
1 Introduction
2 Getting Started
   2.1 Obtaining a copy of MATAMG
   2.2 Installing MATAMG
   2.3 Running MATAMG
3. Change Parameter
4. Example Programs
   4.1 5-point stencil Poisson
   4.2 9-point stencil Poisson
5. How to cite MATAMG
6. Bug Reports

1. INTRODUCTION
===============
MATAMG stands for MATLAB Algebraic Multigrid. It is MATLAB toolbox
designed to solve a linear system with algebraic multigrid algorithms. 
MATAMG support classical algebraic multigrid(AMG) interpolation, 
adaptive AMG(aAMG) interpolation and  bootstrap AMG(BAMG) 
interpolation. MATAMG is designed to allow the user to choose several 
kinds of relaxations, coarsening approaches and adaptive processes as well. 

2. GETTING STARTED
==================

2.1 OBTAINING A COPY OF MATAMG
==============================
MATAMG is available at http://grandmaster.colorado.edu/~parkmh/MATAMG/MATAMG.html. 
It can be obtained by e-mailing the author directly(parkmh@colorado.edu).

2.2 INSTALLING MATAMG
=====================
To install MATAMG  toolbox:

1. Unzip the distribution file and move the directory/folder
  to a convenient location. The archive will create a matamg 
  folder and all program files and documents will be unpacked
  into this folder. 

2. Once unpacking is done, you need to add the folder to 
  your MATLAB path. Simply run "install.m" file. After adding the path,
  you can run MATAMG whenever you want. 

2.3 RUNNING MATAMG
==================
To run MATAMG, make sure matamg folder is in MATLAB path. You can 
check this by typing the command path at the command line. 
>>path

First, you need to set option file. Run amgset.m at the commend line.
>>option = amgset;

To see all parameter names and their possible values, with default shown in {}, run
"amgset.m" with no input arguments.
>>amgset

Once option structure is generated, you can run amg.
>>v = amg(A,v,f,option);
where v is initail vector and f is right-hand side one. 

The convergnece history, complexity and amg OPTION are recored in rst file in the current
directory. 

3. CHANGE PARAMETER
For example, if you want to use 5 levels and rBAMG interpolation, use as follows
>> option = amgset(option,'Coarsest',5,'IntpType','rbamg');

4. EXAMPLE PROGRAMS
=====================
In this section, we introduce two examples solved by AMG algorithm.
The problem size is 64x64 and the boundary condition is dirichlet on unit square.

4.1 5- Point Stencil Poisson
=====================
>> load 5poisson
>> v = rand(length(b),1);
>> option = amgset;
>> v = amg(A,v,b,option);

4.2 9-Point Stencil Poisson
===================== 
>> load 9poisson
>> v = rand(length(b),1);
>> option = amgset;
>> v = amg(A,v,b,option);


5. HOW TO CITE MATAMG
=====================
Please cite the following document:
 Minho Park, MATAMG: an Algebraic Multigrid toolbox for MATAMG, 2009. 
Software available at http://grandmaster.colorado.edu/~parkmh/MATAMG/MATAMG.html
The bibtex format is 
@Manual{Park09,
  author =       {Minho Park},
  title =        {{MATAMG}: a {AMG} toolbox for {MATAMG}},
  year =         {2009},
  note =         {Software available at \url{http://grandmaster.colorado.edu/~parkmh/MATAMG/MATAMG.html}}
}

6. BUG REPORTS
===============
If you find a bug, either in the code or the documentation,
we are interested in hearing about it. Please send a report to :
Minho Park <parkmh@colorado.edu>
  Website : http://grandmaster.colorado.edu/~parkmh

Minho Park
University of Colorado at Boulder
Department of Applied Mathematics

updated [2010.10.21]. 
