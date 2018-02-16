METROPOLIS_HASTINGS.ZIP
-----------------------------------------------------------------------------------------------------
*Created by:                  Date:              Comment:
 Felipe Uribe Castillo        August 2011        Final project algorithms
 felipe@dagobah:~$
*Mail:
 felipeuribecastillo@gmx.com
*University:
 National university of Colombia at Manizales. Civil engineering Dept.
-----------------------------------------------------------------------------------------------------
This folder contains 3 basic programs and a function which is needed to run the programs.
Below, a brief description of them is presented,
1. metropolis.m
This program develops a very basic example, for the sampling of functions by means of Metropolis algorithm; showing the correlograms and the histogram of the generated samples. Additionally, it should be noted that this algorithm considers only symmetric proposals PDFs.
It does not need the "MH_routine.m" function.

2. metropolis_hastings.m
This program develops 4 examples, for the sampling of complex functions by means of Metropolis-Hastings algorithm, showing the correlograms and the histograms of the generated samples. In this case the proposals PDF its no longer symmetric. Additionally, the burn-in period, the lag period and the Geweke test have been implemented.
It needs the "MH_routine.m" function.

3. metropolis_hastings2.m
This program develops 1 example, for the sampling of a bivariate Gaussian PDF by means of Metropolis-Hastings algorithm, showing the correlograms and the histograms of the generated samples, and the function with its contours and marginals PDF. Additionally, the burn-in period, the lag period and the Geweke test have been implemented.
It needs the "MH_routine.m" function.

4. MH_routine.m (function)
This function is required in most attached programs, since it runs the Metropolis-Hastings algorithm. Through it, the required samples of the target distribution are generated. This function does not have many comments and clarifications, as they appear in each of the other programs.
-----------------------------------------------------------------------------------------------------
Thanks for the download, I hope that they can serve those who are beginning with the study of MCMC algorithms. Comments, clarifications and corrections will be gladly received...

(c)2012 Felipe Uribe.
Civil Engineering Student.
http://felipegamma.blogspot.com/
-----------------------------------------------------------------------------------------------------
References:
1."Markov chain Monte Carlo and Gibbs sampling"
   B.Walsh ----- Lecture notes for EEB 581, version april 2004.
   http://web.mit.edu/~wingated/www/introductions/mcmc-gibbs-intro.pdf
2."The Metropolis-Hastings algorithm"
   Dan Navarro & Amy Perfors.
   COMPSCI 3016: Computational cognitive science. University of Adelaide.
   http://www.psychology.adelaide.edu.au/personalpages/staff/danielnavarro/ccs2011.html
3."An introduction to MCMC for machine learning"
   C.Andrieu, N.De Freitas, A.Doucet & M.I.Jordan.
   Machine Learning, 50, 5-43, 2003.
   http://www.cs.ubc.ca/~nando/papers/mlintro.pdf