MSThesis
========

**Masters Thesis:** Approximating eigenvalues of large stochastic matrices<br/>
**Author:** William DeMeo<br/>
**Advisor:** Jonathan Goodman<br/>
**Institution:** Courant Institute of Mathematical Sciences, NYU<br/>
**Description:** A dissertation submitted in partial satisfaction of the requirements for the degree of
Master of Science in Mathematics in the Graduate Division of the New York University.


Files in the repository
-----------------------
The thesis is in the file **msthesis.pdf**.  

All of the LaTeX source code is contained in a single file, **msthesis.tex**, which can be compiled with the following command:

    pdflatex msthesis.tex
    
This command will produce the file msthesis.pdf, assuming you have a resonably good installation of the TeXLive package or its equivalent. 

The computer program files are in the **src** directory.  I have not tested them recently.  If you try them and encounter a problem, *please [submit an issue](https://github.com/williamdemeo/MSThesis/issues)*.


Abstract
--------
The rate at which a Markov chain converges to a given probability distribution
has long been an active area of research. Well known bounds on this rate of
convergence involve the subdominant eigenvalue of the chain's underlying
transition probability matrix. However, many transition probability matrices are
so large that we are unable to store even a vector of the matrix in fast
computer memory. Thus, traditional methods for approximating eigenvalues are
rendered useless. 

In this paper we demonstrate that, if the Markov chain is reversible, and we
understand the structure of the chain, we can derive the coefficients of the
traditional Lanczos algorithm without storing a single vector. We do this by
considering the variational properties of observables on the chain's state
space. In the process we present the classical theory which relates the
information contained in the Lanczos coefficients to the eigenvalues of the
Markov chain. 
