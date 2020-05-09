# Baer-Nunziato
Solve the Baer-Nunziato equations using second order (TVD) path conservative finite volume method. The code uses either LLF (Rusanov) riemann solver or the more recent HLLEM riemann solver. To compile the code run the following command:

* g++ -Wall -O2 hype.cc

The code depends on Boost multiarray class, which can be easily installed on most linux systems. For example on
* Ubuntu, Debian, Linux Mint: sudo apt install libboost-all-dev
* Fedora, CentOS, Redhat: sudo dnf install boost-devel

## References
1. [A new efficient formulation of the HLLEM Riemann solver for general conservative and non-conservative hyperbolic systems](https://www.sciencedirect.com/science/article/pii/S0021999115006786), Michael Dumbser, Dinshaw S. Balsara, Journal of Computational Physics, Volume 304, 1 January 2016, Pages 275-319
