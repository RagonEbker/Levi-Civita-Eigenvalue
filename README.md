# Levi-Civita Poweriteration
This library is an implementation of the poweriteration algorithm for the Levi-Civita field (or the Puiseux field). More specificaly it was created to see if we can calculate the eigenvectors and eigenvalues of the laplacian of a Levi-Civita graph.
## Usage
### Create a graph example
So far we can test the algorithm with laplacian graph examples and through companion matrices of LC polynomials. Of course it is also possible to enter matrix examples by hand. 
To create a laplacian graph example one can use the file `random_graph_example.py` and modify it accordingly.
### Calculate Eigenvalues
The calculation of the eigenvalues can work in two different ways.
1. With a relatively slow taylor expansion in `power_iteration_taylor_example.py`
2. With faster methods that calculate inverse and square roots with the help of newtons method, in `power_iteration_lc_example`.
## Acknowledgment
This implementation came to life through a joint project and many fruitfull discussions with Anna Muranova and Max Schmidt.
