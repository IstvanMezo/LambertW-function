# LambertW-function
This is a C++ code for the complex Lambert W function

Use:
Make an initial complex number by using the <complex> standard C++ library:
    complex<double> z{ 1.73, 2.118 };
Then call the LambertW function as follows:
    cout << LambertW(z,k) <<'\n';
to print out the value of W at the point z on the k-th branch. k must be an integer. If k is left out, it is understood to be 0 (principal branch).
  
The declaration of LambertW reads as
    complex<double> LambertW(complex<double> z, int k = 0);

Example can be found in the source code.
