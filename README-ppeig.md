## How to use the ppeig branch to compute eigenvalues of a square matrix. 
The two main functions in this branch are `tridiag` and `qrtrd`.

The `tridiag` function takes a square symmetric matrix as input and reduces the matrix to tri-diagonal form.

An example input file to tridiag: 
```
4 4
1 -1 2 1 
-1 2 0 1 
2 0 3 3
1 1 3 4 
4
1 2 3 3
```
(TODO: in the future, remove the last b vector which is not part of the input). 

You can test it by first changing the second paramter in execYaoProtocol to `tridiag` in the file test_linear_system.c 
and then do 
```
bin/test_linear_system 1234 1 linalg_example.in cholesky 0 54 &
bin/test_linear_system 1234 2 linalg_example.in cholesky 0 54
```


The `qrtrd` function takes a symmetric tridiagonal matrix (represetned by two vectors : its diagonal and off-diagonal), 
and applies the QR algorithm with Wilkinson shift to reduce it to diagonal form. An example input file is 

```
4 
1 2 3 4 
3 
1 2 3
```

You can test it by calling commands below

```
bin/test_linear_system 1234 1 qrtrd_example.in cholesky 0 54 & 
bin/test_linear_system 1234 2 qrtrd_example.in cholesky 0 54
```

## Examples 

FIXME

## TODOs

-- combine the two functions into one QR algorithm. 

-- test on larger size matrices and report timings. 

-- Provide examples and usage data. 
