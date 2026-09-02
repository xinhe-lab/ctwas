# Estimate s in `susie_rss` Model Using Regularized LD

The estimated s gives information about the consistency between the z
scores and LD matrix. A larger \\(s\\) means there is a strong
inconsistency between z scores and LD matrix. The “null-mle” method
obtains mle of \\(s\\) under \\(z | R \~ N(0,(1-s)R + s I)\\), \\(0 \< s
\< 1\\). The “null-partialmle” method obtains mle of \\(s\\) under
\\(U^T z | R \~ N(0,s I)\\), in which \\(U\\) is a matrix containing the
of eigenvectors that span the null space of R; that is, the eigenvectors
corresponding to zero eigenvalues of R. The estimated \\(s\\) from
“null-partialmle” could be greater than 1. The “null-pseudomle” method
obtains mle of \\(s\\) under pseudolikelihood \\(L(s) =
\\prod\_{j=1}^{p} p(z\_j | z\_{-j}, s, R)\\), \\(0 \< s \< 1\\).

## Usage

``` r
estimate_s_rss(z, R, n, r_tol = 1e-08, method = "null-mle")
```

## Arguments

  - z:
    
    A p-vector of z scores.

  - R:
    
    A p by p symmetric, positive semidefinite correlation matrix.

  - n:
    
    The sample size. (Optional, but highly recommended.)

  - r\_tol:
    
    Tolerance level for eigenvalue check of positive semidefinite matrix
    of R.

  - method:
    
    a string specifies the method to estimate \\(s\\).

## Value

A number between 0 and 1.
