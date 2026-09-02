# Get a subset of weights, by group, context or type.

Get a subset of weights, by group, context or type.

## Usage

``` r
subset_weights(weights, names, select_by = c("group", "context", "type"))
```

## Arguments

  - weights:
    
    a list of pre-processed prediction weights.

  - names:
    
    names of groups, contexts, or types to be selected.

  - select\_by:
    
    select weights by group, context, or type.

## Value

a list of selected weights.
