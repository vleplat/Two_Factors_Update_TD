# Two Factors Update Algorithms for Tensor Decompositions

This repository contains the MATLAB implementation of the two-factor update algorithms developed in the paper:

> **Two Factors Update Algorithms for Tensor Decompositions**  
> Valentin Leplat, Anastasia Sozykina, Igor Vorona, Salman Ahmadi-Asl, and Anh-Huy Phan.

The main idea is to update two factor matrices simultaneously instead of updating the factors one by one. The repository currently includes implementations and numerical tests for:

- the Canonical Polyadic Decomposition (CPD);
- the Block-Term Decomposition (BTD) with ranks \((L_r,L_r,1)\).

The paper also discusses extensions of this idea to the Tucker decomposition and the general Block-Term Decomposition.

## Requirements

The experiments require MATLAB and Tensorlab.

A copy of Tensorlab is included in:

```text
Libraries/tensorlab_2016-03-28/
```

The scripts add the repository folders to the MATLAB path automatically. They should therefore be run from the root directory of the repository.

## Repository structure

```text
Libraries/
    External libraries used by the numerical experiments, including Tensorlab.

functions/
    MATLAB implementations of the proposed algorithms and supporting routines.

main_cpd.m
    Synthetic experiment for the Canonical Polyadic Decomposition.

main_ll1.m
    Synthetic experiments for the BTD with ranks (L_r,L_r,1), together with
    a comparison against the `ll1` algorithm from Tensorlab.
```

## Running the experiments

Open MATLAB in the root directory of the repository.

For the CPD experiment, run:

```matlab
main_cpd
```

For the BTD experiments reported in Section 5 of the paper, run:

```matlab
main_ll1
```

The numerical parameters are specified directly in the corresponding scripts. Random seeds are fixed where needed for reproducibility.

## ADMM update order

The inner ADMM routines follow the order used in the convergence analysis:

```text
X  ->  Z  ->  T
```

where:

- `X` is the non-smooth structured update;
- `Z` is the smooth quadratic update;
- `T` is the dual update.

## Citation

The paper is currently under review. Please use the following temporary citation:

```bibtex
@unpublished{leplat_two_factors_update,
  title  = {Two Factors Update Algorithms for Tensor Decompositions},
  author = {Leplat, Valentin and Sozykina, Anastasia and Vorona, Igor
            and Ahmadi-Asl, Salman and Phan, Anh-Huy},
  note   = {Under review}
}
```

The bibliographic information will be updated after publication.

## Acknowledgements

The baseline algorithms used in the numerical experiments are due to their respective authors. In particular, the reference BTD implementation is provided by Tensorlab.
