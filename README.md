# Two_Factors_Update_TD
A family of new algorithms to compute tensor decompositions.

This MATLAB software reproduces the results from the following paper:

```
@misc{2fac_updates_2025,
      title={Two Factors Update for Tensor Decompositions}, 
      author={Valentin Leplat and Anastasia Sozykina and Igor Vorona and Salman Ahmadi-Asl and Anh Huy Phan},
      year={2025},
      eprint={},
      archivePrefix={arXiv},
      primaryClass={cs.LG}
}
```

## Acknowledgements

The baseline algorithms used in the manuscript are courtesy of their respective authors.

## Content - MatLab framework "2_fac_updates_CPD_and_LL1"
 
 - /Libraries : contains helpful libaries, such as Tensorlab (https://tensorlab.net) 
 
 - /Results : contains some results presented in Section 5 of the paper.

 - /functions : contains helpful files, the MatLab implementations of the Algorithms developped in the paper, and MatLab routines to run the demos.
   
   
### Test files
 
 Test files are available. To proceed, open and start one of the following files:

- main_cpd.m : the code allowing to compute the CPD with our new algorithm of a synthetic input tensor
- main_ll1.m : the code allowing to compute the CPD with our new algorithm of a synthetic input tensor, and compared with the "ll1" algorithm of Tensorlab, see section 5 of the paper.

