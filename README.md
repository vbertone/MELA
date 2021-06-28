[![Circle CI](https://circleci.com/gh/vbertone/MELA.svg?style=svg)](https://circleci.com/gh/vbertone/MELA)

![alt text](https://raw.githubusercontent.com/vbertone/mela/master/logos/logomela.png "Logo MELA")

# MELA: a Mellin Evolution LibrAry

MELA is a Mellin-space based code for parton distribution function
(PDF) and fragmentation function (FF) DGLAP evolution.
MELA also allows to compute DIS structure functions in
different mass schemes up to NNLO in QCD.

## Download

MELA can be downloaded directly from the github repository:

https://github.com/vbertone/MELA

For the last development version you can clone the master code:

```Shell
git clone https://github.com/vbertone/MELA.git
```

For the latest tag:

```Shell
git tag -l
git checkout tags/tag_name
```

## Installation 

Checkout the code and compile the it using the
standar procedure:

```Shell
cd MELA
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=/your/installation/path/ ..
make && make install
```

## References

- V. Bertone, S. Carrazza, E.R. Nocera, *Reference results for time-like evolution up to O(α_s^3)*, [arXiv:1501.00494](http://arxiv.org/abs/1501.00494)

## Contact Information

Maintainers: Valerio Bertone, Stefano Carrazza, Emanuele Roberto Nocera

