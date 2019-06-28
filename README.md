# Matlab code for the computation of coverage probability and rate coverage of a C-V2X network modeled by a Cox process driven by Poisson line process
---
This repository contains the matlab scripts used to generate the results of the paper [https://arxiv.org/pdf/1901.09236.pdf](https://arxiv.org/pdf/1901.09236.pdf) presented in Section-V. 


Run 'main.m' to compute coverage probability/rate coverage through Monte Carlo simulation of the network. 

Run 'main_thy.m' to evaluate coverage probability/rate coverage using the theoretical expressions derived in the paper. Please note that this code is provided for Rayleigh fading channels.

For general Nakagami-$m$ fading with fading parameters $m_1$ and $m_{20}$, please run 'coverage_prob_theoretical.m' to compute the coverage probability using the analytical expressions.

Email to vishnucr@vt.edu for further questions/issues.
