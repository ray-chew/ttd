A simple benchmarking of calculations involving a tensor-train decomposition (TTD) against a simple-tensor. This is 
based on the Xerus library[1], and on the nearest-neighbour interaction system (NNIS) example [published 
here](https://libxerus.org/cascade/), and mentioned in-depth in the paper by Gelß et. al.[2].

* *nnis.cxx* is the TTD example. The complete source code is published on the [Xerus 
website](https://libxerus.org/cascade/).

* *nnis_noTT.cxx* is my attempt at using the Xerus library to solve the simple-tensor format of the same system.
* The accuracy are similar. See *plot_noTT_long.pdf* or, for the TT-format vs the simple-tensor,
    Residual: 0.00134 vs 0.00133,
    Norm:  0.976103 vs 0.976128.
* However, the TTD calculations are about 500 times faster.
    |&nbsp;    | TT-format| simple-tensor|
    | ---- | :----: | :----: |
    | Run 1 | 1.751 | 874.2 | 
    | Run 2 | 1.568 | 887.6 |
    | Run 3 |  1.550 | 911.4 |
    | **Avg. time** |  **1.623s** | **891.1s** |
* The advantage of the tensor-train decomposition for large arrays is significant.


*References*:
[1] Huber, B. & Wolf, S. Xerus - A General Purpose Tensor Library. *https://libxerus.org*, (2014–2018).
[2] Gelß, P., Klus, S., Matera, S., & Schütte, C. Nearest-Neighbor Interaction Systems in the Tensor-Train Format. 
*Journal of Computational Physics 341*, (2017), 140-162.
[3] Oseledets, I. V. Tensor-Train Decomposition. *SIAM Journal on Scientific Computing 33*, 5 (2011), 2295-2317.


