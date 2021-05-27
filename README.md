<div style="text-align:center"><img src="resources/icons/000.png" width="110" height="110"/></div>

# OM-Data library

**OM-Data** is a currently developing library of elementary functions, mostly aimed at data processing and analysis in **OpenMusic**. The library includes functions for data metrics, classification, sorting, and handling, as well as some functions for midicent manipulation.
These include `K-means`, `DTW`, `NNS`, `Markov-build`, `Markov-run`, `List-moments`, `Euclidean-distance`, `Segment-seq`, `Chroma-count`, `IC-vector`, among others.

A list of example patches are included, demonstrating possible musical applications of some these functions.

Since **OM-Data** is still in its early stages, only beta releases will be made at the moment.

For info on how to install external libraries in **OM** visit: https://support.ircam.fr/docs/om-libraries/main/co/OM-Libraries.html 

Once the library is installed and loaded into the **OM** workspace, the example patches will be available in `/Help/Import Tutorial Patches/Libraries/OM-Data`

```.
OM-Data
├── 1-Classification
│   ├── NNS
│   ├── K-Means
│   ├── DTW
│   ├── DTW-align
│   ├── Optimal Sorting
│   ├── KDTree
│   ├── KNN
│   ├── Markov-build
│   └── Markov-run
│
├── 2-Metrics
│   ├── Euclidean-distance
│   ├── List-moments
│   ├── List-median
│   ├── List-mode
│   ├── List-Zscore
│   ├── List-covariance
│   ├── List-correlation
│   ├── List-Histogram
│   └── Plot-points
│
├── 3-List-processing
│   ├── X-interpolation
│   ├── Posn-map
│   ├── Pick-random
│   ├── List-wrap
│   ├── List-fold
│   ├── List-symmetries
│   ├── List-frames
│   ├── List-quantize
│   ├── List-mod
│   ├── Nested-position
│   ├── Nested-nth
│   ├── Deep-replace
│   ├── Deep-nth
│   ├── N-occurances
│   ├── Unique-seq
│   └── Unique-scramble
│
├── 4-Midicents
│   ├── NRT
│   ├── Distortion
│   ├── Fill-range
│   ├── Shift-posn
│   ├── Chroma-count
│   ├── IC-vector
│   ├── IC-cycle
│   ├── Mc-clip
│   └── Mc-wrap
│
├── 5-Rhythm
│   ├── Rhythmicon
│   ├── Risset-rhythm
│   └── Euclidean-rhythm
│
├── 6-Score
│   ├── Segment-seq
│   ├── Extract-channel
│   ├── Multi-join
│   └── Get-transients
│
└── 7-Other
    ├── L-system
    ├── 2D-Turtle
    ├── Make-sieve
    └── 2D-Turtle

 ```