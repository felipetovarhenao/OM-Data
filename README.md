<div style="text-align:center"><img src="resources/icons/000.png" width="150" height="150"/></div>

# OM-Data library

**OM-Data** is a library of elementary functions, partially aimed at data modelling and analysis in **OpenMusic**. The library includes functions for data metrics, classification, and processing, as well as some higher-level functions for specific musical operations.
These include `K-means`, `DTW`, `NNS`, `KDTree`, `KNN`, `Markov-build`, `Markov-run`,  `Segment-seq`, `Score-filter`, `Get-transients`, `Chroma-count`, `IC-vector`, among many others.

A list of example patches are included, demonstrating possible musical applications for some of these functions.

**OM-Data** is still in development.

For info on how to install external libraries in **OM** visit: https://support.ircam.fr/docs/om-libraries/main/co/OM-Libraries.html 

Once the library is installed and loaded into the **OM** workspace, the example patches will be available in `/Help/Import Tutorial Patches/Libraries/OM-Data`

Updated content tree of the **OM-Data** library:
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
│   ├── Harmonic-distr
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
│   ├── Score-filter
│   ├── Segment-seq
│   ├── Extract-channel
│   ├── Multi-join
│   ├── Granulate
│   └── Get-transients
│
└── 7-Other
    ├── L-system
    ├── 2D-Turtle
    ├── Make-sieve
    ├── Recursivus
    └── Vieru-seq

 ```