#|--------------------------------- * 
|     OM-Data library functions     |
|   [www.felipe-tovar-henao.com]    |
|               2021                |
* --------------------------------- *
|#

(in-package :om)

(compile&load (om-relative-path '("sources") "OM-Data-functions")) ; includes the functions in 'sources' folder

(om::fill-library 
    '(
        ("1-Classification" nil nil (
            NNS
            K-means
            DTW
            DTW-align
            KDTree
            KNN
            Optimal-sorting
            Markov-build
            Markov-run
            ) nil)
        ("2-Metrics" nil nil (
            Euclidean-distance
            List-moments
            List-median
            List-mode
            List-Zscore
            ) nil)
        ("3-List processing" nil nil (
            X-interpolation
            List-quantize
            List-mod
            Posn-map
            Pick-random
            List-wrap
            List-fold
            Nested-position
            Nested-nth
            N-occurances
            Unique-seq
            Unique-scramble
            List-frames
            Deep-nth
            Deep-replace
            Histogram
            ) nil)
        ("4-Midicents" nil nil (
            Distortion
            Fill-range
            Shift-posn
            Chroma-count
            IC-vector
            IC-cycle
            Mc-clip
            Mc-wrap
            ) nil)
        ("5-Rhythm" nil nil (
            Rhythmicon
            Euclidean-rhythm
            Risset-rhythm
            ) nil)
        ("6-Score" nil nil (
            Segment-seq
            Extract-channel
            Multi-join
            Get-transients
            ) nil)
        ("7-Other" nil nil (
            L-system
            2D-Turtle
            ) nil)
    )
)

(print "
    ***************************
        OM-Data successfully loaded!

       [www.felipe-tovar-henao.com]
                        2021
    ***************************
")
