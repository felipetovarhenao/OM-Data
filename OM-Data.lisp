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
            Optimal-sorting
            ) nil)
        ("2-Metrics" nil nil (
            Euclidean-distance
            List-moments
            List-Zscore
            ) nil)
        ("3-List processing" nil nil (
            X-interpolation
            List-quantize
            List-mod
            Posn-map
            Nth-wrand
            List-wrap
            List-fold
            Nested-position
            Nested-nth
            N-occurances
            Unique-seq
            ) nil)
        ("4-Midicents" nil nil (
            Distortion
            Fill-range
            Shift-posn
            Chroma-count
            IC-vector
            Mc-clip
            ) nil)))