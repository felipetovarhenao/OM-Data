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
            Optimal-sorting
            ) nil)
        ("2-Metrics" nil nil (
            Euclidean-distance
            List-mean
            List-variance
            List-stdev
            List-Zscore
            ) nil)
        ("2-List processing" nil nil (
            X-interpolation
            List-quantize
            List-mod
            Posn-map
            Nth-wrand
            List-wrap
            List-fold
            ) nil)
        ("3-Midicents" nil nil (
            Distortion
            Fill-range
            Shift-posn
            ) nil)
    )
) ; 