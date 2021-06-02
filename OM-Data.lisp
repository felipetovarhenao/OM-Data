#|--------------------------------- * 
|     OM-Data library functions     |
|   [www.felipe-tovar-henao.com]    |
|               2021                |
* --------------------------------- *
|#

(in-package :om)

; (compile&load (om-relative-path '("sources") "OM-Data-functions" "omd-classification")) ; includes the functions in 'sources' folder

(mapcar #'(lambda (file) (compile&load 
                        (make-pathname :directory (append (pathname-directory *load-pathname*) '("sources")) 
                                       :name file)))
      '(
          "unlisted-functions"
          "classification"
          "metrics"
          "list-processing"
          "midicents"
          "rhythm"
          "score"
          "misc"
        ))

(om::fill-library 
    '(
        ("1-Classification" nil nil (
            NNS
            K-means
            DTW
            DTW-align
            KDTree
            KNN
            Markov-build
            Markov-run
            Optimal-sorting
            ) nil)
        ("2-Metrics" nil nil (
            Euclidean-distance
            List-moments
            List-median
            List-mode
            List-Zscore
            List-covariance
            List-correlation
            Histogram
            Plot-points
            ) nil)
        ("3-List processing" nil nil (
            X-interpolation
            Posn-map
            Pick-random
            List-wrap
            List-fold
            List-symmetries
            List-frames
            List-quantize
            List-mod
            Nested-position
            Nested-nth
            Deep-nth
            Deep-replace
            N-occurances
            Unique-seq
            Unique-scramble
            ) nil)
        ("4-Midicents" nil nil (
            NRT
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
        ("7-Miscellaneous" nil nil (
            L-system
            2D-Turtle
            Vieru-seq
            Make-sieve
            ) nil)
    )
)

(print "
***************************
OM-Data library - 2021
[www.felipe-tovar-henao.com]
***************************
")
