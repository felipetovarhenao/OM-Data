#|--------------------------------- * 
|     OM-Data library functions     |
|   [www.felipe-tovar-henao.com]    |
|               2021                |
* --------------------------------- *
|#

(in-package :om)

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
          "misc"))

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
            Rep-filter
            Unique-scramble
            ) nil)
        ("4-Midicents" nil nil (
            NRT
            Harmonic-distr
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
            Score-filter
            Segment-seq
            Extract-channel
            Multi-join
            Get-transients
            Granulate
            Harmonic-mapping
            ) nil)
        ("7-Miscellaneous" nil nil (
            L-system
            2D-Turtle
            Vieru-seq
            Make-sieve
            Recursivus
            ) nil)
    )
)

(doc-library "
OM-Data is a library of elementary functions, partially aimed at data modelling and analysis in OpenMusic. The library includes functions for data metrics, classification, and processing, as well as some higher-level functions for specific musical operations.
These include <i>K-means</i>, <i>DTW</i>, <i>NNS</i>, <i>KDTree</i>, <i>KNN</i>, <i>Markov-build</i>, <i>Markov-run</i>,  <i>Segment-seq</i>, <i>Score-filter</i>, <i>Get-transients</i>, <i>Chroma-count</i>, <i>IC-vector</i>, among many others.
<br><br>
A list of example patches are included, demonstrating possible musical applications for some of these functions.
<br><br>
Once the library is installed and loaded into the OM workspace, the example patches will be available in <i>/Help/Import Tutorial Patches/Libraries/OM-Data</i>
<br><br>
OM-Data is still in development.
" 
    (find-library "OM-Data"))

(om::set-lib-release 1.0)

(om-print "
***************************
OM-Data library v1.0
[www.felipe-tovar-henao.com]
2021
***************************
")