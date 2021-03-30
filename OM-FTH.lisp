; OM-FTH library
; list all the functions written in the 'sources' subfolder
(in-package :om) ; loads om library

(compile&load (om-relative-path '("sources") "FTH-functions")) ; includes the functions in 'sources' folder

(om::fill-library 
    '(
        ("1-Analysis" nil nil (
            Euclidean-distance
            NNS
            K-means
            List-mean
            List-variance
            List-path
            ) nil)
        ("2-List processing" nil nil (
            X-interpolation
            List-quantize
            List-mod
            Posn-map
            Nth-wrand
            ) nil)
        ("3-Midicents" nil nil (
            Distortion
            Fill-range
            Shift-posn
            ) nil)
    )
) ; 