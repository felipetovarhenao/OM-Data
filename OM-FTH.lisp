; OM-FTH library
; list all the functions written in the 'sources' subfolder
(in-package :om) ; loads om library

(compile&load (om-relative-path '("sources") "FTH-functions")) ; includes the functions in 'sources' folder

(om::fill-library 
    '(
        ("1-List processing" nil nil (
            Euclidean-distance
            NNS
            K-means
            X-interpolation
            List-path
            List-quantize
            List-mod
            List-mean
            Posn-map
            Nth-wrand
            ) nil)
        ("2-Midicents" nil nil (
            Distortion
            Fill-range
            Shift-posn
            ) nil)
    )
) ; 