; OM-FTH library
; list all the functions written in the 'sources' subfolder
(in-package :om) ; loads om library

(compile&load (om-relative-path '("sources") "FTH-functions")) ; includes the functions in 'sources' folder

(om::fill-library 
    '(
        ("1-Basic functions" nil nil (
            euclidean-distance
            get-posn
            closest
            thin
            matrix-nth
            ) nil) 
        ("2-List processing" nil nil (
            Euc-distance
            NNS
            List-path
            List-quantize
            List-mod
            Posn-map
            ) nil)
        ("3-Midicents" nil nil (
            Distortion
            Fill-range
            Shift-posn
            ) nil)
    )
) ; 