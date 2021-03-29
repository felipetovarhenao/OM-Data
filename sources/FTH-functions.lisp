; OM-FTH functions
; by Felipe Tovar-Henao

(in-package :om)

; -------------- FUNCTIONS -------------------

(defun get-posn (n source)
    (setq i 0)
    (setq out nil)
    (while (< i (length source))
        (if (eq (nth i source) n)
            (setq out (flat (list out i)))
        )
        (setq i (+ i 1))
    )
    (remove nil out)
)

(defun matrix-nth (matrix posn-list)
    (setq out (copy-list matrix))
    (loop for p in posn-list do
        (setq out (nth p out))
    )
    out
)

(defun ic-vector (mc-list &optional mod-n?)
    (if (equal mod-n? nil)
        (setq mod-n? 12)
    )
    (setq pc-set (loop for n in mc-list collect (mod n (* 100 mod-n?))))
    (stable-sort pc-set #'<)
    (setq i 0)
    (setq set-size (length pc-set))
    (setq out-vector '(0 0 0 0 0 0 0))
    (while (< i set-size)
        (setq pc-a (nth i pc-set))
        (setq j (+ 1 i))
        (while (< j set-size)
            (setq pc-b (nth j pc-set))
            (setq mc-i (mod (- pc-b pc-a) 1200))
            (if (> mc-i 600) 
                (setq mc-i (- 600 (mod mc-i 600)))
            )
            (setq v-pos (om-round (/ mc-i 100)))
            (setf (nth v-pos out-vector) (+ 1(nth v-pos out-vector)))
            (print out-vector)
            (setq j (+ j 1))
        )
        (setq i (+ i 1))
    )
    (cdr out-vector)
)

(defun depth (list)
    (if (listp list)
        (+ 1 (reduce #'max (mapcar #'depth list) :initial-value 0)) 0))

; -------------- METHODS ---------------------
(defmethod! Distortion ((mc-list list) (dist number))
    :initvals '((100 200 300 400 500) 1.125)
    :indoc '("midicent list" "distortion index")
    :icon 000
    :doc "Takes a list of midicents and applies spectral distortion (compression/expansion) to it, according to a distortion index."
    (setq fqlist (mc->f mc-list))
    (stable-sort fqlist #'<)
    (setq fq0 (car fqlist))
    (loop for fq in fqlist collect
        (* fq0 (expt (/ fq fq0) dist))
    )
)
; --------------- Euclidean-distance ---------------
(defmethod! Euclidean-distance ((a-list list) (b-list list) (weights list))
    :initvals '((0 1 2 3) (1 2 3 4) nil)
    :indoc '("list" "list of lists" "list (optional)")
    :icon 000
    :doc "Computes the Euclidean distance from main-list to other-lists

    NOTE: All lists must have the same length.
    "
    (if (equal weights nil)
        (setq weights (repeat-n 1.0 (length a-list)))
    )
    (setq weights (mapcar #'(lambda (input-list) (/ input-list (list-max weights))) weights))
    (setq l-depth (depth b-list))
    (cond 
        ((eq l-depth 1) (sqrt (reduce #'+ 
            (loop for a in a-list and b in b-list and w in weights collect 
                (* w(expt (- b a) 2))
            ))))
        ((> l-depth 1) 
            (loop for b in b-list collect
                (Euclidean-distance a-list b weights)
            ))))

; --------------- NNS ---------------
(defmethod! NNS ((main-list list) (other-lists list) (weights list))
    :initvals '((0 1 2 3) ((0 1 2 3) (1 2 3 4) (2 3 4 5)) nil)
    :indoc '("list" "list of lists" "list (optional)")
    :icon 000
    :doc "Sorts the lists based on the Nearest-Neighbor Search algorithm.
    
    NOTE: All lists must have the same length."
    (setq nns-list nil) ; create empty nns-list
    (setq positions nil) ; create empty positions list
    (setq distances 
        (loop for l in other-lists collect 
            (list l (Euclidean-distance main-list l weights))
        )
    )
    (stable-sort distances #'< :key #'second)
    (append (list main-list) (car (mat-trans distances)))
)

; --------------- List-path ---------------
(defmethod! List-path ((st-list list) (other-lists list) (weights list))
    :initvals '((0 1 2 3) ((0 1 2 3) (1 2 3 4) (2 3 4 5)) nil)
    :indoc '("list (initial)" "list of lists" "list (optional)")
    :icon 000
    :doc "Sorts all the lists such that the element-wise distance between resulting adjacent lists is minimal. Helpful for sorting chords (mc lists) for optimal voice leading.
    
    NOTE: All lists must have the same length."  
    (setq neighbors (NNS st-list other-lists weights))
    (setq remaining (copy-list neighbors))
    (setq output nil)
    (loop for n in neighbors do
        (setq output (append output (list (car remaining))))
        (setq remaining (cdr remaining))
        (setq remaining (NNS (car remaining) (cdr remaining) weights))
    )
    output
)

; --------------- List-quantize ---------------
(defmethod! List-quantize ((a number) (b-list list) (accuracy number))
    :initvals '(2.5 (0 1 2 3 4 5 6) 1)
    :indoc '("source (list)" "target (list)" "accuracy (optional)")
    :icon 000
    :doc "Approximates/quantizes the elements from the source list to the closest elements in the target list" 
    (if (equal accuracy nil)
        (setq accuracy 1)
    )
    (setq accuracy (clip accuracy 0 1))
    (setq distances nil)
    (loop for b in b-list collect
        (setq distances(append distances (list (abs (- b a)))))
    )
    (setq sorted-distances (copy-list distances))
    (stable-sort sorted-distances #'<)
    (+ (* accuracy (nth (nth 0 (get-posn (car sorted-distances) distances)) b-list)) (* a (- 1 accuracy)))
)

(defmethod! List-quantize ((a-list list) (b-list list) (accuracy number))
    (setq l-depth (depth a-list))
    (cond 
        ((eq l-depth 1)
            (mapcar #'(lambda (input) 
                (List-quantize input b-list accuracy)) a-list))
        ((> l-depth 1)
            (loop for a in a-list collect
                (List-quantize a b-list accuracy)
            ))))

; --------------- List-mod ---------------
(defmethod! List-mod ((input-list list) (n number))
    :initvals '((-3 -2 -1 0 1 2 3) 2)
    :indoc '("list" "mod-n")
    :icon 000
    :doc "Applies sign-preserving modulo arithmetic to input-list"  
    (setq output nil)
    (loop for i in input-list collect
        (setq output (append output (list
            (cond 
                ((< i 0) (* (mod (abs i) n) (/ i (abs i))))
                ((>= i 0) (mod i n))
            )
        )))
    )
    output)

(defmethod! List-mod ((input-list list) (n list))
    (mapcar (lambda (input) 
        (List-mod input-list input)) n))

; --------------- Fill-range ---------------
(defmethod! Fill-range ((chord-list list) (lower-bound number) (upper-bound number))
    :initvals '((3600 5200 6700 7000) 3600 9000)
    :indoc '("midicent list" "range lower bound" "range upperbound")
    :icon 000
    :doc "Fills the specified range with the midicents from chord-list" 
    (setq chord-list (remove-dup (List-mod chord-list 1200) #'eq 1))
    (stable-sort chord-list #'<)
    (setq base-chord (copy-list chord-list))
    (setq o 1)
    (setq offset (* 1200 (om-round (/ lower-bound 1200))))
    (setq numoctaves (+ (om-round (/ (abs (- upper-bound lower-bound)) 1200.0)) 1))
    (while (<= o numoctaves)
        (setq chord-list (append chord-list (om+ base-chord (+ (* o 1200) offset))))
        (setq o (+ o 1))
    )
    chord-list
    (band-filter chord-list (list (list lower-bound upper-bound)) 'pass)
)
; --------------- Shift-posn ---------------
(defmethod! Shift-posn ((chord-list list) (n-step number))
    :initvals '((3600 5200 6700 7000) 1)
    :indoc '("midicent list" "chordal step")
    :icon 000
    :doc "shifts a collection of midicents by n steps along itself, assuming octave equivalence between pitches." 
    (setq filled-range (Fill-range chord-list 0 20000))
    (loop for note in chord-list collect
        (nth (+ (position note filled-range) n-step) filled-range)
    )
)

(defmethod! Shift-posn ((chord-list list) (n-step list))
    (mapcar #'(lambda (input) 
        (Shift-posn chord-list input)) n-step))

; --------------- Posn-map ---------------
(defmethod! Posn-map ((l list))
    :initvals '(((0 0) (0 1) (1 0) (1 0) (0 0) (1 1)))
    :indoc '("list")
    :icon 000
    :doc "Returns a list of element positions, based on possible repetitions" 
    (setq thin-l nil)
    (setq out nil)
    (loop for x in l do
        (cond (
            (equal (member x thin-l :test 'equal) nil) 
            (setq thin-l (append thin-l (list x))))
        )
        (setq out (append out (list (position x thin-l :test 'equal))))
    ) 
    out)

