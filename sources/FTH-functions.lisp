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

(defun closest (a b-list)
    (setq distances (loop for b in b-list and n from 0 to (- (length b-list) 1) collect 
        (list (abs (- b a)) b n)
    ))
    (stable-sort distances #'< :key #'first)
    (setq distances (car distances))
    (list (second distances) (third distances))
)

(defun depth (list)
    (if (listp list)
        (+ 1 (reduce #'max (mapcar #'depth list) :initial-value 0)) 0))

(defun random-list (dims minval maxval)
    (setq v-out nil)
    (dotimes (x dims)
        (setq v-out (append v-out (list (om-random minval maxval))))
    ) v-out)

(defun k-smart (data k)
    (setq k-centroids (list (nth-random data)))
    (setq data (remove (car k-centroids) data :test 'equal))
    (setq x 1)
    (while (< x k)
        (setq w-prob nil)
        (loop for d in data do
            (setq distances (loop for kc in k-centroids collect
                (expt (Euclidean-distance d kc nil) 2)))
            (setq w-prob (append w-prob (list (list-min distances))))
        )
        (setq new-ck (Nth-wrand data w-prob 1))
        (setq data (remove new-ck data :test 'equal))
        (setq k-centroids (append k-centroids new-ck))
        (setq x (+ x 1)))
    k-centroids)

(defun position-from-nested (a b-list)
    (if (atom a)
        (position a b-list)
        (if (eq (depth a) 1)
            (progn
                (loop for x in a collect
                    (position x b-list)
                )
            )
            (loop for x in a collect
                (position-from-nested x b-list)
            ))))

(defun nth-from-nested (a b-list)
    (if (atom a)
        (nth a b-list)
        (if (eq (depth a) 1)
            (loop for x in a collect
                    (nth x b-list))
            (loop for x in a collect
                (nth-from-nested x b-list)
            ))))

(defun mix (a b f)
    (if (and (atom a) (atom b))
        (+ (* a (- 1 f)) (* b f))
        (loop for x in a and y in b collect
            (mix x y f))))

(defun nested-mix (a b f)
    (if (or (eq (+ (depth a) (depth b)) 2) (eq (atom a) (atom b)))
        (mix a b f)
        (mapcar #'(lambda (in1 in2) (nested-mix in1 in2 f)) a b)))

; -------------- M E T H O D S ---------------------

; -------------- Distortion ---------------------
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
    (setq nns-list nil) 
    (setq positions nil)
    (setq distances 
        (loop for l in other-lists collect 
            (list l (Euclidean-distance main-list l weights))
        )
    )
    (stable-sort distances #'< :key #'second)
    (car (mat-trans distances))
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
        ; (setq remaining (cdr remaining))
        (setq remaining (NNS (car remaining) (cdr remaining) weights))
    )
    (append (list st-list) output)
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

; --------------- List-mean ---------------
(defmethod! List-mean ((l list))
    :initvals '((0 1 2 3))
    :indoc '("list")
    :icon 000
    :doc "Computes list-wise mean" 
    (if (> (depth l) 1)
        (loop for i in l collect
            (if (not (eq i nil))
                (List-mean (remove nil i))
            )
        )
        (/ (reduce #'+ l) (length l))
    )
)

; --------------- Nth-wrand ---------------
(defmethod! Nth-wrand ((data list) (weights list) (times integer))
    :initvals '(((1 1 1) (2 2 2) (3 3 3) (4 4 4)) (1 2 3 4) 1)
    :indoc '("list" "list" "integer")
    :icon 000
    :doc "Weighted nth random"
    (setq datasize (length data))
    (setq weights (om/ weights (list-max weights)))
    (setq weights (om-round (om-scale weights (* datasize 5.0) (list-min weights) (list-max weights))))
    (setq options (flat (loop for d in data and w in weights collect
        (repeat-n d w)) 1))
    (repeat-n (nth-random options) times) 
)

;--------------- K-means ---------------
(defmethod! K-means ((data list) (k integer) (weights list))
    :initvals '(((0 1 0) (-3 -1 2) (4 0 9) (-3 -5 -1) (0 4 -3) (2 1 -4)) 2 nil)
    :indoc '("list" "k (integer)" "weights (optional)")
    :icon 000
    :doc "Unsupervised data clustering algorithm.
    
    NOTE: All data items must have the same size. Weights are optional" 
    
    ; clip k to suitable range
    (setq k (clip k 1 (- (length data) 1))) 

    ; copy data and add label slots
    (setq labeled-data 
        (loop for d in data collect (list d nil)))

    ; initialize k-centroids and sort by first element
    (setq k-centroids (k-smart data k))
    (stable-sort k-centroids #'< :key #'first)

    ; convergence flag for loop
    (setq convergence-flag nil)

    ; ----- K_MEANS routine ----
    (while (eq convergence-flag nil)

        ; keep a history of last k-centroids
        (setq pk-centroids (copy-list k-centroids))

        ; assign a k-centroid to each data item based on proximity
        (loop for ld in labeled-data and ld-pos from 0 to (length labeled-data) do
            (setq nearest-k (car (NNS (car ld) k-centroids weights)))
            (setq ck (position nearest-k k-centroids :test 'equal))
            (setf (nth ld-pos labeled-data) (list (car ld) ck)))
        
        ; update k-centroids 
        (setq t-labeled-data (mat-trans labeled-data))
        (loop for ck from 0 to (- k 1) do
            (setq current-k nil)
            (loop for ld in labeled-data do
                (if (eq ck (second ld))
                    (setq current-k (append current-k (list (first ld))))
                )
            )
            (setq current-k (mat-trans current-k))
            (if (not (equal current-k nil))
                (progn (setq new-centroid (List-mean current-k))
                    (setf (nth ck k-centroids) new-centroid))
            )
        )
        ; stop loop if centroids do not change
        (setq convergence-flag (equal k-centroids pk-centroids))
    )

    ; group data by classes
    (stable-sort labeled-data #'< :key #'second)
    (setq output (loop for n from 0 to (- k 1) collect nil))
    (loop for ld in labeled-data do
        (setf (nth (second ld) output) (append (nth (second ld) output) (list (first ld)))))
    output)

;--------------- X-interpolation ---------------
(defmethod! X-interpolation ((a-list list) (b-list list) (traj list))
    :initvals '(
        ((7200 7700 8100) (6200 6700 7100) (7600 8100 7200) (5300 5900 5000) (7900 7200 7600) (5700 5000 5300) (7100 6400 6700) (6000 6500 6900)) 
        (4500 5700 6402 6900 7286 7602 7868 8100 8304 8486) 
        (0.0 1.0 0.0))
    :indoc '("list" "list" "list")
    :icon 000
    :doc "Cross-interpolation between elements of list A and elements of list B, defined by a given interpolation trajectory/path list. The trajectory list must have at least two values, all between 0.0 and 1.0, which represent the normalized percentage of linear interpolation from list A to list B.

    List A and B do not need to have the same length, but list-b must be of depth 1. List A can be a nested list"

    (setq numpts (length a-list))
    (setq traj (om/ traj (list-max (om-abs traj))))
    (setq traj (om-clip traj 0.0 1.0))
    (setq interp-pts (nth-value 2 (om-sample traj numpts)))
    (setq a-table (remove-dup (flat (copy-list a-list)) 'eq 1))
    (stable-sort a-table #'<)
    (setq ab-matrix nil)
    (loop for a in a-table do 
        (setq b-target (closest a b-list))
        (setq ab-matrix (append ab-matrix (list (list a (car b-target))))))
    (setq ab-positions (position-from-nested a-list (first (mat-trans ab-matrix))))
    (setq nested-target (nth-from-nested ab-positions (second (mat-trans ab-matrix))))
    (loop for a in a-list and nt in nested-target and f in interp-pts collect
        (nested-mix a nt f)))  

#| 
    TODO:
        - DTW
        - KDTree
        - Chroma count
        - moments (stdev mean mode median)
        - sort-data by

 |#