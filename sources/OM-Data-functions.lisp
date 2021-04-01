#|--------------------------------- * 
|     OM-Data library functions     |
|   [www.felipe-tovar-henao.com]    |
|               2021                |
* --------------------------------- *
|#

(in-package :om)

; -------------- (HIDDEN) FUNCTIONS -------------------
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
    :doc "Applies spectral distortion to a list of midicents, given a distortion index d. In other words, it expands or compresses the intervals in relation to the lowest midicent value in the input list."
    (if (eq (depth mc-list) 1)
        (progn
            (setq fqlist (mc->f mc-list))
            (stable-sort fqlist #'<)
            (setq fq0 (list-min fqlist))
            (setq output (loop for fq in fqlist collect
                (* fq0 (expt (/ fq fq0) dist))))
            (f->mc output))
        (loop for mc-l in mc-list collect
            (Distortion mc-l dist))))

; --------------- Euclidean-distance ---------------
(defmethod! Euclidean-distance ((a-list list) (b-list list) (weights list))
    :initvals '((0 1 2 3) (1 2 3 4) nil)
    :indoc '("list" "list of lists" "list (optional)")
    :icon 000
    :doc "Computes the Euclidean distance from one list to a list of lists.

    NOTE: All lists must have the same length."
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
                (Euclidean-distance a-list b weights)))))

; --------------- NNS ---------------
(defmethod! NNS ((main-list list) (other-lists list) (weights list))
    :initvals '((0 1 2 3) ((0 1 2 3) (1 2 3 4) (2 3 4 5)) nil)
    :indoc '("list" "list of lists" "list (optional)")
    :icon 000
    :doc "Sorts the lists based on the exhaustive nearest neighbor seach algorithm, using Euclidean distance as the sorting measurement.
    
    NOTE: All lists must have the same length."
    (setq nns-list nil) 
    (setq positions nil)
    (setq distances 
        (loop for l in other-lists collect 
            (list l (Euclidean-distance main-list l weights))
        )
    )
    (stable-sort distances #'< :key #'second)
    (car (mat-trans distances)))

; --------------- Optimal-sorting ---------------
(defmethod! Optimal-sorting ((st-list list) (other-lists list) (weights list))
    :initvals '((0 0 0 2) ((0 1 2 3) (1 2 3 4) (2 3 4 5)) nil)
    :indoc '("list (initial)" "list of lists" "list (optional)")
    :icon 000
    :doc "Sorts a list of lists such that the distance between adjacent lists is optimally minimized, given a starting list.
    
    NOTE: All lists must have the same length."  
    (setq neighbors (NNS st-list other-lists weights))
    (setq remaining (copy-list neighbors))
    (setq output nil)
    (loop for n in neighbors do
        (setq output (append output (list (car remaining))))
        (setq remaining (NNS (car remaining) (cdr remaining) weights))
    )
    (append (list st-list) output))

; --------------- List-quantize ---------------
(defmethod! List-quantize ((a number) (b-list list) (accuracy number))
    :initvals '(2.5 (0 1 2 3 4 5 6) nil)
    :indoc '("source (list)" "target (list)" "accuracy (optional)")
    :icon 000
    :doc "Approximates/quantizes the values from the source list to the closest elements in the target list, given an normalized accuracy percentage (optional) between 0.0 and 1.0" 
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
                (List-quantize a b-list accuracy)))))

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
                ((>= i 0) (mod i n)))))))
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
    (setq base-chord (copy-list chord-list))
    (setq o 0)
    (setq offset (* 1200 (values (floor (/ lower-bound 1200)))))
    (setq numoctaves (values (ceiling (/ (abs (- upper-bound lower-bound)) 1200.0))))
    (setq output nil)
    (while (<= o numoctaves)
        (setq output (append output (om+ base-chord (+ (* o 1200) offset))))
        (setq o (+ o 1))
    )
    (stable-sort output #'< )
    (band-filter output (list (list lower-bound upper-bound)) 'pass))

; --------------- Shift-posn ---------------
(defmethod! Shift-posn ((chord-list list) (n-step number))
    :initvals '((3600 5200 6700 7000) 1)
    :indoc '("midicent list" "chordal step")
    :icon 000
    :doc "shifts a collection of midicents by n steps along itself, assuming octave equivalence between pitches." 
    (setq filled-range (Fill-range chord-list 0 20000))
    (loop for note in chord-list collect
        (nth (+ (position note filled-range) n-step) filled-range)))

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

(defmethod! List-moments ((data list) (moments list))
    :initvals '((0 1 2 3) (0 1 2 3))
    :indoc '("list" "list")
    :icon 000
    :numouts 1
    :outdoc '("mean st-dev skewness kurtosis")
    :doc "Computes the statistical moments of a list of values: population mean, population standard deviation, skewness and kurtosis" 
    (if (eq (depth data) 1)
        (if (> (length data) 1)
            (progn 
                (setq mean (/ (reduce #'+ data) (* 1.0 (length data))))
                (setq st-dev (sqrt 
                    (/ (reduce #'+ 
                        (loop for x in data collect
                            (expt (- x mean) 2))) 
                        (length data))))
                (setq skewness
                    (/ (reduce #'+
                        (loop for x in data collect
                            (expt (- x mean) 3)
                        )) 
                        (* (expt st-dev 3) (length data) 1)))
                (setq kurtosis
                    (/ (reduce #'+
                        (loop for x in data collect
                            (expt (- x mean) 4))) 
                        (* (expt st-dev 4) (length data))))

                (posn-match (list mean st-dev skewness kurtosis) moments))
            (posn-match (list (car data) 0 nil nil) moments))
            (loop for d in data collect (List-moments d moments))))

;--------------- List-Zscore ---------------
(defmethod! List-Zscore ((x number) (l list))
    :initvals '((0 3 6) (0 1 2 3 4 5 6))
    :indoc '("list" "list")
    :icon 000
    :doc "Computes the standard score (a.k.a. z-score) value of a list"  
    (setq m-list (List-moments l '(0 1)))
    (setq mu (first m-list))
    (setq sigma (second m-list))
    (/ (- mu x) sigma))

(defmethod! List-Zscore ((x-list list) (l list))
    (mapcar #'(lambda (input) (List-Zscore input l)) x-list))

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
    (repeat-n (nth-random options) times))

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
                (progn 
                    (setq new-centroid (flat (List-moments current-k '(0))))
                    (setf (nth ck k-centroids) new-centroid))))

        ; stop loop if centroids do not change
        (setq convergence-flag (equal k-centroids pk-centroids)))

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
    (stable-sort b-list #'<)
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
    (setq ab-positions (Nested-position a-list (first (mat-trans ab-matrix))))
    (setq nested-target (Nested-nth ab-positions (second (mat-trans ab-matrix))))
    (loop for a in a-list and nt in nested-target and f in interp-pts collect
        (nested-mix a nt f)))  

; -------------- Nested-position ---------------------
(defmethod! Nested-position ((a list) (b-list list))
    :initvals '((("a") "b" (("c") "d") ((("e")))) ("a" "b" "c" "d" "e"))
    :indoc '("list" "integer")
    :icon 000
    :doc "Finds the position in list B for every element in list A. List B must be of depth 1"
    (if (eq (depth a) (depth b-list))
        (loop for x in a collect
            (position x b-list :test 'equal))
        (loop for x in a collect
            (Nested-position x b-list))))

(defmethod! Nested-position ((a number) (b-list list))
    (position a b-list :test 'equal))

(defmethod! Nested-position ((a string) (b-list list))
    (position a b-list :test 'equal))

(defmethod! Nested-position ((a symbol) (b-list list))
    (position a b-list :test 'equal))
    
; -------------- Nested-nth ---------------------
(defmethod! Nested-nth ((a list) (b-list list))
    :initvals '(((0) 1 ((2) 3) (((4)))) ("a" "b" "c" "d" "e"))
    :indoc '("list" "integer")
    :icon 000
    :doc "Finds the nth element in list B for every position in list A. List B must be of depth 1"
    (if (eq (depth a) 1)
        (loop for x in a collect
            (nth x b-list))
        (loop for x in a collect
            (Nested-nth x b-list))))

(defmethod! Nested-nth ((a number) (b-list list))
    (nth a b-list))

(defmethod! Nested-nth ((a string) (b-list list))
    (nth a b-list))

(defmethod! Nested-nth ((a symbol) (b-list list))
    (nth a b-list))

;--------------- List-wrap---------------
(defmethod! List-wrap ((val number) (lower-bound number) (upper-bound number))
    :initvals '(8 2 6)
    :indoc '("list" "number" "number")
    :icon 000
    :doc "Wraps the values of a list around a given range"
    (setq range (abs (- upper-bound lower-bound)))
    (setq out val)
    (cond
        (
            (< val lower-bound)
            (progn
                (setq diff (mod (abs (- lower-bound val)) range))
                (setq out (- upper-bound diff))))
        (
            (>= val upper-bound)
            (progn
                (setq diff (mod (abs (- upper-bound val)) range))
                (setq out (+ lower-bound diff)))))
    out
)

(defmethod! List-wrap ((val list) (lower-bound number) (upper-bound number))
    (mapcar #'(lambda (input) (List-wrap input lower-bound upper-bound)) val))

;--------------- List-fold---------------
(defmethod! List-fold ((val number) (lower-bound number) (upper-bound number))
    :initvals '(8 2 6)
    :indoc '("list" "number" "number")
    :icon 000
    :doc "Wraps the values of a list around a given range"
    (setq range (abs (- upper-bound lower-bound)))
    (setq out val)
    (cond
        (
            (< val lower-bound)
            (setq diff (abs (- lower-bound val)))
            (setq mode (evenp (nth-value 0 (om// diff range))))
            (setq diff (mod diff range))
            (if (eq mode t)
                (setq out (+ lower-bound diff))
                (setq out (- upper-bound diff))))
        (
            (>= val upper-bound)
            (setq diff (abs (- upper-bound val)))
            (setq mode (evenp (nth-value 0 (om// diff range))))
            (setq diff (mod diff range))
            (if (eq mode t)
                (setq out (- upper-bound diff))
                (setq out (+ lower-bound diff)))))
    out)

(defmethod! List-fold ((val list) (lower-bound number) (upper-bound number))
    (mapcar #'(lambda (input) (List-fold input lower-bound upper-bound)) val))

;--------------- Chroma-count ---------------
(defmethod! Chroma-count ((mc list))
    :initvals '((6000 6400 6700 7000 7200))
    :indoc '("list")
    :icon 000
    :doc "Recursively computes the chroma count of a list of midicents"
    (setq chroma-vector (repeat-n 0 12))
    (if (eq (depth mc) 1)
        (progn 
            (setq pc-list (om/ (List-mod mc 1200) 100)) 
            (loop for pc in pc-list do    
                (setq posn (nth-value 0 (round pc)))
                (setf (nth posn chroma-vector) (+ (nth posn chroma-vector) 1)))
            chroma-vector)
        (mapcar #'(lambda (input) (Chroma-count input)) mc)))

;--------------- IC-vector ---------------
(defmethod! IC-vector ((mc list))
    :initvals '((6000 6400 6700 7000 7200))
    :indoc '("list")
    :icon 000
    :doc "Recursively computes the interval vector of a list of midicents"
    (if (and (eq (depth mc) 1) (> (length mc) 1))
        (progn 
            (setq icv (repeat-n 0 7))
            (setq times 0)
            (setq pc-list (om-round (om/ (List-mod mc 1200) 100) 1))
            (loop for pc in pc-list do   
                (loop for x from times to (- (length mc) 1) do
                    (setq interval (abs (- (nth x pc-list) pc)))
                    (setq posn (List-fold (nth-value 0 (round interval)) 0 6))
                    (setf (nth posn icv) (+ (nth posn icv) 1)))
                (setq times (+ times 1)))
            (cdr icv))
        (mapcar #'(lambda (input) (IC-vector input)) mc)))

(defun dtw-path-cost (a-list b-list)
    ; build cost matrix
    (setq matrix (make-array (list (length a-list) (length b-list))))
    (loop for a in a-list and x from 0 to (length a-list) do
        (loop for b in b-list and y from 0 to (length b-list) do
            (if (> x 0) 
                (setq leftval (aref matrix (- x 1) y))
                (setq leftval nil))
            (if (> y 0) 
                (setq topval (aref matrix x (- y 1)))
                (setq topval nil)) 
            (if (and (> y 0) (> x 0)) 
                (setq diagval (aref matrix (- x 1) (- y 1)))
                (setq diagval nil))
            (setq options (remove nil (list leftval topval diagval)))
            (if (not (eq options nil))
                (setq cost (list-min options))
                (setq cost 0))
            (setf (aref matrix x y) (+ (abs (- a b)) cost))))

    (setq x (- (length a-list) 1)) 
    (setq y (- (length b-list) 1))
    (setq path-posn (list (list x y)))
    (setq path-cost (aref matrix x y))
    (while (or (> x 0) (> y 0))
         (if (> x 0) 
                (setq leftval (aref matrix (- x 1) y))
                (setq leftval nil))
            (if (> y 0) 
                (setq topval (aref matrix x (- y 1)))
                (setq topval nil)) 
            (if (and (> y 0) (> x 0)) 
                (setq diagval (aref matrix (- x 1) (- y 1)))
                (setq diagval nil))
            (setq options (list leftval topval diagval))
            (setq minval (list-min options))
            (cond
                (
                    (eq minval leftval)
                    (setq x (- x 1)))
                (
                    (eq minval topval)
                    (setq y (- y 1)))
                (
                    (eq minval diagval)
                    (progn 
                        (setq x (- x 1))
                        (setq y (- y 1)))))
            (setq path-cost (+ path-cost minval))
            (setq path-posn (append path-posn (list (list x y)))))
    (list path-cost (reverse path-posn))
    ; path-cost
)

; (defmethod! DTW ((list-a list) (list-b list))
;     :initvals '((0 1 2 3 4 5) ((1 1 2 3 5) (0 0 1 2 3 3 4 5) (1 3 4 5)))
;     :indoc '("list" "list of lists")
;     :icon 000
;     :doc "Sorts the lists in second input using Dynamic Time Warping"
;     (cond
;         (
;             (eq (depth list-b) 1)
;             (dtw-path-cost list-a list-b)
;         )
;         (
;             (eq (depth list-b) 2)
;             (setq dtw-list 
;                 (loop for b in list-b collect
;                     (list (dtw-path-cost list-a b) b)))
;             (stable-sort dtw-list #'< :key #'first)
;             (second (mat-trans dtw-list)))))

(defmethod! DTW ((list-a list) (list-b list))
    :initvals '((0 1 2 3 4 5) ((1 1 2 3 5) (0 0 1 2 3 3 4 5) (1 3 4 5)))
    :indoc '("list" "list of lists")
    :icon 000
    :numouts 2
    :doc "Sorts the lists in second input using Dynamic Time Warping"
    (cond
        (
            (eq (depth list-b) 1)
            (progn
                (setq out (dtw-path-cost list-a list-b))
                (values-list (list (first out) (second out)))
            )
        )
        (
            (eq (depth list-b) 2)
            (setq dtw-list nil)
                (loop for b in list-b do
                    (setq out (list (append (dtw-path-cost list-a b) (list b))))
                    (setq dtw-list (append dtw-list out)))
            (stable-sort dtw-list #'< :key #'first)
            (setq out (mat-trans dtw-list))
            (values-list (list (third out) (second out)))
        )
    )
)

(defmethod! DTW-align ((a list) (b list) (posn list))
    :initvals '((0 1 2 3 4 2) (0 2 3 4 1) ((0 0) (1 0) (2 1) (3 2) (4 3) (5 4)))
    :indoc '("list" "list" "list")
    :icon 000
    :doc "Works in combination with DTW, aligning the values of the main time series to the others. Use the left input of DTW as in1, and the outputs of DTW as in2 and in3, respectively."
    (if (and (eq (depth b) 1) (eq (depth posn) 2))
        (loop for p in posn collect
            (list (nth (first p) a) (nth (second p) b)))
        (mapcar #'(lambda (input1 input2) (DTW-align a input1 input2)) b posn)))
#| 
    TODO:
        - DTW
        - KDTree
        - Fundamental frequency estimator
        - Inharmonicity
        - sort-data by

 |#