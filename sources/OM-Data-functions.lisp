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
        (if (equal (nth i source) n)
            (setq out (flat (list out i))))
        (setq i (+ i 1)))
    (remove nil out))

(defun closest (a b-list)
    (setq distances (loop for b in b-list and n from 0 to (- (length b-list) 1) collect 
        (list (abs (- b a)) b n)))
    (stable-sort distances #'< :key #'first)
    (setq distances (car distances))
    ; list of element and position
    (list (second distances) (third distances)))

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
        (setq new-ck (Pick-random data w-prob 1))
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

;--------------- dtw-path-cost ---------------
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
            (setf (aref matrix x y) (+ (Euclidean-distance a b nil) cost))))
    
    ; find optimal path
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
    (list path-cost (reverse path-posn)))

; -------------- M E T H O D S ---------------------
#| 
; -------------- Distortion ---------------------
(defmethod! Distortion ((mc-list list) (dist number))
    :initvals '((100 200 300 400 500) 1.125)
    :indoc '("midicent list" "distortion index")
    :icon 000
    :doc "Applies spectral distortion to a list of midicents, given a distortion index d. In other words, it expands or compresses the intervals in relation to the lowest midicent value in the input list.
    
    Example:  
    (distortion '(100 200 300 400 500) 1.125) => (100 212 326 438 550)
    "
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
 |#

; -------------- Distortion ---------------------
(defmethod! Distortion ((mc-list list) (dist number) &optional (mc-fund nil))
    :initvals '((6000 6400 6700 7200) 1.125 nil)
    :indoc '("list or score sequence" "distortion index" "integer")
    :icon 000
    :doc "Applies spectral distortion to a list of midicents, given a distortion index d. In other words, it expands or compresses the intervals in relation to the lowest midicent value in the input list.
    
    Examples:  
    (distortion '(6000 6400 6700 7200) 1.125) => (6000 6450 6788 7350)
    (distortion '(6000 6400 6700 7200) 1.125 6000) => (6000 6450 6788 7350)
    (distortion '(6000 6400 6700 7200) 1.125 7200) => (5850 6300 6638 7200)
    "
    (if (equal mc-fund nil)
        (progn 
            (setq mc-fund (list-min (flat mc-list)))
            (setq fq0 (mc->f mc-fund)))
        (setq fq0 (mc->f mc-fund)))
    (if (eq (depth mc-list) 1)
        (progn
            (setq fqlist (mc->f mc-list))
            ; (stable-sort fqlist #'<)
            (if (equal fq0 nil)
                (setq fq0 (list-min fqlist))
            )
            (setq output (loop for fq in fqlist collect
                (* fq0 (expt (/ fq fq0) dist))))
            (f->mc output))
        (loop for mc-l in mc-list collect
            (Distortion mc-l dist mc-fund))))

(defmethod! Distortion ((self chord-seq) (dist number) &optional (mc-fund nil))
    (make-instance 'chord-seq :lmidic (Distortion (lmidic self) dist mc-fund) :lonset (lonset self) :ldur (ldur self) :lvel (lvel self) :loffset (loffset self) :lchan (lchan self)))

(defmethod! Distortion ((m-seq multi-seq) (dist number) &optional (mc-fund nil))
    (make-instance 'multi-seq :chord-seqs (loop for self in (chord-seqs m-seq) collect
        (make-instance 'chord-seq :lmidic (Distortion (lmidic self) dist mc-fund) :lonset (lonset self) :ldur (ldur self) :lvel (lvel self) :loffset (loffset self) :lchan (lchan self))))) 

(defmethod! Distortion ((self voice) (dist number) &optional (mc-fund nil))
    (setq new-chords (loop for ch in (chords self) collect 
        (make-instance 'chord :lmidic (Distortion (lmidic ch) dist mc-fund) :lvel (lvel ch) :loffset (loffset ch) :ldur (ldur ch) :lchan (lchan ch))))
    (make-instance 'voice :tree (tree self) :chords new-chords :tempo (tempo self) :legato (legato self) :ties (ties self)))

(defmethod! Distortion ((self poly) (dist number) &optional (mc-fund nil))
    (setq voice-list 
        (loop for v in (voices self) collect
            (Distortion v dist mc-fund)))
    (make-instance 'poly :voices voice-list))

; --------------- Euclidean-distance ---------------
(defmethod! Euclidean-distance ((a-list list) (b-list list) (weights list))
    :initvals '((0 1 2 3) ((1 2 3 4) (3 4 6 8)) nil)
    :indoc '("list" "list of lists" "list (optional)")
    :icon 000
    :doc "Computes the Euclidean distance from one list to a list of lists.

    NOTE: All lists must have the same length.
    
    Example:  
    (euclidean-distance '(0 1 2 3) '((1 2 3 4) (3 4 6 8)) nil)) => (2.0 7.6811457)
    "
    (if (equal weights nil)
        (setq weights (repeat-n 1.0 (length a-list))))
    (setq weights (mapcar #'(lambda (input-list) (/ input-list (list-max weights))) weights))
    (setq l-depth (depth b-list))
    (cond 
        ((eq l-depth 1) (sqrt (reduce #'+ 
            (loop for a in a-list and b in b-list and w in weights collect 
                (* w(expt (- b a) 2))))))
        ((> l-depth 1) 
            (loop for b in b-list collect
                (Euclidean-distance a-list b weights)))))

(defmethod! Euclidean-distance ((a number) (b number) (weights list))
    (abs (- b a)))

(defmethod! Euclidean-distance ((a number) (b-list list) (weights list))
    (loop for b in b-list collect (abs (- b a))))

; --------------- NNS ---------------
(defmethod! NNS ((main-list list) (other-lists list) (weights list))
    :initvals '((0 1 2 3) ((1 2 3 4) (0 0 2 3) (2 3 4 5)) nil)
    :indoc '("list" "list of lists" "list (optional)")
    :icon 000
    :doc "Sorts the lists based on the exhaustive nearest neighbor seach algorithm, using Euclidean distance as the sorting measurement.
    
    NOTE: All lists must have the same length.
    Example:  
    (NNS '(0 1 2 3) '((1 2 3 4) (0 0 2 3) (2 3 4 5)) nil) => ((0 0 2 3) (1 2 3 4) (2 3 4 5))
    "
    (setq nns-list nil) 
    (setq positions nil)
    (setq distances 
        (loop for l in other-lists collect 
            (list l (Euclidean-distance main-list l weights))))
    (stable-sort distances #'< :key #'second)
    (car (mat-trans distances)))

(defmethod! NNS ((a number) (list-b list) (weights list))
    (setq distances (loop for b in list-b collect (list b (Euclidean-distance a b nil))))
    (stable-sort distances #'< :key #'second)
    (car (mat-trans distances)))

; --------------- Optimal-sorting ---------------
(defmethod! Optimal-sorting ((st-list list) (other-lists list) (weights list))
    :initvals '((0 0 0 2) ((0 1 2 3) (2 3 4 5) (1 2 3 4)) nil)
    :indoc '("list (initial)" "list of lists" "list (optional)")
    :icon 000
    :doc "Sorts a list of lists such that the distance between adjacent lists is optimally minimized, given a starting list.
    
    NOTE: All lists must have the same length.
    
    Example:
    (optimal-sorting '(0 0 0 2) '((0 1 2 3) (2 3 4 5) (1 2 3 4)) nil) => ((0 0 0 2) (0 1 2 3) (1 2 3 4) (2 3 4 5))
    "  
    (setq neighbors (NNS st-list other-lists weights))
    (setq remaining (copy-list neighbors))
    (setq output nil)
    (loop for n in neighbors do
        (setq output (append output (list (car remaining))))
        (setq remaining (NNS (car remaining) (cdr remaining) weights)))
    (append (list st-list) output))

; --------------- List-quantize ---------------
(defmethod! List-quantize ((a number) (b-list list) (accuracy number))
    :initvals '((2.5 4.03) (0 1 2 3 4 5 6) 1.0)
    :indoc '("source (list)" "target (list)" "accuracy (optional)")
    :icon 000
    :doc "Approximates/quantizes the values from the source list to the closest elements in the target list, given an normalized accuracy percentage (optional) between 0.0 and 1.0.
    
    Example:
    (list-quantize '(2.5 4.03) '(0 1 2 3 4 5 6) 0.5) => (2.25 4.0150004)
    (list-quantize '(2.5 4.03) '(0 1 2 3 4 5 6) 1.0) => (2.0 4.0)
    " 
    (if (equal accuracy nil)
        (setq accuracy 1.0))
    (setq accuracy (clip accuracy 0 1))
    (setq distances nil)
    (loop for b in b-list collect
        (setq distances (append distances (list (abs (- b a))))))
    (setq sorted-distances (copy-list distances))
    (stable-sort sorted-distances #'<)
    (+ (* accuracy (nth (nth 0 (get-posn (car sorted-distances) distances)) b-list)) (* a (- 1 accuracy)))
    )

(defmethod! List-quantize ((a-list list) (b-list list) (accuracy number))
    (mapcar #'(lambda (input) 
        (List-quantize input b-list accuracy)) a-list))

; --------------- List-mod ---------------
(defmethod! List-mod ((input-list list) (n number))
    :initvals '((-3 -2 -1 0 1 2 3) 2)
    :indoc '("list" "mod-n")
    :icon 000
    :doc "Applies sign-preserving modulo arithmetic to input-list.
    
    Example:
    (list-mod '(-3 -2 -1 0 1 2 3) 2) => (-1 0 -1 0 1 0 1)
    "  
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
    :doc "Fills the specified range with the midicents from chord-list.
    
    Example:
    (fill-range '(3600 5200 6700 7000) 6000 7200) => (6000 6400 6700 7000 7200)
    " 
    (setq chord-list (remove-dup (List-mod chord-list 1200) #'eq 1))
    (setq base-chord (copy-list chord-list))
    (setq o 0)
    (setq offset (* 1200 (values (floor (/ lower-bound 1200)))))
    (setq numoctaves (values (ceiling (/ (abs (- upper-bound lower-bound)) 1200.0))))
    (setq output nil)
    (while (<= o numoctaves)
        (setq output (append output (om+ base-chord (+ (* o 1200) offset))))
        (setq o (+ o 1)))
    (stable-sort output #'< )
    (band-filter output (list (list lower-bound upper-bound)) 'pass))

; --------------- Shift-posn ---------------
(defmethod! Shift-posn ((chord-list list) (n-step number))
    :initvals '((3600 5200 6700 7000) 1)
    :indoc '("midicent list" "chordal step")
    :icon 000
    :doc "Shifts a collection of midicents by n steps along itself, assuming octave equivalence between pitches.
    
    Example:
    (shift-posn '(3600 5200 6700 7000) '(1 2 3)) => ((4000 5500 7000 7200) (4300 5800 7200 7600) (4600 6000 7600 7900))
    " 
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
    :doc "Returns a list of element positions, based on possible repetitions.
    
    Example:
    (posn-map '((0 0) (0 1) (1 0) (1 0) (0 0) (1 1))) => (0 1 2 2 0 3)
    " 
    (setq thin-l nil)
    (setq out nil)
    (loop for x in l do
        (cond (
            (equal (member x thin-l :test 'equal) nil) 
            (setq thin-l (append thin-l (list x)))))
        (setq out (append out (list (position x thin-l :test 'equal))))) 
    out)

;--------------- List-moments ---------------
(defmethod! List-moments ((data list) (moments list))
    :initvals '((5 2 3 4 5 6) (0 1 2 3))
    :indoc '("list" "list")
    :icon 000
    :numouts 1
    :outdoc '("mean st-dev skewness kurtosis")
    :doc "Computes the statistical moments of a list of values: population mean, population standard deviation, skewness and kurtosis. 
    The right inlet specifies which moments to put out, by index. For instance, (0 2) outputs the mean and skewness only. The list (0 1 2 3) outputs all four moments.
    
    Example:
    (list-moments '(5 2 3 4 5 6) '(0 1 2 3)) => (4.1666665 1.3437096 -0.30531597 1.8482844)
    " 
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
                            (expt (- x mean) 3))) 
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
    :doc "Computes the standard score (a.k.a. z-score) value of a list based on another.

    Example: 
    (list-zscore '(0 3 6) '(0 1 2 3 4 5 6)) => (1.5 0.0 -1.5)
    "  
    (setq m-list (List-moments l '(0 1)))
    (setq mu (first m-list))
    (setq sigma (second m-list))
    (/ (- mu x) sigma))

(defmethod! List-Zscore ((x-list list) (l list))
    (mapcar #'(lambda (input) (List-Zscore input l)) x-list))

; --------------- Pick-random ---------------
(defmethod! Pick-random ((in-list list) (weights list) (times integer))
    :initvals '(((1 1 1) (2 2 2) (3 3 3) (4 4 4)) (1 2 3 4) 1)
    :indoc '("list" "list" "integer")
    :icon 000
    :doc "Chooses a random element with weighted probabilities, N number of times.
    
    Example:
    (pick-random '((1 1 1) (2 2 2) (3 3 3) (4 4 4)) '(1 2 3 4) 3) => ((4 4 4) (4 4 4) (4 4 4))
    "
    (if (equal weights nil)
        (setq weights (repeat-n 1 (length in-list)))
    )
    (setq sum-w (reduce-tree weights #'+))
    (setq out nil)

    (loop for n from 1 to times do
        (setq rand (random (* 1.0 sum-w)))
        (setq stop nil)
        (setq index 0)
        (while (and (eq stop nil) (< index (length weights)))
            (setq w (nth index weights))
            (if (< rand w)
                (progn 
                    (setq out (append out (list (nth index in-list))))
                    (setq stop t)))
            (setq rand (- rand w))
            (setq index (+ index 1))))
    out)

;--------------- K-means ---------------
(defmethod! K-means ((data list) (k integer) (weights list))
    :initvals '(((0 1 0) (-3 -1 2) (4 0 9) (-3 -5 -1) (0 4 -3) (2 1 -4)) 2 nil)
    :indoc '("list" "k (integer)" "weights (optional)")
    :icon 000
    :doc "Unsupervised data clustering algorithm.
    
    NOTE: All data items must have the same size. Weights are optional.
    
    Example:
    (k-means '((0 1 0) (-3 -1 2) (4 0 9) (-3 -5 -1) (0 4 -3) (2 1 -4)) 2 nil) => (((0 1 0) (-3 -1 2) (-3 -5 -1) (0 4 -3) (2 1 -4)) ((4 0 9)))
    " 
    
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
                    (setq current-k (append current-k (list (first ld))))))
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

    List A and B do not need to have the same length, but list-b must be of depth 1. List A can be a nested list.
    "
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
    :doc "Finds the position in list B for every element in list A. List B must be of depth 1.
    
    Example: 
    (nested-position '((\"a\") \"b\" ((\"c\") \"d\") (((\"e\")))) '(\"a\" \"b\" \"c\" \"d\" \"e\")) => ((0) 1 ((2) 3) (((4))))
    "
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
    :doc "Finds the nth element in list B for every position in list A. List B must be of depth 1.
    
    Example:
    (nested-nth '((0) 1 ((2) 3) (((4)))) '(\"a\" \"b\" \"c\" \"d\" \"e\")) => ((\"a\") \"b\" ((\"c\") \"d\") (((\"e\"))))
    "
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
    :initvals '((-2 -1 0 1 2 3 4 5 6) 0 3)
    :indoc '("list" "number" "number")
    :icon 000
    :doc "Wraps the values of a list around a given range.
    
    Example:
    (list-wrap '(-2 -1 0 1 2 3 4 5 6) 0 3) => (1 2 0 1 2 0 1 2 0)
    "
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
    out)

(defmethod! List-wrap ((val list) (lower-bound number) (upper-bound number))
    (mapcar #'(lambda (input) (List-wrap input lower-bound upper-bound)) val))

;--------------- List-fold---------------
(defmethod! List-fold ((val number) (lower-bound number) (upper-bound number))
    :initvals '(8 2 6)
    :indoc '("list" "number" "number")
    :icon 000
    :doc "Folds the values of a list around a given range.
    
    Example:
    (list-fold '(-2 -1 0 1 2 3 4 5 6) 0 3) => (2 1 0 1 2 3 2 1 0)
    "
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
    :doc "Recursively computes the chroma count of a list of midicents.
    
    Example: 
    (chroma-count '(6000 6400 6700 7000 7200)) => (2 0 0 0 1 0 0 1 0 0 1 0)
    "
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
    :doc "Recursively computes the interval vector of a list of midicents.
    
    Example:
    (ic-vector '(6000 6400 6700 7000 7200)) => (0 2 2 2 2 1)
    "
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

;--------------- DTW ---------------
(defmethod! DTW ((list-a list) (list-b list))
    :initvals '((0 1 2 3 4 5) ((1 1 2 3 5) (0 0 1 2 3 3 4 5) (1 3 4 5)))
    :indoc '("list" "list of lists")
    :icon 000
    :numouts 2
    :outdoc '("sorted lists" "sorted")
    :doc "Sorts the lists in second input using Dynamic Time Warping.
    
    Example:
    (dtw '(0 1 2 3 4 5) '(1 1 2 3 5)) => [ 9 ((0 0) (1 0) (1 1) (2 2) (3 3) (4 4) (5 4)) ]
    "
    (cond
        (
            (eq (depth list-a) (depth list-b))
            (progn
                (setq out (dtw-path-cost list-a list-b))
                (values-list (list (first out) (second out)))))
        (
            (< (depth list-a) (depth list-b))
            (setq dtw-list nil)
            (loop for b in list-b do
                (setq out (list (append (dtw-path-cost list-a b) (list b))))
                (setq dtw-list (append dtw-list out)))
            (stable-sort dtw-list #'< :key #'first)
            (setq out (mat-trans dtw-list))
            (values-list (list (third out) (second out))))))

;--------------- DTW-align ---------------
(defmethod! DTW-align ((a list) (b list) (posn list))
    :initvals '((0 1 2 3 4 2) (0 2 3 4 1) ((0 0) (1 0) (2 1) (3 2) (4 3) (5 4)))
    :indoc '("list" "list" "list")
    :icon 000
    :doc "Works in combination with DTW, aligning the values of the main time series to the others. Use the left input of DTW as in1, and the outputs of DTW as in2 and in3, respectively.
    
    Example: 
    (dtw-align '(0 1 2 3 4 2) '(0 2 3 4 1) '((0 0) (1 0) (2 1) (3 2) (4 3) (5 4))) => ((0 0) (1 0) (2 2) (3 3) (4 4) (2 1))
    "
    (if (eq (depth a) (depth b))
        (loop for p in posn collect
            (list (nth (first p) a) (nth (second p) b)))
        (mapcar #'(lambda (input1 input2) (DTW-align a input1 input2)) b posn)))

;--------------- N-occurances ---------------
(defmethod! N-occurances ((x number) (l list))
    :initvals '(2 (0 1 2 1 4 3 2 2))
    :indoc '("item" "list")
    :icon 000
    :doc "Counts the number of occurances of an element in a list.

    Example:
    (n-occur 2 '(0 1 2 1 4 3 2 2)) => 3
    "
    (setq counter 0)
    (loop for y in l do
        (if (numberp y)
            (if (eq x y) 
                (setq counter (+ counter 1)))
            (setq counter (+ counter (N-occurances x y)))))
    counter)

(defmethod! N-occurances ((x list) (l list))
    (setq counter 0)
    (loop for y in l do
        (cond 
            (
                (eq (depth y) (depth x))
                (if 
                    (equal x y) 
                    (setq counter (+ counter 1))))
            (
                (> (depth y) (depth x))
                (setq counter (+ counter (N-occurances x y))))))
    counter)

;--------------- Unique-seq ---------------
(defmethod! Unique-seq ((l list))
    :initvals '(((1 2) (3 3) (4 2) (4 2) (2 2) (2 2)))
    :indoc '("list")
    :icon 000
    :doc "Removes repetitions between consecutive elements.

    Example:
    (unique-seq '((1 2) (3 3) (4 2) (4 2) (2 2) (2 2))) => ((1 2) (3 3) (4 2) (2 2))
    "
    (setq output (list (car l)))
    (loop for x from 1 to (- (length l) 1) do
        (if (not (equal (nth x l) (nth (- x 1) l)))
            (setq output (append output (list (nth x l))))))
    output)

;--------------- Mc-clip ---------------
(defmethod! Mc-clip ((mc number) (lower-bound number) (upper-bound number))
    :initvals '(5500 6000 7200)
    :indoc '("list" "number" "number")
    :icon 000
    :doc "Constrains a list of of midicents to a given range, using octave equivalence.
    
    Example: 
    (mc-clip 5500 6000 7200) => 6700
    "
    (setq lowdif (- 1200 (mod (abs (- lower-bound mc)) 1200)))
    (setq hidif (- 1200 (mod (abs (- upper-bound mc)) 1200)))
    (setq out mc)
    (cond
        (
            (< mc lower-bound)
            (setq out (+ lower-bound lowdif)))
        (
            (>= mc upper-bound)
            (setq out (- upper-bound hidif))))
    out)

(defmethod! Mc-clip ((mc-list list) (lower-bound number) (upper-bound number))
    (mapcar #'(lambda (input) (Mc-clip input lower-bound upper-bound)) mc-list))

;--------------- Mc-wrap ---------------
; (defmethod! Mc-wrap ((mc number) (lower-bound number) (upper-bound number))
;     :initvals '(5500 6000 7200)
;     :indoc '("list" "number" "number")
;     :icon 000
;     :doc ""
;     (setq range (abs (- upper-bound lower-bound)))
;     (setq numoct (* 1200 (max 1 (nth-value 0 (ceiling range 1200)))))
;     (setq lowdif (mod (abs (- mc lower-bound)) numoct))
;     (setq hidif (mod (abs (- mc upper-bound)) numoct))
;     (setq out mc)
;     (
;         cond
;         (
;             (< mc lower-bound)
;             (progn 
;                 (setq out (- upper-bound hidif))
;                 (if    
;                     (> hidif range)
;                     (setq out (+ out 1200)))))
;         (
;             (>= mc upper-bound)
;             (progn 
;                 (setq out (+ lower-bound lowdif))
;                 (if    
;                     (> lowdif range)
;                     (setq out (- out 1200))))))
;     out)

;--------------- Mc-wrap ---------------
(defmethod! Mc-wrap ((mc number) (lower-bound number) (upper-bound number))
    :initvals '(5500 6000 7200)
    :indoc '("number or list" "number" "number")
    :icon 000
    :doc "Wraps a list of of midicents around a given range, using octave equivalence.
    
    Example: 
    (mc-wrap 5900 6000 7200) => 7100
    "
    (setq range (abs (- upper-bound lower-bound)))
    (setq oct-offset (* 1200 (nth-value 0 (floor lower-bound 1200))))
    (setq oct-range (* 1200 (nth-value 0 (ceiling range 1200))))
    (setq oct-range-floor (* 1200 (nth-value 0 (floor range 1200))))
    (setq mc (+ oct-offset (abs (nth-value 1 (om// (- mc oct-offset) oct-range)))))
    (cond 
        (  
            (< mc lower-bound)
            (setq mc (+ mc oct-range-floor)))
        (
            (> mc upper-bound)
            (setq mc (- mc oct-range-floor))))
    mc)

(defmethod! Mc-wrap ((mc-list list) (lower-bound number) (upper-bound number))
    (mapcar #'(lambda (input) (Mc-wrap input lower-bound upper-bound)) mc-list))

;--------------- Ic-cycle ---------------
(defmethod! Ic-cycle ((mc-start number) (ic number) (lower-bound number) (upper-bound number) (n-times number))
    :initvals '(6300 (500 200) 6000 6700 5)
    :indoc '("number" "number or list" "number" "number" "number")
    :icon 000
    :doc "Builds an interval cycle given a starting pitch, a interval class, a lower and upper bound (range), and a number of iterations.

    Example:
    (ic-cycle 6300 '(500 200) 6000 6700 5) => (6300 6800 7000 6300 6500)
    "
    (setq out nil)  
    (loop for x from 0 to (- n-times 1) collect
        (setq out (append out (list (+ mc-start (* x ic)))))
    )
    (mc-wrap out lower-bound upper-bound))

(defmethod! Ic-cycle ((mc-start number) (ic list) (lower-bound number) (upper-bound number) (n-times number))
    (setq out (list mc-start))  
    (loop for x from 0 to (- n-times 2) collect
        (if (eq (depth ic) 1)
            (progn 
                (setq mc-start (+ mc-start (nth (mod x (length ic)) ic)))
                (setq out (append out (list mc-start))))))
    (mc-wrap out lower-bound upper-bound))

;--------------- Unique-scramble ---------------
(defmethod! Unique-scramble ((a-list list) (times integer))
    :initvals '((0 1 2) 4)
    :indoc '("list" "integer")
    :icon 000
    :doc "Performs a series of random permutations such that no element appears consecutively in the same position.

    Example:
    (unique-scramble '(0 1 2) 4) => ((0 1 2) (2 0 1) (1 2 0) (0 1 2))
    "
    (setq out (list a-list))
    (setq current (copy-list a-list))
    (loop for i from 1 to (- times 1) do
        (setq unique nil)
        (while (eq unique nil)
            (setq dups 0)
            (setq scrambled (permut-random current))
            (loop for c in current and s in scrambled do
                (if (eq c s)
                    (setq dups (+ dups 1))))
            (if (eq dups 0)
                (setq unique t)))
        (setq out (append out (list scrambled)))
        (setq current (copy-list scrambled)))
    out)

;--------------- Euclidean-rhythm ---------------
(defmethod! Euclidean-rhythm ((numbeats integer) (period integer) (rotation integer))
    :initvals '(5 13 0)
    :indoc '("integer" "integer" "integer")
    :icon 000
    :doc "euclidean-rhythm, as the name suggests, computes the Euclidean rhythm for a given period k and number of beats/pulses n, where n < k. The output can be interpreted as a list of intervals or distances between between pulses
    "
    (setq a numbeats)
    (setq b (- period numbeats))
    (setq binary-seq (list (repeat-n 1 a) (repeat-n 0 b)))
    (while (> (min a b) 1)
        (setq binary-seq (mapcar #'(lambda (input) (remove nil (flat input))) (mat-trans binary-seq)))
        (setq binary-seq (group-list binary-seq (sort-list (list a b)) 'linear))
        (setq a (length (car binary-seq)))
        (setq b (length (second binary-seq))))
    (setq binary-seq (flat binary-seq))
    (setq val 1)
    (setq out nil)
    (loop for n in binary-seq do
        (if (eq n 0)
            (setq val (+ val 1))
            (progn
                (setq out (append out (list val)))
                (setq val 1))))
    (rotate (cdr (append out (list val))) rotation))

;--------------- Rhythmicon ---------------
(defmethod! Rhythmicon ((base-dur number) (subdivisions list) (times integer))
    :initvals '(3000 '(1 2 3 4 5 6 7) 4)
    :indoc '("number" "list" "integer")
    :icon 000
    :doc "Outputs a rhythmicon as a MULTI-SEQ, given a fundamental duration (ms), a list of subdivisions, and a number of repetitions. For each rhythmic partial, the corresponding pitch is automatically assigned, using 3300 as the fundamental.
    "
    (setq onsets nil)
    (setq pitches nil)
    (setq durations nil)
    (setq vels (om-scale subdivisions 120 30))
    (setq f0 (mc->f 3300))
    (setq output nil)
    (setq velocities nil)
    (setq max-onset (* base-dur times))
    (loop for m in subdivisions and v in vels do
        (setq numnotes (nth-value 0 (ceiling (* times m))))
        (setq durs (repeat-n (/ base-dur m) numnotes))
        (setq pre-onsets (om-clip (dx->x 0 durs) 0 max-onset))
        (setq onsets (append onsets (list pre-onsets)))
        (setq durations (append durations (list (x->dx pre-onsets))))
        (setq pitches (append pitches (list (repeat-n (f->mc (* f0 m)) (* times (floor m))))))
        (setq velocities (append velocities (list (repeat-n v numnotes)))))
    (make-instance 'multi-seq :chord-seqs (reverse (loop for o in onsets and p in pitches and d in durations and v in velocities collect 
        (make-instance 'chord-seq :lmidic p :lonset o :ldur d :lvel v)))))

;--------------- Extract-channel ---------------
(defmethod! Extract-channel ((self chord-seq) (midi-chan integer))
    :initvals '((make-instance 'chord-seq :lmidic '((6000 7200) (6400 7000) 6700) :lchan '((1 1) (2 1) 3)) 1)
    :indoc '("chord-seq or multi-seq" "integer")
    :icon 000
    :doc "
        Extracts the notes with the specified midi-channel from a CHORD-SEQ or MULTI-SEQ 
    "
    (setq channels (lchan self))
    (setq cents nil) (setq onsets nil) (setq durations nil) (setq velocities nil) (setq offsets nil)
    (loop for chans in channels and i from 0 to (- (length channels) 1) do
        (setq current_chord nil)
        (setq current_dur nil)
        (setq current_vel nil)
        (setq current_offset nil)
        (loop for ch in chans and j from 0 to (- (length chans) 1) do
            (if (eq ch midi-chan)
                (progn
                    (setq current_chord (append current_chord (list (nth j (nth i (lmidic self))))))
                    (setq current_dur (append current_dur (list (nth j (nth i (ldur self))))))
                    (setq current_vel (append current_vel (list (nth j (nth i (lvel self))))))
                    (setq current_offset (append current_offset (list (nth j (nth i (loffset self)))))))))
        (if (not (eq current_chord nil))
            (progn 
                (setq cents (append cents (list current_chord)))
                (setq onsets (append onsets (list (nth i (lonset self)))))
                (setq durations (append durations (list current_dur)))
                (setq velocities (append velocities (list current_vel)))
                (setq offsets (append offsets (list current_offset))))))
    (make-instance 'chord-seq :lmidic cents :lonset onsets :ldur durations :lvel velocities :loffset offsets :lchan midi-chan))

(defmethod! Extract-channel ((self multi-seq) (midi-chan integer))
    (setq chord-seq-list (inside self))
    (make-instance 'multi-seq :chord-seqs (mapcar #'(lambda (input) (Extract-channel input midi-chan)) chord-seq-list)))

;--------------- Segment-seq ---------------
(defmethod! Segment-seq ((self chord-seq) (time-pt-list list) (samp-dur number) &optional (detection-mode '0) (clip-mode '0))
    :initvals '((make-instance 'chord-seq :lmidic '((6000 7200) (6400 7000) 6700) :lchan '((1 1) (2 1) 3)) (0 1000) 250 '0 '0)
    :indoc '("chord-seq or multi-seq" "list or number" "number" "menu" "menu")
    :icon 000
    :menuins '((3 (("detect onsets" '0) ("detect onsets and durations" '1))) (4 (("no clipping" '0) ("clip onsets" '1) ("clip durations" '2) ("clip onsets and durations" '3))))
    :doc "Extracts a segment from a CHORD-SEQ, given a list of time points and duration for all segments.
    "
    (setq cents nil) (setq onsets nil) (setq durations nil) (setq velocities nil) (setq offsets nil) (setq chans nil)
    (setq seq-onsets (lonset self))
    (loop for time-pt in time-pt-list do
        (setq end-time-pt (+ time-pt samp-dur))

        (loop for s-onset in seq-onsets and s-dur in (ldur self) and i from 0 to (- (length seq-onsets) 1) do

            (setq end-pts (om+ s-onset s-dur))
            (setq st-pts (om- end-pts s-dur))
            (setq time-points (mat-trans (list st-pts end-pts)))

            (setq current_chord nil)
            (setq current_dur nil)
            (setq current_vel nil)
            (setq current_offset nil)
            (setq current_chan nil)

            (loop for tp in time-points and j from 0 to (- (length time-points) 1) do
                (if 
                    (if (eq detection-mode '1)
                        (or
                            (and ; detect onset within range
                                (>= (first tp) time-pt) ; onset after st
                                (< (first tp) end-time-pt)) ; onset before end
                            (and ; detect ending within range
                                (> (second tp) time-pt) ; outset after st
                                (<= (second tp) end-time-pt)) ; and outset before end
                            (and  ; detect note during range
                                (<= (first tp) time-pt) ; onset before st
                                (>= (second tp) end-time-pt) ; and onset after end
                            )
                        )
                        (and ; detect onset within range only
                            (>= (first tp) time-pt)
                            (< (first tp) end-time-pt)))
                    (progn 
                        (setq maxdur (- end-time-pt (first tp)))
                        (setq out-dur (nth j (nth i (ldur self))))
                        (if (or (eq clip-mode 2) (eq clip-mode 3))
                            (setq out-dur (min out-dur maxdur)))

                        (setq current_chord (append current_chord (list (nth j (nth i (lmidic self))))))
                        (setq current_dur (append current_dur (list out-dur)))
                        (setq current_vel (append current_vel (list (nth j (nth i (lvel self))))))
                        (setq current_chan (append current_chan (list (nth j (nth i (lchan self))))))
                        (setq current_offset (append current_offset (list (nth j (nth i (loffset self)))))))))
            (if (not (eq current_chord nil))
                (progn 
                    (setq cents (append cents (list current_chord)))
                    (setq out-onset (nth i (lonset self)))
                    (if (or (eq clip-mode '1) (eq clip-mode '3))
                        (setq out-onset (max out-onset time-pt)))
                    (setq onsets (append onsets (list out-onset)))
                    (setq durations (append durations (list current_dur)))
                    (setq velocities (append velocities (list current_vel)))
                    (setq chans (append chans (list current_chan)))
                    (setq offsets (append offsets (list current_offset)))))))
    (setq min-onset (list-min (flat onsets)))
    (make-instance 'chord-seq :lmidic cents :lonset (om- onsets min-onset) :ldur durations :lvel velocities :loffset offsets :lchan chans))

(defmethod! Segment-seq ((self chord-seq) (time-pt number) (samp-dur number) &optional (detection-mode '0) (clip-mode '0))
    (Segment-seq self (list time-pt) samp-dur detection-mode clip-mode))

(defmethod! Segment-seq ((self multi-seq) (time-pt number) (samp-dur number) &optional (detection-mode '0) (clip-mode '0))
    (make-instance 'multi-seq :chord-seqs 
        (loop for seq in (inside self) collect
            (Segment-seq seq (list time-pt) samp-dur detection-mode clip-mode))))

(defmethod! Segment-seq ((self multi-seq) (time-pt list) (samp-dur number) &optional (detection-mode '0) (clip-mode '0))
    (make-instance 'multi-seq :chord-seqs 
        (loop for seq in (inside self) collect
            (Segment-seq seq time-pt samp-dur detection-mode clip-mode))))

;--------------- List-frames ---------------
(defmethod! List-frames ((in-list list) (size integer) &optional (hop 1))
    :initvals '((0 1 2 3 4 5 6) 2 1)
    :indoc '("list" "integer" "integer")
    :icon 000
    :doc "Parses a list into smaller lists (frames), given a frame length and a hop size.

    Examples:

    (list-frames '(0 1 2 3 4 5 6) 2 1)  =>  ((0 1) (1 2) (2 3) (3 4) (4 5) (5 6))
    (list-frames '(0 1 2 3 4 5 6) 2 2)  =>  ((0 1) (2 3) (4 5))
    "
    (setq indices (arithm-ser 0 (- size 1) 1))
    (setq x 0)
    (setq out nil)
    (setq max-index (- (+ (length in-list) 1) size))
    (while (< x max-index)
        (setq out (append out (list (posn-match in-list (om+ indices x) ))))
        (setq x (+ x (max 1 hop))))
    out)

;--------------- Markov-build ---------------
(defmethod! Markov-build ((data list) (order integer))
    :initvals '((0 2 1 3 2 4 3 5 4 6 5 4 3 2 1 0) 1)
    :indoc '("list" "integer")
    :icon 000
    :doc "Takes a list and computes a Nth-order transition probability matrix of the elements (states) in the list. MARKOV-BUILD is meant to be used along with MARKOV-RUN.
    "
    (setq data (List-frames data order 1))
    (setq thin-data (reverse (remove-dup data 'equal 1)))
    (setq dims (length thin-data))
    (setq matrix (mat-trans (append (list thin-data) (repeat-n (repeat-n 0 dims) dims))))
    (loop for thin-item in thin-data and row from 0 to (- (length thin-data) 1) do
        (loop for item in data and pos from 0 to (- (length data) 2) do
            (if (equal item thin-item)
                (progn
                    (setq next-item (nth (+ pos 1) data))
                    (setq col (position next-item thin-data :test 'equal))
                    (setf (nth (+ col 1) (nth row matrix)) (+ 1 (nth (+ col 1) (nth row matrix))))))))
    matrix
)

(defmethod! Markov-run ((matrix list) (iterations integer) &optional (initial nil) (mode '1))
    :initvals '('(((0) 0 0 1 0 0) ((1) 1/2 0 0 1/2 0) ((2) 0 2/3 0 0 1/3) ((3) 0 0 1 0 0) ((4) 0 0 0 1 0)) 5 nil '1)
    :indoc '("list" "integer" "atom or list" "menu")
    :icon 000
    :menuins '((3 (("no reset" '0) ("allow reset" '1))))
    :doc "Takes a transition probability matrix and outputs a sequence of a given length. MARKOV-RUN is meant to be used along with MARKOV-BUILD.
    "
    (setq states (car (mat-trans matrix)))
    (setq output nil)

    ; initialize first state
    (if (equal initial nil)
        (setq current-state (nth-random states))
        (progn
            (setq positions (get-posn initial (car (mat-trans states))))
            (if (eq positions nil)
                (error "
                    
                    INITIAL STATE MUST BELONG TO INPUT MATRIX.
                ")
                (setq current-state (nth (nth-random positions) states)))))
    (loop for rep from 0 to (- iterations 1) do
        (setq output (append output (list (car current-state))))
        (setq row-index (position current-state states :test 'equal))
        (setq weights (cdr (nth row-index matrix)))
        (if (> (reduce-tree weights #'+) 0)
            (progn 
                (setq choice (car (pick-random states weights 1)))
                (setq current-state choice))
            (if (equal mode '0)
                (progn
                    (setq rep iterations)
                    (print "MARKOV-RUN: Sequence reached dead end before desired length")
                )
                (setq current-state (nth-random states)))))
    output)

; (defmethod! Covariance ((data list))
;     :initvals '(((10 3 2) (-30 1 2) (0 0 1) (45 0 3) (-50 3 2)))
;     :indoc '("list")
;     :icon 000
;     :doc "Covariance"
;     (setq mat-center (flat (list-moments (mat-trans data) '(0))))
;     (setq datasize (length data))

; )
(defmethod! Histogram ((data list) (bins integer))
    :initvals '((0 5 2 8 4.5 9 1 0.3 4 5 0.4 -0.3 5 13) 5)
    :indoc '("list" "integer")
    :icon 000
    :doc "Histogram"
    (setq min-val (list-min data))
    (setq max-val (list-max data))
    (setq bin-size (/ (- max-val min-val) bins))
    (stable-sort data #'<)
    (setq histo nil)
    (loop for b from 1 to bins do
        (setq lower (+ min-val (* bin-size (- b 1))))
        (setq upper (+ lower bin-size))
        (setq counter 0)
        (loop for d in data collect
            (if (and (>= d lower) (< d upper))
                (setq counter (+ counter 1))))
        (setq histo (append histo (list (om-make-point lower counter)))))
    (make-instance 'bpf :point-list histo))
    
; ;--------------- PCA ---------------
; (defmethod! PCA ((data list))
;     :initvals '(((10 3 2) (-30 1 2) (0 0 1) (45 0 3) (-50 3 2)))
;     :indoc '("list")
;     :icon 000
;     :doc "PCA"
;     (mat-trans (loop for dim-row in (mat-trans data) collect 
;         (om- dim-row (car (list-moments dim-row '(0))))
;     ))

; )
#| 
    TODO:
        - KDTree
        - Granulate
        - PCA
 |#