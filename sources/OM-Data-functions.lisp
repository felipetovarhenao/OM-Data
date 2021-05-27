#|--------------------------------- * 
|     OM-Data library functions     |
|   [www.felipe-tovar-henao.com]    |
|               2021                |
* --------------------------------- *
|#

(in-package :om)

; -------------- (UNLISTED) FUNCTIONS -------------------
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

; --------------- string-rewrite ---------------
(defun string-rewrite (axiom rules iterations)
    (setq axiom (flat (list axiom) 1))
    (loop for n from 0 to (- iterations 1) do
        (loop for a in axiom and i from 0 to (- (length axiom) 1) do
            (loop for r in rules do
                (if (equal a (first r))
                    (progn
                        (setf (nth i axiom) (second r))))))
        (setq axiom (flat axiom 1)))
    axiom)

; -------------- M E T H O D S ---------------------

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
    :initvals '((0 0) ((1 1) (-1 1) (0 -1) (1 0) (2 1) (-1 -2)) nil)
    :indoc '("list" "list of lists" "list (optional)")
    :icon 000
    :numouts 2
    :doc "Sorts the lists based on the exhaustive nearest neighbor seach algorithm, using Euclidean distance as the sorting measurement.
    
    NOTE: All lists must have the same length.
    Example:  
    (NNS '(0 1 2 3) '((1 2 3 4) (0 0 2 3) (2 3 4 5)) nil) => ((0 0 2 3) (1 2 3 4) (2 3 4 5))
    "
    (setq nns-list nil) 
    (setq positions nil)
    (setq distances 
        (loop for l in other-lists and p from 0 to (- (length other-lists) 1) collect 
            (list l p (Euclidean-distance main-list l weights))))
    (stable-sort distances #'< :key #'third)
    (setq distances (mat-trans distances))
    (values-list (list (first distances) (second distances))))

#| (defmethod! NNS ((a number) (list-b list) (weights list))
    (setq out (NNS (list a) (mat-trans (list list-b)) weights))
    (values-list (list (flat (nth-value 0 out)) (nth-value 1 out)))) |#

; --------------- Optimal-sorting ---------------
(defmethod! Optimal-sorting ((st-list list) (other-lists list) (weights list))
    :initvals '((0 0 0 2) ((0 1 2 3) (2 3 4 5) (1 2 3 4)) nil)
    :indoc '("list (initial)" "list of lists" "list (optional)")
    :icon 000
    :numouts 2
    :doc "Sorts a list of lists such that the distance between adjacent lists is optimally minimized, given a starting list.
    
    NOTE: All lists must have the same length.
    
    Example:
    (optimal-sorting '(0 0 0 2) '((0 1 2 3) (2 3 4 5) (1 2 3 4)) nil) => ((0 0 0 2) (0 1 2 3) (1 2 3 4) (2 3 4 5))
    "  
    (setq neighbors (nth-value 0 (NNS st-list other-lists weights)))
    (setq remaining (copy-tree neighbors))
    (setq output nil)
    (loop for n in neighbors do
        (setq output (append output (list (car remaining))))
        (setq remaining (nth-value 0 (NNS (car remaining) (cdr remaining) weights))))
    (setq positions (nested-position output other-lists))
    (setq output (append (list st-list) output))
    (values-list (list output positions)))

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
    (setq sorted-distances (copy-tree distances))
    (stable-sort sorted-distances #'<)
    (+ (* accuracy (nth (nth 0 (get-posn (car sorted-distances) distances)) b-list)) (* a (- 1 accuracy))))

(defmethod! List-quantize ((a-list list) (b-list list) (accuracy number))
    (mapcar #'(lambda (input) 
        (List-quantize input b-list accuracy)) a-list))

(defmethod! List-quantize ((a-list chord-seq) (b-list list) (accuracy number))
    (setq mc (lmidic a-list))
    (setq mc (mapcar #'(lambda (input) 
        (List-quantize input b-list accuracy)) mc))
    (make-instance 'chord-seq 
        :lmidic mc
        :lonset (lonset a-list)
        :ldur (ldur a-list)
        :lvel (lvel a-list)
        :loffset (loffset a-list)
        :lchan (lchan a-list)))

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
    (setq chord-list (remove-dup (List-mod chord-list 1200) 'equal 1))
    (setq base-chord (copy-tree chord-list))
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
    (setq filled-range (Fill-range chord-list 0 12700))
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
    :outdoc '("mean st-dev skewness kurtosis")
    :doc "Computes the statistical moments of a list of values: population mean, population standard deviation, skewness and kurtosis. 
    The right inlet specifies which moments to put out, by index. For instance, (0 2) outputs the mean and skewness only. The list (0 1 2 3) outputs all four moments.
    
    Example:
    (list-moments '(5 2 3 4 5 6) '(0 1 2 3)) => (4.1666665 1.3437096 -0.30531597 1.8482844)
    " 
    (if (eq (depth data) 1)
        (if (> (length data) 1)
            (progn 
                (setq skewness nil) (setq kurtosis nil)
                (setq mean (/ (reduce #'+ data) (* 1.0 (length data))))
                (setq st-dev (sqrt 
                    (/ (reduce #'+ 
                        (loop for x in data collect
                            (expt (- x mean) 2))) 
                        (length data))))
                (if (> (abs st-dev) 0)
                    (progn 
                        (setq skewness
                            (/ (reduce #'+
                                (loop for x in data collect
                                    (expt (- x mean) 3))) 
                                (* (expt st-dev 3) (length data) 1)))
                        (if (> (abs skewness) 0)
                            (setq kurtosis
                                (/ (reduce #'+
                                    (loop for x in data collect
                                        (expt (- x mean) 4))) 
                                    (* (expt st-dev 4) (length data))))))
                    (progn (setq skewness nil) (setq kurtosis nil)))
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

#| ;--------------- K-means ---------------
(defmethod! K-means ((data list) (k integer) (weights list))
    :initvals '(((0 1 0) (-3 -1 2) (4 0 9) (-3 -5 -1) (0 4 -3) (2 1 -4)) 2 nil)
    :indoc '("list" "k (integer)" "weights (optional)")
    :icon 000
    :numouts 2
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
        (setq pk-centroids (copy-tree k-centroids))
        ; assign a k-centroid to each data item based on proximity
        (loop for ld in labeled-data and ld-pos from 0 to (length labeled-data) do
            (setq nearest-k (car (nth-value 0 (NNS (car ld) k-centroids weights))))
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
                    (setq new-centroid (flat (List-moments current-k (list 0))))
                    (setf (nth ck k-centroids) new-centroid))))

        ; stop loop if centroids do not change
        (setq convergence-flag (equal k-centroids pk-centroids)))

    ; group data by classes
    (stable-sort labeled-data #'< :key #'second)
    (setq output (loop for n from 0 to (- k 1) collect nil))
    (loop for ld in labeled-data do
        (setf (nth (second ld) output) (append (nth (second ld) output) (list (first ld)))))
    (values-list (list output (nested-position output data)))) |#

;--------------- K-means ---------------
(defmethod! K-means ((data list) (k integer) (weights list))
    :initvals '(((0 1 0) (-3 -1 2) (4 0 9) (-3 -5 -1) (0 4 -3) (2 1 -4)) 2 nil)
    :indoc '("list" "k (integer)" "weights (optional)")
    :icon 000
    :numouts 2
    :doc "Unsupervised data clustering algorithm.
    
    NOTE: All data items must have the same size. Weights are optional.
    
    Example:
    (k-means '((0 1 0) (-3 -1 2) (4 0 9) (-3 -5 -1) (0 4 -3) (2 1 -4)) 2 nil) => (((0 1 0) (-3 -1 2) (-3 -5 -1) (0 4 -3) (2 1 -4)) ((4 0 9)))
    " 
    
    ; clip k to suitable range
    (setq k (clip k 1 (- (length data) 1))) 

    ; copy data and add label slots
    (setq labeled-data 
        (loop for d in data and n from 0 to (- (length data) 1) collect (list d nil n)))

    ; initialize k-centroids and sort by first element
    (setq k-centroids (k-smart data k))
    (stable-sort k-centroids #'< :key #'first)

    ; convergence flag for loop
    (setq convergence-flag nil)

    ; ----- K_MEANS routine ----
    (while (eq convergence-flag nil)

        ; keep a history of last k-centroids
        (setq pk-centroids (copy-tree k-centroids))
        ; assign a k-centroid to each data item based on proximity
        (loop for ld in labeled-data and ld-pos from 0 to (length labeled-data) do
            (setq nearest-k (car (nth-value 0 (NNS (car ld) k-centroids weights))))
            (setq ck (position nearest-k k-centroids :test 'equal))
            (setf (nth ld-pos labeled-data) (list (car ld) ck (third ld))))
        
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
                    (setq new-centroid (flat (List-moments current-k (list 0))))
                    (setf (nth ck k-centroids) new-centroid))))

        ; stop loop if centroids do not change
        (setq convergence-flag (equal k-centroids pk-centroids)))

    ; group data by classes
    (stable-sort labeled-data #'< :key #'second)
    (setq output (repeat-n nil k))
    (setq positions (repeat-n nil k))
    (loop for ld in labeled-data do
        (setf (nth (second ld) output) (append (nth (second ld) output) (list (first ld))))
        (setf (nth (second ld) positions) (append (nth (second ld) positions) (list (third ld)))))
    (values-list (list output positions)))

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
    (setq a-table (remove-dup (flat (copy-tree a-list)) 'eq 1))
    (stable-sort a-table #'<)
    (setq ab-matrix nil)
    (loop for a in a-table do 
        (setq b-target (closest a b-list))
        (setq ab-matrix (append ab-matrix (list (list a (car b-target))))))
    (setq ab-positions (Nested-position a-list (first (mat-trans ab-matrix))))
    (setq nested-target (Nested-nth ab-positions (second (mat-trans ab-matrix))))
    (loop for a in a-list and nt in nested-target and f in interp-pts collect
        (nested-mix a nt f)))  

(defmethod! X-interpolation ((a-list chord-seq) (b-list list) (traj list))
    (setq mc (lmidic a-list))
    (setq mc (mapcar #'(lambda (input) 
        (X-interpolation input b-list traj)) mc))
    (make-instance 'chord-seq 
        :lmidic mc
        :lonset (lonset a-list)
        :ldur (ldur a-list)
        :lvel (lvel a-list)
        :loffset (loffset a-list)
        :lchan (lchan a-list)))

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
    (setq current (copy-tree a-list))
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
        (setq current (copy-tree scrambled)))
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
(defmethod! Rhythmicon ((base-dur number) (subdivisions list) (times integer) mode)
    :initvals '(3000 '(1 2 3 4 5 6 7) 4 'divisive)
    :indoc '("number" "list" "integer" "menu")
    :icon 000
    :menuins '((3 (("div" 'div) ("mult" 'mult))))
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
    (setq period (list-lcm subdivisions))
    (loop for m in subdivisions and v in vels do
        (cond
            (
                (equal mode 'div)
                (progn
                    (setq numnotes (nth-value 0 (ceiling (* m times))))
                    (setq durs (repeat-n (/ base-dur m) numnotes))
                    (setq pre-onsets (om-clip (dx->x 0 durs) 0 max-onset))
                    (setq onsets (append onsets (list pre-onsets)))
                    (setq durations (append durations (list (x->dx pre-onsets))))
                    (setq pitches (append pitches (list (repeat-n (f->mc (* f0 m)) numnotes))))
                    (setq velocities (append velocities (list (repeat-n v numnotes))))))
            (
                (equal mode 'mult)
                (progn
                    (setq numnotes (* (/ period m) times))
                    (setq durs (repeat-n (* base-dur m) numnotes))
                    (setq pre-onsets (dx->x 0 durs))
                    (setq onsets (append onsets (list pre-onsets)))
                    (setq durations (append durations (list (x->dx pre-onsets))))
                    (setq pitches (append pitches (list (repeat-n (f->mc (* f0 m)) numnotes))))
                    (setq velocities (append velocities (list (repeat-n v numnotes))))))))
    (make-instance 'multi-seq :chord-seqs (reverse (loop for o in onsets and p in pitches and d in durations and v in velocities collect 
        (make-instance 'chord-seq 
            :lmidic p 
            :lonset o 
            :ldur d 
            :lvel v)))))

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
    matrix)

;--------------- Markov-run ---------------
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
;--------------- Histogram ---------------
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
    (make-instance 'bpf :point-list histo :decimal 2))

;--------------- Multi-join---------------
(defmethod! Multi-join ((seqs list) &optional (mode '0) (concat-offset nil))
    :initvals '(nil '0 nil)
    :indoc '("multi-seq or poly" "join mode" "list")
    :icon 000
    :doc "Joins a list of score objects, either by merging or concatenating them." 
    :menuins '((1 (("concat" '0) ("merge" '1))))
    (if (and (equal mode '0) (equal concat-offset nil))
        (setq concat-offset 
            (cdr (dx->x 0 (loop for s in seqs collect
                (list-max (lonset s)))))))
    (setq out (car seqs))
    (setq seq-list (cdr seqs))
    (loop for s in seq-list and i from 0 to (- (length seq-list) 1) do
        (if (equal mode '0)
            (setq out (concat out s (nth i concat-offset)))
            (setq out (merger out s))))
    out)

(defmethod! Multi-join ((self multi-seq) &optional (mode '0) (concat-offset nil))
    (setq seqs (chord-seqs self))
    (Multi-join seqs mode concat-offset))

(defmethod! Multi-join ((self poly) &optional (mode '0) (concat-offset nil))
    (setq seqs (voices self))
    (Multi-join seqs mode concat-offset))

;--------------- Get-transients ---------------
(defmethod! Get-transients ((self chord-seq) (threshold number))
    :initvals '(nil 0.05)
    :indoc '("sequence" "number")
    :icon 000
    :doc "Outputs a list of onsets corresponding to detected transients in a score object." 
    (setq energy (list 0))
    (setq velocities (lvel self))
    (setq min-vel (list-min velocities))
    (setq velocities (append (list (list min-vel) (list min-vel)) velocities))
    (setq onsets (lonset self))
    (setq out nil)
    (loop for v in velocities and i from 0 to (- (length velocities) 1) do
        (setq vel (log (reduce #'+ v)))
        (if (> i 0)
            (progn 
                (setq val (max 0 (- vel pvel)))
                (if (> val threshold)
                    (setq out (append out (list (nth (- i 2) onsets)))))))
        (setq pvel vel))
    out)

(defmethod! Get-transients ((self multi-seq) (threshold number))
    (Get-transients (Multi-join self 1) threshold))

; --------------- L-system ---------------
(defmethod! L-system ((axiom number) (rules list) (generations integer))
    :initvals '('(x f) '((x (x + y f + + y f - f x - - f x f x - y f +)) (y (- f x + y f y f + + y f + f x - - f x - y))) 3)
    :indoc '("atom" "list" "integer")
    :icon 000
    :doc "Outputs a deterministically generated sequence of elements, given an axiom, a list of production rules, and a number of generations." 
    (string-rewrite axiom rules generations))

(defmethod! L-system ((axiom string) (rules list) (generations integer))
    (string-rewrite axiom rules generations))

(defmethod! L-system ((axiom list) (rules list) (generations integer))
    (string-rewrite axiom rules generations))
#| 
; --------------- 1D-Turtle ---------------
(defmethod! 1D-Turtle ((lsys list) (mag-rules list) (sign-rules list)(memory-rules list) &optional (x 0))
    :initvals '(([ x [ + f ] f - ] x [ + f ] f [ + [ f - ] f [ + x [ + f ] f ] x [ + f ] f ] [ f - ] f [ + x [ + f ] f ] x [ + f ] f) '((f 1)) '((+ 1) (- -1)) '(([ 1) (] 0)) 0)
    :indoc '("list" "list" "list" "list" "number")
    :icon 000
    :doc "Turtle"
    (setq y 0)
    (setq size 0)
    (setq out (list y))
    (setq memory nil)
    (loop for s in lsys do
        (loop for mr in mag-rules do
            (if (equal s (first mr))
                (progn
                    (setq y (+ y (* (second mr) size)))
                    (setq out (append out (list y))))))
        (loop for sr in sign-rules do
            (if (equal s (first sr))
                (progn
                    (setq size (+ size (second sr))))))
        (loop for mr in memory-rules do
            (if (equal s (first mr))
                (progn
                    (if (equal (second mr) 1)
                        (setq memory (append memory (list (list y size))))
                    )
                    (if (equal (second mr) 0)
                        (progn 
                            (setq state (car (last memory)))
                            (setq y (first state))
                            (setq size (second state))
                            (setq memory (butlast memory))))))))
    out) 
 |#

; --------------- 2D-Turtle ---------------
(defmethod! 2D-Turtle ((lsys list) (mag-rules list) (theta-rules list) (memory-rules list) &optional (theta 0))
    :initvals '(([ f - f ] + f [ f - f ] + f [ f - f ] + f [ f - f ] + f [ f - f ] + f [ f - f ] + f) '((f 1)) '((+ 60) (- -60)) '(([ 1) (] 0)) 0)
    :indoc '("list" "list" "list" "list" "number")
    :icon 000
    :doc "2D Turtle graphics"
    (setq x 0)
    (setq y 0)
    (setq mag 0)
    (setq memory nil)
    (setq out (list (om-make-point x y)))
    (loop for s in lsys do
        (loop for mr in mag-rules do
            (if (equal s (first mr))
                (progn
                    (setq mag (second mr))
                    (setq x (+ x (* mag (cos (deg->rad theta)))))
                    (setq y (+ y (* mag (sin (deg->rad theta)))))
                    (setq out (append out (list (om-make-point x y)))))))
        (loop for tr in theta-rules do
            (if (equal s (first tr))
                (progn
                    (setq theta (+ theta (second tr))))))
        (loop for mr in memory-rules do
            (if (equal s (first mr))
                (progn 
                    (if (equal (second mr) 1)
                        (setq memory (append memory (list (list x y theta mag))))
                    )
                    (if (equal (second mr) 0)
                        (progn
                            (setq state (car (last memory)))
                            (setq x (first state))
                            (setq y (second state))
                            (setq theta (third state))
                            (setq mag (fourth state))
                            (setq memory (butlast memory))))))))
    (make-instance 'bpc :point-list out))

; --------------- Risset rhythm---------------
(defun riss-onsets (onsets totdur period v mode)
    (setq totdur (list-max onsets))
    (setq v-exp (expt 2 v))
    (setq tl nil)
    (loop for x from 0 to (- v-exp 1) do
        (setq tl (append tl (om+ (butlast onsets) (om* totdur x)))))
    (setq tl (om/ (append tl (list (* v-exp totdur))) (* v-exp totdur)))
    (cond 
        (
            (equal mode 0)
            (setq te (om-log (om+ 1 tl) 2)))
        (
            (equal mode 1)
            (setq te (om- (om^ 2 tl) 1))))
    (om* te period))

(defun riss-rates (te totdur period v mode)
    (setq normalized (om/ te period))
    (if (equal mode 0)
        (setq normalized (om^ 2 (om+ v normalized)))
        (progn 
            (setq minval (expt 2 v))
            (setq maxval (expt 2 (+ v 1)))
            (setq normalized (om-scale (om* (om+ (om-log (om+ normalized 1) 2) 1) minval) maxval minval))))
    (setq rates (om* (/ (* totdur (log 2)) period) normalized))
    rates)

(defun riss-amps (te period voices v mode)
    (setq index-list (om/ te period))
    (if (equal mode 1)
        (setq index-list (om- 1 index-list)))
    (setq amp (om+ 0.5 (om* -0.5 (om-cos (om* (* pi 2) (om+ (/ v voices) (om/ index-list voices))))))))

;--------------- Risset-rhyhtm ---------------
(defmethod! Risset-rhythm ((self chord-seq) (speed number) (voices integer) (rep integer) &optional (onset-mode '0) (mc-mode '2))
    :initvals '((make-instance 'chord-seq :lmidic '(6000 6200 6400 6700) :lonset '(0 250) :ldur '(250)) 2.5 4 4 '0 '2)
    :indoc '("chord-seq" "number" "integer" "integer" "menu")
    :icon 000
    :menuins '((4 (("accelerando" '0) ("ritardando" '1))) (5 (("preserve pitch" '0) ("period-wise transp." '1) ("chord-wise transp." '2))))
    :doc "Risset rhythm."

    (setq mc-list (lmidic self))
    (setq onsets (lonset self))
    (setq durs (ldur self))
    (setq vels (lvel self))
    (setq offsets (loffset self))
    (setq chans (lchan self))
    (setq leg (legato self))

    (setq out nil)
    (setq totdur (list-max onsets))
    (setq period (* totdur speed))
    (setq fq-list (mc->f mc-list))

    (loop for v from 0 to (- voices 1) do
        (setq v-exp (expt 2 v))
        (setq new-onsets (riss-onsets onsets totdur period v onset-mode))
        (setq rates (riss-rates new-onsets totdur period v onset-mode))
        (setq amps (riss-amps new-onsets period voices v onset-mode))
        (setq new-mc nil)
        (setq new-vels nil)
        (setq new-durs nil)
        (loop for f in (flat (repeat-n fq-list v-exp) 1) and v in (flat (repeat-n vels v-exp) 1) and d in (flat (repeat-n durs v-exp) 1) and a in (butlast amps) and r in (butlast rates) and n from 0 to (- (length rates) 1) do
            (setq new-durs (append new-durs (list (om/ d r))))
            (if (equal mc-mode 1)
                (progn
                    (setq r-index (* (nth-value 0 (floor (/ n (length durs)))) (length durs)))
                    (setq r (nth r-index rates))))
            (if (> mc-mode 0)
                (setq f (om* r f)))
            (setq new-mc (append new-mc (list (f->mc f))))
            (setq new-vels (append new-vels (list (om* v a)))))
        (setq te (flat 
            (loop for x from 0 to (- rep 1) collect
                (om+ (butlast new-onsets) (* x period))
            ) 1))
        (setq te (append te (list (* rep (list-max new-onsets)))))
        (setq seq (make-instance 'chord-seq 
            :lmidic (flat (repeat-n new-mc rep) 1) 
            :lonset te 
            :ldur (flat (repeat-n new-durs rep) 1)  
            :lvel (flat (repeat-n new-vels rep) 1)
            :loffset (flat (repeat-n offsets (* rep v-exp)) 1)
            :lchan (flat (repeat-n chans (* rep v-exp)) 1)
            :legato leg))
        (setq out (append out (list seq))))
    (reverse out))

;--------------- List-median ---------------
(defmethod! List-median ((l list))
    :initvals '((4 2 5 6 1 2 0 -4))
    :indoc '("list")
    :icon 000
    :doc "Computes the median of a list of numeric values"
    (if (eq (depth l) 1)
        (progn 
            (setq size (length l))
            (stable-sort l #'<)
            (if (oddp size)
                (progn
                    (setq index (nth-value 0 (floor (/ size 2))))
                    (nth index l))
                (progn
                    (setq index (/ size 2))
                    (* (+ (nth (- index 1) l) (nth index l)) 0.5))))
        (mapcar #'(lambda (input) (List-median input)) l)))

;--------------- List-mode ---------------
(defmethod! List-mode ((l list))
    :initvals '((4 2 5 6 1 2 0 -4))
    :indoc '("list")
    :icon 000
    :doc "Computes the mode of a list of numeric values"
    (if (eq 1 (depth l))
        (progn 
            (setq counts
                (loop for x in l collect
                    (n-occurances x l)))
            (if (eq (list-max counts) (list-min counts))
                nil
                (remove-dup (posn-match l (get-posn (list-max counts) counts)) 'equal 1)))
        (mapcar #'(lambda (input) (List-mode input)) l)))

;--------------- nil-matrix ---------------
(defun nil-matrix (dims)
    (setq out nil)
    (loop for d in dims do
        (setq out (repeat-n out d)))
    out)

;--------------- Deep-replace ---------------
(defmethod! Deep-replace ((data list) (path list) new)
    :initvals '((((0 1) (2 3)) ((4 5 6) (7 8 (9)))) (1 1 2) "foo")
    :indoc '("list" "list" "list or atom")
    :icon 000
    :doc "Replaces any element in a list with a new element, given a list of positions corresponding to each level in the input list."
    (setq data-copy (copy-tree data))
    (setq levels (list data-copy))
    (loop for p in path do
        (setq x (nth p data-copy))
        (setq levels (append levels (list x)))
        (setq data-copy x))
    (setq levels (reverse levels))
    (setq path (append (list 0) (reverse path)))
    (setf (nth 0 levels) new)
    (loop for p in path and l in levels and n from 0 to (length levels) do
        (if (and (listp l) (> n 0))
            (setf (nth p l) new)
            (setq l new))
        (setq new l))
    (car (last levels)))

;--------------- KDTree ---------------
(defmethod! KDTree ((data-set list))
    :initvals '(((0 5) (-2 -3) (4 0) (5 -3) (5 7) (0.3 6.5)))
    :indoc '("list")
    :icon 000
    :numouts 2
    :outdoc '("nodes" "list")
    :doc "Computes a k-dimensional tree"
    (setq transp-data (mat-trans data-set))
    (setq nodes (list-median transp-data))
    (setq tree (nil-matrix (repeat-n 2 (length nodes))))
    (loop for d in data-set do
        (setq branch nil)
        (loop for n in nodes and i from 0 to (- (length nodes) 1) do
            (setq val (nth i d))
            (if (< val n)
                (setq pos 0)
                (setq pos 1))
            (setq branch (append branch (list pos))))
        (setq update (append (deep-nth tree branch) (list d)))
        (setq tree (deep-replace tree branch update)))
    (values-list (list nodes tree))) 

;--------------- Deep-nth ---------------
(defmethod! Deep-nth ((data list) (path list))
    :initvals '((((0 1) (2 3)) ((4 5 6) (7 8 (9)))) (1 1 2))
    :indoc '("list" "list")
    :icon 000
    :doc "Gets the nth element of a nested list, given a list of positions corresponding to each level in the input list."
    (setq out nil)
    (if (eq 1 (depth path))
        (progn 
            (loop for p in path do
                (setq data (nth p data)))
        data)
        (mapcar #'(lambda (input) (deep-nth data input)) path)))

;--------------- KNN ---------------
(defmethod! KNN ((data list) (k integer) (tree-nodes list) (kd-tree list) &optional (weights nil))
    :initvals '((0.5 0.2) 1 (0 1) ((((-2 -2)) ((-1 1))) (((0 0)) ((2 2) (1 1)))) nil)
    :indoc '("list" "integer" "list" "list" "list")
    :icon 000
    :doc "Finds the K-nearest neighbors within a given KDTree."
    (if (equal (depth data) 1)
        (progn
            (setq branch nil)
            (loop for n in tree-nodes and i from 0 to (- (length tree-nodes) 1) do
                (setq val (nth i data))
                (if (< val n)
                    (setq pos 0)
                    (setq pos 1))
                (setq branch (append branch (list pos))))
            (setq neighbors (deep-nth kd-tree branch))
            (setq out (first-n (nth-value 0 (NNS data neighbors weights)) k))
            (if (eq k 1)
                (setq out (car out)))
            out)
        (mapcar #'(lambda (input) (KNN input k tree-nodes kd-tree)) data)))

;--------------- List-covariance ---------------
(defmethod! List-covariance ((data list))
    :initvals '(((-2 2) (0 0) (1 -1) (-1 1) (-4 3) (-3 4)))
    :indoc '("list")
    :icon 000
    :doc "Computes the covariance of a given set of data points"
    (setq means (flat (list-moments (mat-trans data) '(0))))
    (setq out
        (reduce #'+
            (loop for d in data collect
                (reduce #'* 
                    (loop for x in d and m in means collect
                        (- x m))))))
    (/ out (- (length data) 1)))

;--------------- List-correlation ---------------
(defmethod! List-correlation ((data list))
    :initvals '(((-2 2) (0 0) (1 -1) (-1 1) (-4 3) (-3 4)))
    :indoc '("list")
    :icon 000
    :doc "Computes the correlation of a given set of data points"
    (setq trans-data (mat-trans data))
    (setq variances (loop for td in trans-data collect
        (sqrt (List-covariance (mat-trans (list td td))))))
    (/ (List-covariance data) (reduce #'* variances)))

;--------------- Plot-points ---------------
(defmethod! Plot-points ((points list))
    :initvals '(((-2 2) (0 0) (1 -1) (-1 1) (-4 3) (-3 4)))
    :indoc '("list")
    :icon 000
    :doc "Plots a list of data points in a bpc-lib or 3dc-lib, depending on the dimensionality of the data."
    (setq bp-list nil)
    (setq dims (length (mat-trans points)))
    (loop for p in points do
        (if (equal dims 2)
            (progn 
                (setq omp (repeat-n (om-make-point (first p) (second p)) 2))
                (setq bp-list (append bp-list (list (make-instance 'bpc :point-list omp :decimals 2))))))
        (if (equal dims 3)
            (progn
                (setq omp (repeat-n (make-3dpoint :x (first p) :y (second p) :z (third p)) 2))
                (setq bp-list (append bp-list (list (make-instance '3dc :point-list omp :decimals 2)))))))
    (if (equal dims 2)
        (setq out (make-instance 'bpc-lib :bpf-list bp-list)))
    (if (equal dims 3)
        (setq out (make-instance '3dc-lib :bpf-list bp-list)))
    out)

;--------------- Neo-Riemannian Transformations ---------------
(defun triad-posn (mc)
    (setq pc-set (nth-value 1 (om// mc 1200)))
    (setq pc-set (om- pc-set (list-min pc-set)))

    (setq maj0 (list 0 400 700))
    (setq maj1 (list 800 0 300))
    (setq maj2 (list 500 900 0))

    (setq min0 (list 0 300 700))
    (setq min1 (list 900 0 400))
    (setq min2 (list 500 800 0))

    (setq maj0-diff (x-diff pc-set maj0))
    (setq maj1-diff (x-diff pc-set maj1))
    (setq maj2-diff (x-diff pc-set maj2))

    (setq min0-diff (x-diff pc-set min0))
    (setq min1-diff (x-diff pc-set min1))
    (setq min2-diff (x-diff pc-set min2))

    (setq major-diff (list maj0-diff maj1-diff maj2-diff))
    (setq minor-diff (list min0-diff min1-diff min2-diff))

    (setq maj-triads (list maj0 maj1 maj2))
    (setq min-triads (list min0 min1 min2))

    (setq diff-list (list major-diff minor-diff))
    (setq triad-list (list maj-triads min-triads))

    (setq triad-type nil)
    (setq out nil)

    (loop for diffs in diff-list and tt from 0 to 1 and triads in triad-list do
        (loop for d in diffs and tri in triads do
            (if (equal d nil)
                (progn 
                    (setq triad-type tt)
                    (loop for pc in pc-set do
                        (loop for n in tri and x from 0 to 2 do
                            (if (equal pc n)
                                (setq out (append out (list x))))))))))
    (list triad-type out))

(defun p-nrt (tri-type posn)
    (setq out nil)
    (if (eq tri-type 0)
        (setq out (posn-match (list 0 -100 0) posn)))
    (if (eq tri-type 1)
        (setq out (posn-match (list 0 100 0) posn)))
    out)

(defun r-nrt (tri-type posn)
    (setq out nil)
    (if (eq tri-type 0)
        (setq out (posn-match (list 0 0 200) posn)))
    (if (eq tri-type 1)
        (setq out (posn-match (list -200 0 0) posn)))
    out)

(defun l-nrt (tri-type posn)
    (setq out nil)
    (if (eq tri-type 0)
        (setq out (posn-match (list -100 0 0) posn)))
    (if (eq tri-type 1)
        (setq out (posn-match (list 0 0 100) posn)))
    out)

(defun apply-nrt (mc nrt-list)
    (setq triad (copy-tree mc))
    (if (atom nrt-list)
        (setq nrt-list (list nrt-list)))
    (loop for tr in nrt-list do
        (setq label (triad-posn triad))
        (setq triad-type (first label))
        (setq posn (second label))
        (cond 
            (
                (equal tr 'p)
                (setq triad (om+ triad (p-nrt triad-type posn))))  
            (
                (equal tr 'r)
                (setq triad (om+ triad (r-nrt triad-type posn))))  
            (
                (equal tr 'l)
                (setq triad (om+ triad (l-nrt triad-type posn))))))
    triad)

(defmethod! NRT ((mc list) (transformations list))
    :initvals '((6000 5500 7600 7200) '(r l (l p) ))
    :indoc '("list" "list")
    :icon 000
    :doc "Performs Neo-riemannian transformations on a starting triadic chord. The possible transformations are p l and r. Compound transformations are specified as lists of these 3 basic transformations."
    (if (eq (depth mc) 1)
        (loop for tr in transformations collect
            (setq mc (apply-nrt mc tr)))
        (mapcar #'(lambda (input) (apply-nrt input transformations)) mc)))

;--------------- Make-sieve ---------------
(defmethod! Make-sieve ((list list) (reps integer) sieve-mode sieve-type &optional (offset '0))
    :initvals '((2 3) 1 'union 'nil 0)
    :indoc '("list" "integer" "menu" "menu" "number")
    :menuins '(
        (2 (("union" 'union) ("diff" 'diff)))
        (3 (("nil" 'nil) ("complement" 'complement))))
    :icon 000
    :doc "Builds N full periods of a sieve, based on a list of integers. Make-sieve is meant to be a compact version of OM's native CRIBLE class and functions"
    (setq list (remove 0 list))
    (setq period (+ offset (* reps (list-lcm list))))
    (setq sieves (loop for l in list collect
        (arithm-ser offset period l)))
    (cond
        (
            (equal sieve-mode 'union)
            (setq out (list-union sieves)))
        (
            (equal sieve-mode 'diff)
            (setq out (list-diff sieves))))
    (if (equal sieve-type 'complement)
        (setq out (list-diff (list (arithm-ser offset period 1) (flat (list out))))))
    (stable-sort out #'<))

(defun list-union (list)
    (setq out (car list))
    (loop for l in (cdr list) do
        (setq out (x-union out l)))
    out)

(defun list-intersect (list)
    (setq out (car list))
    (loop for l in (cdr list) do
        (setq out (x-intersect out l)))
    out)

(defun list-diff (list)
    (setq out (car list))
    (loop for l in (cdr list) do
        (setq out (x-diff out l)))
    out)

(defun list-lcm (list)
    (setq out (car list))
    (loop for l in (cdr list) do
        (setq out (lcm out l)))
    out)

(defun list-gcd (list)
    (setq out (car list))
    (loop for l in (cdr list) do
        (setq out (gcd out l)))
    out)

;--------------- Vieru-sequence ---------------
(defmethod! Vieru-seq ((seq list) (n-tiers integer))
    :initvals '((4 1 2 1 4) 1)
    :indoc '("list" "integer" "integer")
    :icon 000
    :doc "Takes the ascending modular differences between adjacent values. Based on Anatol Vieru's modal sequences.

    Examples:
    (vieru-seq '(4 1 2 1 4) 1) => ((9 1 11 3 0))
    (vieru-seq '(2 1 4 1 2) 3) => ((9 3 7 1 0) (4 4 4 9 9) (0 0 5 0 5))
    "
    (setq mod-n (reduce #'+ seq))
    (setq seq (nth-value 1 (om// (append seq (list (car seq))) mod-n)))
    (setq out nil)
    (loop for x from 1 to n-tiers do
        (setq diff-seq nil)
        (loop for i from 0 to (- (length seq) 2) do
            (setq current (nth i seq))
            (setq next (nth (+ i 1) seq))
            (if (>= next current)
                (setq val (abs (- next current)))
                (setq val (abs (- (+ next mod-n) current))))
            (setq diff-seq (append diff-seq (list val))))
        (setq out (append out (list diff-seq)))
        (setq seq (append diff-seq (list (car diff-seq)))))
    out)

;--------------- List-symmetries ---------------
(defmethod! List-symmetries ((in-list list) mode &optional (tolerance 0.0))
    :initvals '((3 1 1 3 2 1 1 2) 'rotations 0.0)
    :indoc '("list" "menu" "number")
    :outdoc '("found symmetries" "assymmetry factor")
    :menuins '((1 (("permutations" 'permutations) ("rotations" 'rotations))))
    :icon 000
    :numouts 2
    :doc "Searches for symmetric or near-symmetric permutations of a list. The search can be done on all permutations or rotations only
    
    Examples:
    (List-symmetries '(3 1 1 3 2 1 1 2) 'rotations 0.0) => [ ((1 3 2 1 1 2 3 1) (1 2 3 1 1 3 2 1)) (0 0) ]
    (List-symmetries '(3 2 3 2 1) 'permutations 0.0) => [ ((3 2 1 2 3) (2 3 1 3 2)) (0 0) ]"
    (if (> tolerance 0.0)
        (progn 
            (setq sorted-list (sort-list in-list))
            (setq maxdist (reduce-tree (om-abs (om- sorted-list (reverse sorted-list))) #'+))))
    (cond 
        (
            (equal mode 'rotations)
            (setq forms (loop for n from 0 to (- (length in-list) 1) collect (rotate in-list n))))
        (
            (equal mode 'permutations)
            (setq forms (remove-dup (permutations in-list) 'equal 1))))
    (setq out nil)
    (loop for f in forms do
        (setq rf (reverse f))
        (if (> tolerance 0.0)
            (progn 
                (setq dist (/ (reduce-tree (om-abs (om- rf f))  #'+) maxdist))
                (if (<= dist tolerance)
                    (setq out (append out (list (list f dist))))))
            (if (equal f rf)
                (setq out (append out (list (list f 0)))))))
    (stable-sort out #'< :key 'second)
    (setq out (mat-trans out))
    (values-list (list (first out) (second out))))

;--------------- Array (Matrix) operations ---------------
(defun simple-arr-determinant (arr)
    (setq dims (array-dimensions arr))
    (if (equal dims (list 2 2))
        (- (* (aref arr 0 0) (aref arr 1 1)) (* (aref arr 0 1) (aref arr 1 0)))))

(defun list-dims (l)
    (setq out nil)
    (if (listp l)  
        (progn 
            (setq out (append out (list (length l))))
            (if (listp (car l))
                (setq out (append out (list (list-dims (car l))))))))
    (flat out))
(defun list->array (l)
    (make-array (list-dims l) :initial-contents l))

(defun col-row-exclude (arr row col)
    (setq dims (array-dimensions arr))
    (setq out nil)
    (loop for r from 0 to (- (first dims) 1) do
        (if (not (eq r row))
            (progn 
                (setq row nil)
                (loop for c from 0 to (- (second dims) 1) do    
                    (if (not (eq c col))
                        (progn 
                            (setq val (aref arr r c))
                            (setq row (append row (list val))))))
                (setq out (append out (list row))))))
    (list->array out))

(defun arr-determinant (arr)
    (setq dims (array-dimensions arr))
    (setq out nil)
    (if (equal dims (list 2 2))
        (setq out (simple-arr-determinant arr))
        (progn
            (setq first-row (flat (get-arr-row arr 0)))
            (setq signs (om- (om* 2 (nth-value 1 (om// (arithm-ser 1 (length first-row) 1) 2))) 1))
            (setq out (reduce #'+ (loop for x in first-row and n from 0 to (- (length first-row) 1) and s in signs collect
                (* s x (arr-determinant (col-row-exclude arr 0 n))))))))
    out)

(defun get-arr-rows (arr rows)
    (setq rows (flat (list rows)))
    (setq dims (array-dimensions arr))
    (setq out nil)
    (loop for r from 0 to (- (first dims) 1) do
        (if (not (eq nil (member r rows)))
            (progn 
                (setq row nil)
                (loop for c from 0 to (- (second dims) 1) do    
                    (setq val (aref arr r c))
                    (setq row (append row (list val))))
                (setq out (append out (list row))))))
    out)

(defun get-arr-cols (arr cols)
    (setq cols (flat (list cols)))
    (setq dims (array-dimensions arr))
    (setq out nil)
    (loop for r from 0 to (- (first dims) 1) do
        (setq row nil)
        (loop for c from 0 to (- (second dims) 1) do    
            (if (not (eq nil (member c cols)))
                (progn
                    (setq val (aref arr r c))
                    (setq row (append row (list val))))))
        (setq out (append out (list row))))
    out)

(defun covariance-matrix (data)
    (setq trans-data (mat-trans data))
    (loop for d1 in trans-data collect
        (loop for d2 in trans-data collect
            (list-covariance (mat-trans (list d1 d2))))))

(defun identity-matrix (dims)
    (loop for x from 0 to (- dims 1) collect
        (subs-posn (repeat-n 0 dims) x 1)))

(defun dot-product (v1 v2)
    (reduce #'+ (om* v1 v2)))

(defun matrix-mult (m1 m2)
    (setq m2 (mat-trans m2))
    (loop for row in m1 collect
        (loop for col in m2 collect
            (dot-product row col))))

#| (defun chord-autodetect (self threshold)
    (setq seq-dur (list-max (lonset self)))
    (setq seg-dur 1250)
    (setq overlap-factor 2)
    (setq max-dist 3.4641016)
    (setq threshold (* max-dist threshold))
    (setq markers (arithm-ser 0 seq-dur seg-dur (/ seg-dur overlap-factor)))
    (setq out nil)

    (loop for m in markers and n from 0 to (- (length markers) 1) do
        (setq seg (segment-seq self m seg-dur 1 0))
        (setq mc-list (lmidic seg))
        (setq chroma (flat (chroma-count mc-list)))
        (setq chroma (om/ chroma (list-max chroma)))
        (if (> n 0)
            (progn
                (setq dist (Euclidean-distance chroma pchroma nil))
                (if (> dist threshold)
                    (setq out (append out (list m)))
                )
            )
        )
        (setq pchroma (copy-tree chroma))
        (setq pmarker m)
    )
    out
) |#

#| 
    TODO:
        - Granulate
        - PCA
        - Chroma vector with durations and velocities
        - Rhythmic-distribution
 |#