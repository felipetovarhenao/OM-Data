#|--------------------------------- * 
|     OM-Data library functions     |
|   [www.felipe-tovar-henao.com]    |
|               2021                |
* --------------------------------- *
|#

(in-package :om)

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

;--------------- K-means ---------------

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

;--------------- DTW ---------------

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

;--------------- nil-matrix ---------------
(defun nil-matrix (dims)
    (setq out nil)
    (loop for d in dims do
        (setq out (repeat-n out d)))
    out)

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


