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

    Arguments:

    - <main-list>: reference values or list of values.
    - <other-lists>: lists of lists of values to calculate distance from <main-list>.
    - <weights>: list of weight values for each element in <main-list>.
    
    NOTE: All lists must have the same length.
    Example:  
    (NNS '(0 1 2 3) '((1 2 3 4) (0 0 2 3) (2 3 4 5)) nil) => ((0 0 2 3) (1 2 3 4) (2 3 4 5))
    "
    (let*
        (
            (nns-list nil) 
            (positions nil)
            (distances 
                (loop for l in other-lists and p from 0 to (- (length other-lists) 1) collect 
                    (list l p (Euclidean-distance main-list l weights)))))
        (stable-sort distances #'< :key #'third)
        (setf distances (mat-trans distances))
        (values-list (list (first distances) (second distances)))))

(defmethod! NNS ((a number) (list-b list) (weights list))
    (let* 
        (
            (out nil))
        (if (equal (depth list-b) 1)
            (setf list-b (mat-trans (list list-b))))
        (setf out (multiple-value-list (NNS (list a) list-b weights)))
        (values-list (list (flat (first out)) (second out)))))

;--------------- K-means ---------------

(defun k-smart (data k)
    (let*
        (
            (k-centroids (list (nth-random data)))
            (x 1))
        (setf data (remove (car k-centroids) data :test 'equal))
        (while (< x k)
            (let*
                (
                    (w-prob nil)
                    (distances nil)
                    (new-ck nil)
                )
                (loop for d in data do
                    (setf distances (loop for kc in k-centroids collect
                        (expt (Euclidean-distance d kc nil) 2)))
                    (setf w-prob (append w-prob (list (list-min distances))))
                )
                (setf new-ck (Pick-random data w-prob 1))
                (setf data (remove new-ck data :test 'equal))
                (setf k-centroids (append k-centroids new-ck))
                (setf x (+ x 1))))
        k-centroids))

(defmethod! K-means ((data list) (k integer) (weights list))
    :initvals '(((0 1 0) (-3 -1 2) (4 0 9) (-3 -5 -1) (0 4 -3) (2 1 -4)) 2 nil)
    :indoc '("list" "k (integer)" "weights (optional)")
    :icon 000
    :numouts 2
    :outdoc '("clustered data" "clustered positions")
    :doc "Unsupervised data clustering algorithm.
    
    NOTE: All data items must have the same size. Weights are optional.

    Arguments:

    - <data>: a list of lists.
    - <k>: number of assumed classes in <data>.
    - <weights>: list of weight values for each element in <data>.
    
    Example:
    (k-means '((0 1 0) (-3 -1 2) (4 0 9) (-3 -5 -1) (0 4 -3) (2 1 -4)) 2 nil) => (((0 1 0) (-3 -1 2) (-3 -5 -1) (0 4 -3) (2 1 -4)) ((4 0 9)))
    " 
    (let*
        (
            (output (repeat-n nil k))
            (positions (repeat-n nil k))
            (labeled-data nil)
            (k-centroids nil)
            (convergence-flag nil))
        ; clip k to suitable range
        (setf k (clip k 1 (- (length data) 1))) 
        ; copy data and add label slots
        (setf labeled-data 
            (loop for d in data and n from 0 to (- (length data) 1) collect (list d nil n)))
        ; initialize k-centroids and sort by first element
        (setf k-centroids (k-smart data k))
        (stable-sort k-centroids #'< :key #'first)

        ; ----- K_MEANS routine ----
        (while (eq convergence-flag nil)
            (let*
                (
                    ; keep a history of last k-centroids
                    (pk-centroids (copy-tree k-centroids))
                    (nearest-k nil)
                    (ck nil)
                    (t-labeled-data nil))
                ; assign a k-centroid to each data item based on proximity
                (loop for ld in labeled-data and ld-pos from 0 to (length labeled-data) do
                    (setf nearest-k (car (nth-value 0 (NNS (car ld) k-centroids weights))))
                    (setf ck (position nearest-k k-centroids :test 'equal))
                    (setf (nth ld-pos labeled-data) (list (car ld) ck (third ld))))
                
                ; update k-centroids 
                (setf t-labeled-data (mat-trans labeled-data))
                (loop for ck from 0 to (- k 1) do
                    (let*
                        (
                            (current-k nil)
                            (new-centroid nil)
                        )
                        (loop for ld in labeled-data do
                            (if (eq ck (second ld))
                                (setf current-k (append current-k (list (first ld))))))
                        (setf current-k (mat-trans current-k))
                        (if (not (equal current-k nil))
                            (progn 
                                (setf new-centroid (flat (List-moments current-k (list 0))))
                                (setf (nth ck k-centroids) new-centroid)))))
                ; stop loop if centroids do not change
                (setf convergence-flag (equal k-centroids pk-centroids))))
        ; group data by classes
        (stable-sort labeled-data #'< :key #'second)
        (loop for ld in labeled-data do
            (setf (nth (second ld) output) (append (nth (second ld) output) (list (first ld))))
            (setf (nth (second ld) positions) (append (nth (second ld) positions) (list (third ld)))))
        (values-list (list output positions))))

;--------------- DTW ---------------

;--------------- dtw-path-cost ---------------
(defun dtw-path-cost (a-list b-list)
    (let*
        (
            ; build cost matrix
            (matrix (make-array (list (length a-list) (length b-list)))))
        (loop for a in a-list and x from 0 to (length a-list) do
            (loop for b in b-list and y from 0 to (length b-list) do
                (let* 
                    (
                        (leftval nil)
                        (diagval nil)
                        (topval nil)
                        (options nil)
                        (cost nil))
                    (if (> x 0) 
                        (setf leftval (aref matrix (- x 1) y))
                        (setf leftval nil))
                    (if (> y 0) 
                        (setf topval (aref matrix x (- y 1)))
                        (setf topval nil)) 
                    (if (and (> y 0) (> x 0)) 
                        (setf diagval (aref matrix (- x 1) (- y 1)))
                        (setf  diagval nil))
                    (setf options (remove nil (list leftval topval diagval)))
                    (if (not (eq options nil))
                        (setf cost (list-min options))
                        (setf cost 0))
                    (setf (aref matrix x y) (+ (Euclidean-distance a b nil) cost)))))
        (let*
            (
                ; find optimal path
                (x (- (length a-list) 1)) 
                (y (- (length b-list) 1))
                (path-posn (list (list x y)))
                (path-cost (aref matrix x y)))
            (while (or (> x 0) (> y 0))
                (let* 
                    (
                        (leftval nil)
                        (diagval nil)
                        (topval nil) 
                        (options nil)
                        (minval nil))
                    (if (> x 0) 
                        (setf leftval (aref matrix (- x 1) y))
                        (setf leftval nil))
                    (if (> y 0) 
                        (setf topval (aref matrix x (- y 1)))
                        (setf topval nil)) 
                    (if (and (> y 0) (> x 0)) 
                        (setf diagval (aref matrix (- x 1) (- y 1)))
                        (setf diagval nil))
                    (setf options (list leftval topval diagval))
                    (setf minval (list-min options))
                    (cond
                        (
                            (eq minval leftval)
                            (setf x (- x 1)))
                        (
                            (eq minval topval)
                            (setf y (- y 1)))
                        (
                            (eq minval diagval)
                            (progn 
                                (setf x (- x 1))
                                (setf y (- y 1)))))
                    (setf path-cost (+ path-cost minval))
                    (setf path-posn (append path-posn (list (list x y))))))
            (list path-cost (reverse path-posn)))))

(defmethod! DTW ((list-a list) (list-b list))
    :initvals '((0 1 2 3 4 5) ((1 1 2 3 5) (0 0 1 2 3 3 4 5) (1 3 4 5)))
    :indoc '("list" "list of lists")
    :icon 000
    :numouts 2
    :outdoc '("sorted lists" "sorted")
    :doc "Sorts the lists in second input using Dynamic Time Warping.
    
    Arguments:

    - <list-a>: a list of numeric values.
    - <list-b>: a list of lists of numeric values.

    Example:
    (dtw '(0 1 2 3 4 5) '(1 1 2 3 5)) => [ 9 ((0 0) (1 0) (1 1) (2 2) (3 3) (4 4) (5 4)) ]
    "
    (let*
        (
            (out nil)
            (dtw-list))
        (cond
            (
                (eq (depth list-a) (depth list-b))
                (progn
                    (setf out (dtw-path-cost list-a list-b))
                    (values-list (list (first out) (second out)))))
            (
                (< (depth list-a) (depth list-b))
                (setf dtw-list nil)
                (loop for b in list-b do
                    (setf out (list (append (dtw-path-cost list-a b) (list b))))
                    (setf dtw-list (append dtw-list out)))
                (stable-sort dtw-list #'< :key #'first)
                (setf out (mat-trans dtw-list))
                (values-list (list (third out) (second out)))))))

;--------------- DTW-align ---------------
(defmethod! DTW-align ((a list) (b list) (posn list))
    :initvals '((0 1 2 3 4 2) (0 2 3 4 1) ((0 0) (1 0) (2 1) (3 2) (4 3) (5 4)))
    :indoc '("list" "list" "list")
    :icon 000
    :doc "Works in combination with DTW, aligning the values of the main time series to the others. Use the left input of DTW as in1, and the outputs of DTW as in2 and in3, respectively.
    
    Arguments:

    - <a>: times series A.
    - <b>: time series B.
    - <posn>: list of lists specifying the matching positions between <a> and <b>.

    Example: 
    (dtw-align '(0 1 2 3 4 2) '(0 2 3 4 1) '((0 0) (1 0) (2 1) (3 2) (4 3) (5 4))) => ((0 0) (1 0) (2 2) (3 3) (4 4) (2 1))
    "
    (if (eq (depth a) (depth b))
        (loop for p in posn collect
            (list (nth (first p) a) (nth (second p) b)))
        (mapcar #'(lambda (input1 input2) (DTW-align a input1 input2)) b posn)))

;--------------- nil-matrix ---------------
(defun nil-matrix (dims)
    (let*
        (
            (out nil))
        (loop for d in dims do
            (setf out (repeat-n out d)))
        out))

;--------------- KDTree ---------------
(defmethod! KDTree ((data-set list) &optional (max-k 4))
    :initvals '(((4 0) (4 3) (-4 1) (1 3) (1 -3) (0 3) (-2 -2) (-2 3)) 4)
    :indoc '("list" "integer")
    :icon 000
    :numouts 3
    :outdoc '("nodes" "tree" "positions")
    :doc "Computes a non-recursive k-dimensional tree to be used in tandem with KNN.
    
    Arguments:

    - <data-set>: list of lists.

    &optional
    - <max-k>: max number of dimensions. (numbranches = 2^dimensions)

    Output:
    - <nodes>
    - <tree>
    - <positions>
    "
    (let*
        (
            (nodes (first-n (list-median (mat-trans data-set)) max-k))
            (branches (repeat-n nil (expt 2 (length nodes))))
            (positions (copy-tree branches))
            (posn nil))
        (loop for x in data-set and i from 0 to (- (length data-set) 1) do
            (setf posn (branch-posn (first-n x max-k) nodes))
            (setf (nth posn branches) (append (nth posn branches) (list x)))
            (setf (nth posn positions) (append (nth posn positions) (list i)))
        )
        (values-list (list nodes branches positions))))

;--------------- KNN ---------------
(defmethod! KNN ((data list) (tree-nodes list) (kd-tree list) (tree-posns list) &optional (k 1) (weights nil))
    :initvals '((0 0) (1.0 -0.5) (((1 -2) (2 -4)) ((3 -2) (4 -3)) ((-4 4)) ((3 -1) (3 4) (3 4))) ((2 6) (0 4) (5 7) (1 3)) 1 nil)
    :indoc '("list" "list" "list" "list" "integer")
    :icon 000
    :numouts 2
    :outdoc '("nearest neighbors" "neighbors' positions")
    :doc "Finds the K-nearest neighbors within a given KDTree.

    Arguments:

    - <data>: list of lists
    - <tree-nodes>: nodes from output 1 of KDTREE.
    - <kd-tree>: tree from output 2 of KDTREE.
    - <tree-posns>: element positions from output 3 of KDTREE
    
    &optional:
    - <k>: number of desired closest neighbors.
    - <weights>: list of weight values for each element in <data>.
    "
    (let* 
        (
            (neighborhood nil)
            (neighbors nil)
            (flag t)
            (numbranches (expt 2 (length tree-nodes)))
            (branch-id nil)
            (positions nil)
            (main-branch nil)
            (search-count 1)
        )
        (if (equal (depth data) 1)
            (progn
                (setf branch-id (branch-posn data tree-nodes))
                (setf main-branch branch-id)
                (while flag
                    (setf neighborhood (nth branch-id kd-tree))
                    (setf positions (nth branch-id tree-posns))
                    (if (equal neighborhood nil)
                        (progn
                            (setf branch-id 
                                (list-wrap (wedge-sum main-branch search-count) 0 numbranches))
                            (setf search-count (+ search-count 1)))
                        (setf flag nil)))
                (let*
                    (
                        (nns-vals nil)
                        (nns-posn nil))
                    (setf neighbors (multiple-value-list (nns data neighborhood weights)))
                    (setf nns-vals (first-n (car neighbors) k))
                    (setf nns-posn (first-n (posn-match positions (second neighbors)) k))
                    (values-list (list nns-vals nns-posn))))
            (values-list (mat-trans (mapcar #'(lambda (in) (multiple-value-list (KNN in tree-nodes kd-tree tree-posns k weights))) data))))))

;--------------- Markov-build ---------------
(defmethod! Markov-build ((data list) (order integer))
    :initvals '((0 2 1 3 2 4 3 5 4 6 5 4 3 2 1 0) 1)
    :indoc '("list" "integer")
    :icon 000
    :doc "Takes a list and computes a Nth-order transition probability matrix of the elements (states) in the list. MARKOV-BUILD is meant to be used along with MARKOV-RUN.

    Arguments:
    - <data>: list of lists.
    - <order>: integer specifying the order of markov model.
    "
    (let* ()
        (setf data (List-frames data order 1))
        (let*
            (
                (thin-data (reverse (remove-dup data 'equal 1)))
                (dims (length thin-data))
                (matrix (mat-trans (append (list thin-data) (repeat-n (repeat-n 0 dims) dims)))))
            (loop for thin-item in thin-data and row from 0 to (- (length thin-data) 1) do
                (let*
                    (
                        (next-item nil)
                        (col nil)
                    )
                    (loop for item in data and pos from 0 to (- (length data) 2) do
                        (if (equal item thin-item)
                            (progn
                                (setf next-item (nth (+ pos 1) data))
                                (setf col (position next-item thin-data :test 'equal))
                                (setf (nth (+ col 1) (nth row matrix)) (+ 1 (nth (+ col 1) (nth row matrix)))))))))
            matrix)))

;--------------- Markov-run ---------------
(defmethod! Markov-run ((matrix list) (iterations integer) &optional (initial nil) (mode '1))
    :initvals '('(((0) 0 0 1 0 0) ((1) 1/2 0 0 1/2 0) ((2) 0 2/3 0 0 1/3) ((3) 0 0 1 0 0) ((4) 0 0 0 1 0)) 5 nil '1)
    :indoc '("list" "integer" "atom or list" "menu")
    :icon 000
    :menuins '((3 (("no reset" '0) ("allow reset" '1))))
    :doc "Takes a transition probability matrix and outputs a sequence of a given length. MARKOV-RUN is meant to be used along with MARKOV-BUILD.

    Arguments:

    - <matrix>: transition matrix from MARKOV-BUILD.
    - <iterations>: number of output states.
    
    &optional:
    - <initial>: initial state.
    - <mode>: output behavior.
        0 -> allow MARKOV-RUN to stop if 'dead-end' state is reached.
        1 -> if a dead-end state is reached, reset with a random state until reaching the specified number in <iterations>.
    "
    (let*
        (
            (states (car (mat-trans matrix)))
            (output nil)
            (current-state nil))
        ; initialize first state
        (if (equal initial nil)
            (setf current-state (nth-random states))
            (progn
                (setf positions (get-posn initial (car (mat-trans states))))
                (if (eq positions nil)
                    (error "
                        
                        INITIAL STATE MUST BELONG TO INPUT MATRIX.
                    ")
                    (setf current-state (nth (nth-random positions) states)))))
        (loop for rep from 0 to (- iterations 1) do
            (let*
                (
                    (row-index (position current-state states :test 'equal))
                    (weights (cdr (nth row-index matrix)))
                    (choice nil))
                (setf output (append output (list (car current-state))))
                (if (> (reduce-tree weights #'+) 0)
                    (progn 
                        (setf choice (car (pick-random states weights 1)))
                        (setf current-state choice))
                    (if (equal mode '0)
                        (progn
                            (setf rep iterations)
                            (om-print "MARKOV-RUN => Sequence reached dead end before desired length"))
                        (setf current-state (nth-random states))))))
        output))

; --------------- Optimal-sorting ---------------
(defmethod! Optimal-sorting ((st-list list) (other-lists list) (weights list))
    :initvals '((0 0 0 2) ((0 1 2 3) (2 3 4 5) (1 2 3 4)) nil)
    :indoc '("list (initial)" "list of lists" "list (optional)")
    :icon 000
    :numouts 2
    :doc "Sorts a list of lists such that the distance between adjacent lists is optimally minimized, given a starting list.
    
    NOTE: All lists must have the same length.

    Arguments:

    - <st-list>: initial list
    - <other-lists>: list of lists to sort.
    - <weights>: list of weight values for each element in <st-list>.
    
    Example:
    (optimal-sorting '(0 0 0 2) '((0 1 2 3) (2 3 4 5) (1 2 3 4)) nil) => ((0 0 0 2) (0 1 2 3) (1 2 3 4) (2 3 4 5))
    "  
    (let*
        (
            (output nil)
            (neighbors (nth-value 0 (NNS st-list other-lists weights)))
            (remaining (copy-tree neighbors))
            (positions nil))
        (loop for n in neighbors do
            (setf output (append output (list (car remaining))))
            (setf remaining (nth-value 0 (NNS (car remaining) (cdr remaining) weights))))
        (setf positions (nested-position output other-lists))
        (setf output (append (list st-list) output))
        (values-list (list output positions))))


