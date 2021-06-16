#|--------------------------------- * 
|     OM-Data library functions     |
|   [www.felipe-tovar-henao.com]    |
|               2021                |
* --------------------------------- *
|#

(in-package :om)

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
    (setq mc (om-round (x-interpolation mc b-list traj)))
    (make-instance 'chord-seq 
        :lmidic mc
        :lonset (lonset a-list)
        :ldur (ldur a-list)
        :lvel (lvel a-list)
        :loffset (loffset a-list)
        :lchan (lchan a-list)))

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
    (let*   
        (
            (output nil))
        (loop for i in input-list collect
            (setf output (append output (list
                (cond 
                    ((< i 0) (* (mod (abs i) n) (/ i (abs i))))
                    ((>= i 0) (mod i n)))))))
        output))

(defmethod! List-mod ((input-list list) (n list))
    (mapcar (lambda (input) 
        (List-mod input-list input)) n))


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


;--------------- Deep-nth ---------------
(defmethod! Deep-nth ((data list) (path list))
    :initvals '((((0 1) (2 3)) ((4 5 6) (7 8 ("nine")))) (1 1 2))
    :indoc '("list" "list")
    :icon 000
    :doc "Gets the nth element of a nested list, given a list of positions corresponding to each level in the input list."
    (let*
        (
            (out nil)
        )
        (if (eq 1 (depth path))
            (progn 
                (loop for p in path do
                    (setf data (nth p data)))
            data)
            (mapcar #'(lambda (input) (deep-nth data input)) path))))

;--------------- Deep-replace ---------------
(defmethod! Deep-replace ((data list) (path list) new)
    :initvals '((((0 1) (2 3)) ((4 5 6) (7 8 (9)))) (1 1 2) "nine")
    :indoc '("list" "list" "list or atom")
    :icon 000
    :doc "Replaces any element in a list with a new element, given a list of positions corresponding to each level in the input list."
    (let*
        (
            (data-copy (copy-tree data))
            (levels (list data-copy)))
        (loop for p in path do
            (let*
                (
                    (x (nth p data-copy)))
                (setf levels (append levels (list x)))
                (setf data-copy x)))
        (setf levels (reverse levels))
        (setf path (append (list 0) (reverse path)))
        (setf (nth 0 levels) new)
        (loop for p in path and l in levels and n from 0 to (length levels) do
            (if (and (listp l) (> n 0))
                (setf (nth p l) new)
                (setf l new))
            (setf new l))
        (car (last levels))))

;--------------- N-occurances ---------------
(defmethod! N-occurances ((x number) (l list))
    :initvals '(2 (0 1 2 1 4 3 2 2))
    :indoc '("item" "list")
    :icon 000
    :doc "Counts the number of occurances of an element in a list.

    Example:
    (n-occur 2 '(0 1 2 1 4 3 2 2)) => 3
    "
    (let*
        (
            (counter 0))
        (loop for y in l do
            (if (numberp y)
                (if (eq x y) 
                    (setf counter (+ counter 1)))
                (setf counter (+ counter (N-occurances x y)))))
        counter))

(defmethod! N-occurances ((x list) (l list))
    (let*
        (
            (counter 0))
        (loop for y in l do
            (cond 
                (
                    (eq (depth y) (depth x))
                    (if 
                        (equal x y) 
                        (setf counter (+ counter 1))))
                (
                    (> (depth y) (depth x))
                    (setf counter (+ counter (N-occurances x y))))))
        counter))

;--------------- Rep-filter ---------------
(defmethod! Rep-filter ((l list))
    :initvals '(((1 2) (3 3) (4 2) (4 2) (2 2) (2 2)))
    :indoc '("list")
    :icon 000
    :doc "Removes repetitions between consecutive elements.

    Example:
    (rep-filter '((1 2) (3 3) (4 2) (4 2) (2 2) (2 2))) => ((1 2) (3 3) (4 2) (2 2))
    "
    (let*
        (
            (output (list (car l))))
        (loop for x from 1 to (- (length l) 1) do
            (if (not (equal (nth x l) (nth (- x 1) l)))
                (setf output (append output (list (nth x l))))))
        output))

;--------------- Unique-scramble ---------------
(defmethod! Unique-scramble ((a-list list) (times integer))
    :initvals '((0 1 2) 4)
    :indoc '("list" "integer")
    :icon 000
    :doc "Performs a series of random permutations such that no element appears consecutively in the same position.

    Example:
    (unique-scramble '(0 1 2) 4) => ((0 1 2) (2 0 1) (1 2 0) (0 1 2))
    "
    (let*
        (
            (out (list a-list))
            (current (copy-tree a-list)))
        (loop for i from 1 to (- times 1) do
            (let*
                (
                    (unique nil)
                    (scrambled nil))
                (while (eq unique nil)
                    (let*
                        (
                            (dups 0))
                        (setf scrambled (permut-random current))
                        (loop for c in current and s in scrambled do
                            (if (eq c s)
                                (setf dups (+ dups 1))))
                        (if (eq dups 0)
                            (setf unique t))))
                (setf out (append out (list scrambled)))
                (setf current (copy-tree scrambled))))
        out))


