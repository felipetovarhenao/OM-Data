#|--------------------------------- * 
|     OM-Data library functions     |
|   [www.felipe-tovar-henao.com]    |
|               2021                |
* --------------------------------- *
|#

(in-package :om)

; --------------- string-rewrite ---------------
(defun string-rewrite (axiom rules iterations)
    (let*
        (
            (axiom (flat (list axiom) 1)))
        (loop for n from 0 to (- iterations 1) do
            (loop for a in axiom and i from 0 to (- (length axiom) 1) do
                (loop for r in rules do
                    (if (equal a (first r))
                        (progn
                            (setf (nth i axiom) (second r))))))
            (setf axiom (flat axiom 1)))
        axiom))

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

; --------------- 2D-Turtle ---------------
(defmethod! 2D-Turtle ((lsys list) (mag-rules list) (theta-rules list) (memory-rules list) &optional (theta 0))
    :initvals '(([ f - f ] + f [ f - f ] + f [ f - f ] + f [ f - f ] + f [ f - f ] + f [ f - f ] + f) '((f 1)) '((+ 60) (- -60)) '(([ 1) (] 0)) 0)
    :indoc '("list" "list" "list" "list" "number")
    :icon 000
    :doc "2D Turtle graphics"
    (let*
        (
            (x 0)
            (y 0)
            (mag 0)
            (memory nil)
            (out (list (om-make-point x y)))
            (state nil))
        (loop for s in lsys do
            (loop for mr in mag-rules do
                (if (equal s (first mr))
                    (progn
                        (setf mag (second mr))
                        (setf x (+ x (* mag (cos (deg->rad theta)))))
                        (setf y (+ y (* mag (sin (deg->rad theta)))))
                        (setf out (append out (list (om-make-point x y)))))))
            (loop for tr in theta-rules do
                (if (equal s (first tr))
                    (progn
                        (setf theta (+ theta (second tr))))))
            (loop for mr in memory-rules do
                (if (equal s (first mr))
                    (progn 
                        (if (equal (second mr) 1)
                            (setf memory (append memory (list (list x y theta mag))))
                        )
                        (if (equal (second mr) 0)
                            (progn
                                (setf state (car (last memory)))
                                (setf x (first state))
                                (setf y (second state))
                                (setf theta (third state))
                                (setf mag (fourth state))
                                (setf memory (butlast memory))))))))
        (make-instance 'bpc :point-list out)))

;--------------- Make-sieve ---------------
(defmethod! Make-sieve ((list list) (reps integer) sieve-mode sieve-type &optional (offset '0))
    :initvals '((2 3) 1 'union 'nil 0)
    :indoc '("list" "integer" "menu" "menu" "number")
    :menuins '(
        (2 (("union" 'union) ("diff" 'diff)))
        (3 (("nil" 'nil) ("complement" 'complement))))
    :icon 000
    :doc "Builds N full periods of a sieve, based on a list of integers. Make-sieve is meant to be a compact version of OM's native CRIBLE class and functions"
    (let* ()
        (setf list (remove 0 list))
        (let*
            (
                (period (+ offset (* reps (list-lcm list))))
                (sieves (loop for l in list collect
                    (arithm-ser offset period l)))
            )
            (cond
                (
                    (equal sieve-mode 'union)
                    (setf out (list-union sieves)))
                (
                    (equal sieve-mode 'diff)
                    (setf out (list-diff sieves))))
            (if (equal sieve-type 'complement)
                (setf out (list-diff (list (arithm-ser offset period 1) (flat (list out))))))
            (stable-sort out #'<))))

(defun list-union (list)
    (let*
        (
            (out (car list)))
        (loop for l in (cdr list) do
            (setf out (x-union out l)))
        out))

(defun list-intersect (list)
    (let*
        (
            (out (car list)))
        (loop for l in (cdr list) do
            (setf out (x-intersect out l)))
        out))

(defun list-diff (list)
    (let*
        (
            (out (car list)))
        (loop for l in (cdr list) do
            (setf out (x-diff out l)))
        out))

(defun list-lcm (list)
    (let*
        (
            (out (car list)))
        (loop for l in (cdr list) do
            (setf out (lcm out l)))
        out))

(defun list-gcd (list)
    (let*
        (
            (out (car list)))
        (loop for l in (cdr list) do
            (setf out (gcd out l)))
        out))

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
    (let*   
        (
            (out nil)
            (mod-n (reduce #'+ seq)))
        (setf seq (nth-value 1 (om// (append seq (list (car seq))) mod-n)))
        (loop for x from 1 to n-tiers do
            (let*
                (
                    (diff-seq nil))
                (loop for i from 0 to (- (length seq) 2) do
                    (let*
                        (
                            (current (nth i seq))
                            (next (nth (+ i 1) seq))
                            (val nil))
                        (if (>= next current)
                            (setf val (abs (- next current)))
                            (setf val (abs (- (+ next mod-n) current))))
                        (setf diff-seq (append diff-seq (list val)))))
                (setf out (append out (list diff-seq)))
                (setf seq (append diff-seq (list (car diff-seq))))))
        out))

;--------------- Recursivus ---------------
(defmethod! Recursivus ((list-a list) (list-b list) &optional (rotations nil))
    :initvals '((7300 7500 7200 6600 6800 6400 6700 7400) (7400 8200 7700 8000 6900 7100 6500 6300 6600 8500) nil)
    :indoc '("list" "number" "menu")
    :menuins '((2 (("nil" nil) ("rotations" 'rotations))))
    :icon 000
    :doc "Given two numeric lists A and B, outputs a list of possible extensions of A through segments B. More precisely, it searches for segments of B, where the sum of changes in those segments match the change (i.e. interval) between adjacent elements in A. Based on the technique described by Antonio de Sousa Dias in 'Free transposition of Audio Processing Techniques into the Domain of Notes'
    "
    (let*
        (
            (dx-a (x->dx (append list-a (last list-a))))
            (numrot (length list-b))
            (out (mat-trans (list (mat-trans (list list-a)))))
            (forms nil))
        (if (equal rotations nil)
            (setf forms (list (x->dx list-b)))
            (setf forms 
                (loop for r from 0 to (- numrot 1) collect
                    (x->dx (rotate list-b r)))))
        (loop for dx-b in forms do
            (loop for fn from 1 to numrot do
                (let*
                    (
                        (b-seg (first-n dx-b fn))
                        (dx-seg (reduce #'+ b-seg))
                        (mel-seg nil))
                    (loop for a in dx-a and x from 0 to (length dx-a) do
                        (if (equal dx-seg a)
                            (progn
                                (setf mel-seg (dx->x (nth x list-a) b-seg))
                                (if (> (length mel-seg) 2)
                                    (setf mel-seg (butlast mel-seg)))
                                (if (equal nil (position mel-seg (nth x out) :test 'equal))
                                    (setf out (subs-posn out x (append (nth x out) (list mel-seg)))))))))))
        out))