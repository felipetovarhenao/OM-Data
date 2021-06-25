#|--------------------------------- * 
|     OM-Data library functions     |
|   [www.felipe-tovar-henao.com]    |
|               2021                |
* --------------------------------- *
|#

(in-package :om)

; -------------- (UNLISTED) FUNCTIONS -------------------
(defun get-posn (n source)
    (let*
        (
            (i 0)
            (out nil))
        (while (< i (length source))
            (if (equal (nth i source) n)
                (setf out (flat (list out i))))
            (setf i (+ i 1)))
        (remove nil out)))

(defun closest (a b-list)
    (let*
        (
            (distances (loop for b in b-list and n from 0 to (- (length b-list) 1) collect 
                (list (abs (- b a)) b n)))            )
        (stable-sort distances #'< :key #'first)
        (setf distances (car distances))
        (list (second distances) (third distances))))

(defun depth (list)
    (if (listp list)
        (+ 1 (reduce #'max (mapcar #'depth list) :initial-value 0)) 0))

(defun mix (a b f)
    (if (and (atom a) (atom b))
        (+ (* a (- 1 f)) (* b f))
        (loop for x in a and y in b collect
            (mix x y f))))

(defun nested-mix (a b f)
    (if (or (eq (+ (depth a) (depth b)) 2) (eq (atom a) (atom b)))
        (mix a b f)
        (mapcar #'(lambda (in1 in2) (nested-mix in1 in2 f)) a b)))

(defun dec->bin (x n)
    (let*
        (
            (out nil)
            (vals nil))
        (if (> x 0)
            (progn
                (while (> x 0)
                    (setf vals (multiple-value-list (floor x 2)))
                    (setf x (car vals))
                    (setf out (append out (cdr vals))))
                (if (< (length out) n)
                    (setf out (append out (repeat-n 0 (- n (length out)))))))
            (setf out (repeat-n 0 n)))
        out))

(defun bin->dec (seq)
    (let*
        (
            (out 0))
        (loop for x in seq and i from 0 to (- (length seq) 1) do
            (if (equal x 1)
                (setf out (+ out (expt 2 i)))))
        out))

(defun branch-posn (data nodes)
    (bin->dec (loop for d in data and n in nodes collect
        (if (< d n) 0 1))))

(defun wedge-sum (a b)
    (if (eq (mod b 2) 0)
        (setf b (* -1 (nth-value 0 (floor b 2))))
        (setf b (nth-value 0 (floor (+ b 1) 2))))   
    (+ a b))

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

#| 
    TODO:
        - Granulate
        - PCA
        - Chroma vector with durations and velocities
        - Rhythmic-distribution
 |#