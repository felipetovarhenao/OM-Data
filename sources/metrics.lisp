#|--------------------------------- * 
|     OM-Data library functions     |
|   [www.felipe-tovar-henao.com]    |
|               2021                |
* --------------------------------- *
|#

(in-package :om)

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
    (let*
        (
            (l-depth (depth b-list)))
        (if (equal weights nil)
            (setf weights (repeat-n 1.0 (length a-list))))
        (setf weights (mapcar #'(lambda (input-list) (/ input-list (list-max weights))) weights))
        (cond 
            ((eq l-depth 1) (sqrt (reduce #'+ 
                (loop for a in a-list and b in b-list and w in weights collect 
                    (* w (expt (- b a) 2))))))
            ((> l-depth 1) 
                (loop for b in b-list collect
                    (Euclidean-distance a-list b weights))))))

(defmethod! Euclidean-distance ((a number) (b number) (weights list))
    (abs (- b a)))

(defmethod! Euclidean-distance ((a number) (b-list list) (weights list))
    (loop for b in b-list collect (abs (- b a))))

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
    (let*
        (
            (mean nil)
            (st-dev nil)
            (skewness nil)
            (kurtosis nil))
        (if (eq (depth data) 1)
            (if (> (length data) 1)
                (progn 
                    (setf skewness nil) (setq kurtosis nil)
                    (setf mean (/ (reduce #'+ data) (* 1.0 (length data))))
                    (setf st-dev (sqrt 
                        (/ (reduce #'+ 
                            (loop for x in data collect
                                (expt (- x mean) 2))) 
                            (length data))))
                    (if (> (abs st-dev) 0)
                        (progn 
                            (setf skewness
                                (/ (reduce #'+
                                    (loop for x in data collect
                                        (expt (- x mean) 3))) 
                                    (* (expt st-dev 3) (length data) 1)))
                            (if (> (abs skewness) 0)
                                (setf kurtosis
                                    (/ (reduce #'+
                                        (loop for x in data collect
                                            (expt (- x mean) 4))) 
                                        (* (expt st-dev 4) (length data))))))
                        (progn (setf skewness nil) (setf kurtosis nil)))
                    (posn-match (list mean st-dev skewness kurtosis) moments))
                (posn-match (list (car data) 0 nil nil) moments))
            (loop for d in data collect (List-moments d moments)))))

;--------------- List-median ---------------
(defmethod! List-median ((l list))
    :initvals '((4 2 5 6 1 2 0 -4))
    :indoc '("list")
    :icon 000
    :doc "Computes the median of a list of numeric values"
    (let* 
        (
            (size nil)
            (index nil))
        (if (eq (depth l) 1)
            (progn 
                (setf size (length l))
                (stable-sort l #'<)
                (if (oddp size)
                    (progn
                        (setf index (nth-value 0 (floor (/ size 2))))
                        (nth index l))
                    (progn
                        (setf index (/ size 2))
                        (* (+ (nth (- index 1) l) (nth index l)) 0.5))))
            (mapcar #'(lambda (input) (List-median input)) l))))

;--------------- List-mode ---------------
(defmethod! List-mode ((l list))
    :initvals '((4 2 5 6 1 2 0 -4))
    :indoc '("list")
    :icon 000
    :doc "Computes the mode of a list of numeric values"
    (let*
        (
            (counts nil))
        (if (eq 1 (depth l))
            (progn 
                (setf counts
                    (loop for x in l collect
                        (n-occurances x l)))
                (if (eq (list-max counts) (list-min counts))
                    nil
                    (remove-dup (posn-match l (get-posn (list-max counts) counts)) 'equal 1)))
            (mapcar #'(lambda (input) (List-mode input)) l))))
    
;--------------- List-Zscore ---------------
(defmethod! List-Zscore ((x number) (l list))
    :initvals '((0 3 6) (0 1 2 3 4 5 6))
    :indoc '("list" "list")
    :icon 000
    :doc "Computes the standard score (a.k.a. z-score) value of a list based on another.

    Example: 
    (list-zscore '(0 3 6) '(0 1 2 3 4 5 6)) => (1.5 0.0 -1.5)
    "  
    (let*
        (
            (m-list (List-moments l '(0 1)))
            (mu (first m-list))
            (sigma (second m-list)))
        (/ (- mu x) sigma)))

(defmethod! List-Zscore ((x-list list) (l list))
    (mapcar #'(lambda (input) (List-Zscore input l)) x-list))

;--------------- List-covariance ---------------
(defmethod! List-covariance ((data list))
    :initvals '(((-2 2) (0 0) (1 -1) (-1 1) (-4 3) (-3 4)))
    :indoc '("list")
    :icon 000
    :doc "Computes the covariance of a given set of data points"
    (let*
        (
        (means (flat (list-moments (mat-trans data) '(0))))
        (out
            (reduce #'+
                (loop for d in data collect
                    (reduce #'* 
                        (loop for x in d and m in means collect
                            (- x m)))))))
        (/ out (- (length data) 1))))

;--------------- List-correlation ---------------
(defmethod! List-correlation ((data list))
    :initvals '(((-2 2) (0 0) (1 -1) (-1 1) (-4 3) (-3 4)))
    :indoc '("list")
    :icon 000
    :doc "Computes the correlation of a given set of data points"
    (let*
        (
            (trans-data (mat-trans data))
            (variances (loop for td in trans-data collect
                (sqrt (List-covariance (mat-trans (list td td)))))))
        (/ (List-covariance data) (reduce #'* variances))))

;--------------- Plot-points ---------------
(defmethod! Plot-points ((points list))
    :initvals '(((-2 2) (0 0) (1 -1) (-1 1) (-4 3) (-3 4)))
    :indoc '("list")
    :icon 000
    :doc "Plots a list of data points in a bpc-lib or 3dc-lib, depending on the dimensionality of the data."
    (let*
        (
            (bp-list nil)
            (dims (length (mat-trans points))))
        (loop for p in points do
            (let*
                (
                    (omp nil)
                )
                (if (equal dims 2)
                    (progn 
                        (setf omp (repeat-n (om-make-point (first p) (second p)) 2))
                        (setf bp-list (append bp-list (list (make-instance 'bpc :point-list omp :decimals 2))))))
                (if (equal dims 3)
                    (progn
                        (setf omp (repeat-n (make-3dpoint :x (first p) :y (second p) :z (third p)) 2))
                        (setf bp-list (append bp-list (list (make-instance '3dc :point-list omp :decimals 2))))))))
        (if (equal dims 2)
            (setf out (make-instance 'bpc-lib :bpf-list bp-list)))
        (if (equal dims 3)
            (setf out (make-instance '3dc-lib :bpf-list bp-list)))
        out))