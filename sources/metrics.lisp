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
    (make-instance 'bpf :point-list histo :decimals 1))

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