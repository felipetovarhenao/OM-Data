#|--------------------------------- * 
|     OM-Data library functions     |
|   [www.felipe-tovar-henao.com]    |
|               2021                |
* --------------------------------- *
|#

(in-package :om ) 

(defclass! criblex ()
  ((cr-exp :initform '(2 0 18) :initarg :cr-exp :accessor cr-exp)
   (maxn :initform 0 :accessor maxn))
  (:icon 423))

(defmethod initialize-instance ((self criblex)  &rest rest)
    (declare (ignore rest))
    (call-next-method)
    (if 
        (not (stringp (car (cr-exp self))))
        (setf (maxn self) (third (cr-exp self)))
        (setf (maxn self) (apply 'max (mapcar 'maxn (cdr (cr-exp self)))))
    )
    self)

(defmethod! c-union ((self criblex) (criblex criblex) &rest rest)
    :doc "It makes the set-theoretical union of two sieves or cribles" 
    :icon 423
    (let* 
        (
            (c-list (append (list self criblex) rest))
            (newc 
                (make-instance 'criblex
                    :cr-exp (cons "criblex-u" c-list)
                )
            )   
        )
    newc)
)

(defmethod! c-intersection ((self criblex) (criblex criblex) &rest rest)
    :doc "It makes the set-theoretical intersection of two sieves or cribles"
    :icon 423
    (let* 
        (
            (c-list (append (list self criblex) rest))
            (newc 
                (make-instance 'criblex
                    :cr-exp (cons "criblex-i" c-list)
                )
            )
        )
    newc)
)

(defmethod! c-complement ((self criblex))
  :doc "It makes the set-theoretical complement of a sieve or criblex with respect to the total chromatic universe"
  :icon 423
    (let* 
        (
            (newc 
                (make-instance 'criblex
                    :cr-exp (list "criblex-c" self)
                )
            )
        )
    newc)
)

(defun criblex-u (&rest rest)
  (let (rep)
    (loop for item in rest do
          (setf rep (union rep item)))
    rep))

(defun criblex-i (&rest rest)
  (let ((rep (car rest)))
    (loop for item in rest do
          (setf rep (intersection rep item)))
    rep))

(defun criblex-c (criblex)
  (let (total )
    (loop for item in rest do
          (setf rep (union rep item)))
    rep))


(defmethod! revel-criblex ((self criblex))
  :doc "It makes visible the criblex object by giving its corresponding residual classes values. You can connect the revel-criblex object to an n-cercle in order to see a possible circular representation of a given sieve.  "
  :icon 423
  (cond
   ((not (stringp (car (cr-exp self))))
    (loop for i from (second (cr-exp self)) to (third (cr-exp self)) by (first (cr-exp self))
          collect i))
   ((string-equal (car (cr-exp self)) "criblex-c")
    (let ((total (loop for i from 0 to (maxn (second (cr-exp self)))  collect i)))
      (sort (set-difference total (revel-criblex  (second (cr-exp self)))) '<)))
   (t
    (let ((function (interne (string (car (cr-exp self)))))
          (args (cdr (cr-exp self))))
      (sort (apply function (loop for item in  args collect (revel-criblex item))) '<))) ))


(defmethod! revel-criblex ((self list))
  (mapcar 'revel-criblex self))

(defmethod! revel-criblex ((self t))
  nil)

