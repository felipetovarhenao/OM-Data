#|--------------------------------- * 
|     OM-Data library functions     |
|   [www.felipe-tovar-henao.com]    |
|               2021                |
* --------------------------------- *
|#

(in-package :om)

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
    (let* ()
        (if (and (eq (depth mc) 1) (> (length mc) 1))
            (progn 
                (let* 
                    (
                        (icv (repeat-n 0 7))
                        (times 0)
                        (pc-list (om-round (om/ (List-mod mc 1200) 100) 1)))
                    (loop for pc in pc-list do   
                        (loop for x from times to (- (length mc) 1) do
                            (let*
                                (
                                    (interval (abs (- (nth x pc-list) pc)))
                                    (posn (List-fold (nth-value 0 (round interval)) 0 6)))
                                (setf (nth posn icv) (+ (nth posn icv) 1))))
                        (setf times (+ times 1)))
                    (cdr icv)))
            (mapcar #'(lambda (input) (IC-vector input)) mc))))

;--------------- Mc-clip ---------------
(defmethod! Mc-clip ((mc number) (lower-bound number) (upper-bound number))
    :initvals '(5500 6000 7200)
    :indoc '("list" "number" "number")
    :icon 000
    :doc "Constrains a list of of midicents to a given range, using octave equivalence.
    
    Example: 
    (mc-clip 5500 6000 7200) => 6700
    "
    (let* 
        (
            (lowdif (- 1200 (mod (abs (- lower-bound mc)) 1200)))
            (hidif (- 1200 (mod (abs (- upper-bound mc)) 1200)))
            (out mc))
        (cond
            (
                (< mc lower-bound)
                (setf out (+ lower-bound lowdif)))
            (
                (>= mc upper-bound)
                (setf out (- upper-bound hidif))))
        out))

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
    (let*
        (
            (range (abs (- upper-bound lower-bound)))
            (oct-offset (* 1200 (nth-value 0 (floor lower-bound 1200))))
            (oct-range (* 1200 (nth-value 0 (ceiling range 1200))))
            (oct-range-floor (* 1200 (nth-value 0 (floor range 1200)))))
        (setf mc (+ oct-offset (abs (nth-value 1 (om// (- mc oct-offset) oct-range)))))
        (cond 
            (  
                (< mc lower-bound)
                (setf mc (+ mc oct-range-floor)))
            (
                (> mc upper-bound)
                (setf mc (- mc oct-range-floor))))
        mc))

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
    (let*
        (
            (out nil))  
        (loop for x from 0 to (- n-times 1) collect
            (setf out (append out (list (+ mc-start (* x ic))))))
        (mc-wrap out lower-bound upper-bound)))

(defmethod! Ic-cycle ((mc-start number) (ic list) (lower-bound number) (upper-bound number) (n-times number))
    (let*
        (
            (out (list mc-start)))
        (loop for x from 0 to (- n-times 2) collect
            (if (eq (depth ic) 1)
                (progn 
                    (setf mc-start (+ mc-start (nth (mod x (length ic)) ic)))
                    (setf out (append out (list mc-start))))))
        (mc-wrap out lower-bound upper-bound)))

(defmethod! Harmonic-distr ((mc-list list) (mc-fund number) (max-dist integer) &optional (mode '0))
    :initvals '((6000 6300 6700 7200) 4500 200 0)
    :indoc '("list" "number" "integer" "menu")
    :menuins '((3 (("nil" 0) ("incl. fundamental" 1))))
    :icon 000
    :doc "Re-distributes a list of midicents in register, such that it resembles a harmonic series of a given (virtual) fundamental.
    "
    (let* ()   
        (if (< max-dist 15)
            (setf max-dist 15))
        (let*
            (
                (f0 (mc->f mc-fund))
                (pcs (nth-value 1 (om// mc-list 1200)))
                (p-index 1)
                (out nil))
            (while (> (length pcs) 0)
                (let*
                    (
                        (partial (f->mc (* f0 p-index)))
                        (pool (fill-range pcs (- partial 600) (+ partial 600)))
                        (closest-elem (car (closest partial pool)))
                        (dist (abs (- partial closest-elem))))
                    (if (< dist max-dist)
                        (progn 
                            (let* 
                                (
                                    (pc (nth-value 1 (om// closest-elem 1200))))
                                (setf pcs (remove pc pcs :count 1))
                                (setf out (append out (list closest-elem))))))
                    (setf p-index (+ p-index 1))))
            (if (equal mode '1)
                (setf out (append (list mc-fund) out)))
            out)))

;--------------- Voice-leading functions ---------------
(defun connect-chords (chord-a chord-b punish-intervals)
    (let*
        (
            (pc-a (nth-value 1 (om// chord-a 1200)))
            (pc-b (nth-value 1 (om// chord-b 1200)))
            (chord-a (sort-list chord-a))
            (chord-b (sort-list chord-b))
            (dxmaster nil))
        (let*
            (
                (same? (and (equal nil (x-diff pc-a pc-b)) (equal nil (x-diff pc-b pc-a)))))
            (if same?
                chord-a
                (progn 
                    (loop for a in pc-a and x from (- (length pc-a) 1) do
                        (let* 
                            (
                                (dx nil))
                            (loop for b in pc-b and y from (- (length pc-b) 1) do
                                (let*
                                    (
                                        (dx1 (- b a))
                                        (dx2 (- (+ b 1200) a))
                                        (dx3 (- (- b 1200) a))
                                        (dxlist (list dx1 dx2 dx3))
                                        (absdx (om-abs dxlist))
                                        (mindx (list-min absdx)))
                                    (setf dx (append dx (list (posn-match dxlist (get-posn mindx absdx)))))))
                            (setf dxmaster (append dxmaster (list (flat dx))))))
                    (let* 
                        (
                            (transformations (apply #'combinations dxmaster))
                            (options nil))
                        (loop for tr in transformations do
                            (let* 
                                (
                                    (candidate (sort-list (om+ chord-a tr)))
                                    (pc-cand (nth-value 1 (om// candidate 1200)))
                                    (member? nil)
                                    (edit-cost 0)
                                    (parallel-cost 0)
                                    (mov-cost (* 0.5 (motion-cost tr)))
                                    (spread-cost 0)
                                    (opt nil)
                                )
                                (if (equal nil (x-diff pc-b pc-cand))
                                    (progn 
                                        (setf edit-cost (ghisi-edit-distance candidate chord-a))
                                        (setf spread-cost (* 1/3 (span-cost candidate chord-a)))
                                        (setf parallel-cost (* 1/2 (reduce #'+ (parallelisms chord-a candidate punish-intervals))))
                                        (setf opt (list candidate (+ spread-cost mov-cost parallel-cost edit-cost)))
                                        (setf member? (equal nil (member opt options :test 'equal)))
                                        (if member?
                                            (setf options (append options (list opt))))))))
                        (stable-sort options #'< :key 'second)
                        (caar (mat-trans options))))))))

(defun combinations (&rest lists)
    (if (endp lists)
        (list nil)
        (mapcan 
            (lambda (inner-val) 
                (mapcar (lambda (outer-val) (cons outer-val inner-val))
                    (car lists)))
            (apply #'combinations (cdr lists)))))

(defun get-duplicates (l)
    (let* 
        (
            (thin-l (remove-dup l 'equal 1))
            (output nil))
        (loop for x in thin-l do
            (let* 
                (
                    (elems (posn-match l (get-posn x l)))
                )
                (if (> (length elems) 1)
                    (setf output (append output (list elems))))))
        output))

(defun ghisi-edit-distance (a b)
    (let* 
        (
            (pc-a (remove-dup (nth-value 1 (om// a 1200)) 'equal 1))
            (pc-b (remove-dup (nth-value 1 (om// b 1200)) 'equal 1))
            (common (length (x-intersect a b)))
            (largest (max (length b) (length b)))
            (cost-a (* 1/12 (- largest common)))
            (pc-a-size (length pc-a))
            (pc-b-size (length pc-b))
            (cost-b (- (max pc-a-size pc-b-size) (length (x-intersect pc-a pc-b)))))
        (+ cost-a cost-b)))

(defun parallelisms (chord-a chord-b mcdx)
    (let* ()
        (if (atom mcdx)
            (setf mcdx (list mcdx)))
        (let* 
            (
                (voices (mat-trans (list chord-a chord-b)))
                (numvoices (length voices))
                (output (repeat-n 0 (length mcdx))))
            (loop for x from 0 to (- numvoices 1) do
                (let* 
                    (  
                        (a (nth x voices)))   
                    (loop for y from (+ x 1) to (- numvoices 1) do    
                        (let* 
                            (
                                (b (nth y voices))
                                (dx (nth-value 1 (om// (om- b a) 1200)))
                                (parallel? (equal (first dx) (second dx)))
                            )
                            (if parallel?
                                (loop for z in mcdx and k from 0 to (- (length mcdx) 1) do
                                    (if (equal (first dx) z)
                                        (setf (nth k output) (+ (nth k output) 1))
                                    )
                                )
                            )
                        )
                    )
                )
            )
            output)
    )
)

(defun motion-cost (dx)
    (let*
        (
            (netdx (/ (abs (reduce #'+ dx)) 100))
            (common (n-occurances 0 dx)))
        (/ netdx (max common 1))))

(defun span-cost (a b)
    (let*
        (
            (span-a (- (list-max a) (list-min a)))
            (span-b (- (list-max b) (list-min b))))
        (/ (abs (- span-a span-b)) 1200)))

;--------------- Voice-leading methods ---------------
(defmethod! Voice-leading ((st-chord list) (other-chords list) &optional (punish-intervals '(0 700)))
    :initvals '((4800 5500 6000 6400 7100) ((200 500 900) (0 500 900) (200 500 700 1100) (0 400)) (0 700))
    :indoc '("list" "list" "list")
    :icon 000
    :doc "Voice leading
    "
    (let* 
        (
            (a (sort-list st-chord))
            (output nil))
        (loop for b in other-chords do
            (let* 
                (
                    (size-a (length a))
                    (size-b (length b))
                    (largest (max size-a size-b))
                    (new nil))
                (setf b (sort-list b))
                (cond
                    (
                        (> size-a size-b)
                        (setf b (iterate-list (remove nil (flat (mat-trans (pc-group b)))) largest)))
                    (
                        (> size-b size-a)
                        (setf a (iterate-list (remove nil (flat (mat-trans (pc-group a)))) largest))))
                (setf new (connect-chords a b punish-intervals))
                (setf output (append output (list (remove-dup new 'equal 1))))
                (setf a new)))
        (append (list st-chord) output)))

(defun pc-group (l)
    (let*
        (
            (pc-list (nth-value 1 (om// l 1200)))
            (unique-pc (remove-dup pc-list 'equal 1))
            (out (repeat-n nil (length unique-pc)))
            (pos nil))
        (loop for p in l and pc in pc-list do
            (setf pos (position pc unique-pc :test 'equal))
            (setf (nth pos out) (append (nth pos out) (list p))))
        out))




