#|--------------------------------- * 
|     OM-Data library functions     |
|   [www.felipe-tovar-henao.com]    |
|               2021                |
* --------------------------------- *
|#

(in-package :om)

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
    (if (equal mode 'mult)
        (setq period (list-lcm subdivisions))
    )
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
    :doc "Builds a seemingly eternal accelerando/ritardando, based on an input CHORD-SEQ, a speed factor, number of voices, number of repetitions."

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
