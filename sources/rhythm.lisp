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
    :doc "EUCLIDEAN-RHYTHM, as the name suggests, computes the Euclidean rhythm for a given period k and number of beats/pulses n, where n < k. The output can be interpreted as a list of intervals or distances between between pulses.

    Arguments:

    - <numbeats>: number of beats in pattern.
    - <period>: period of pattern.
    - <rotation>: rotation of pattern.
    "
    (let*   
        (
            (a numbeats)
            (b (- period numbeats))
            (binary-seq (list (repeat-n 1 a) (repeat-n 0 b)))
            (val 1)
            (out nil))
        (while (> (min a b) 1)
            (setf binary-seq (mapcar #'(lambda (input) (remove nil (flat input))) (mat-trans binary-seq)))
            (setf binary-seq (group-list binary-seq (sort-list (list a b)) 'linear))
            (setf a (length (car binary-seq)))
            (setf b (length (second binary-seq))))
        (setf binary-seq (flat binary-seq))
        (loop for n in binary-seq do
            (if (eq n 0)
                (setf val (+ val 1))
                (progn
                    (setf out (append out (list val)))
                    (setf val 1))))
        (rotate (cdr (append out (list val))) rotation)))

;--------------- Rhythmicon ---------------
(defmethod! Rhythmicon ((base-dur number) (subdivisions list) (times integer) mode)
    :initvals '(3000 '(1 2 3 4 5 6 7) 4 'divisive)
    :indoc '("number" "list" "integer" "menu")
    :icon 000
    :menuins '((3 (("div" 'div) ("mult" 'mult))))
    :doc "Outputs a rhythmicon as a MULTI-SEQ, given a fundamental duration (ms), a list of subdivisions, and a number of repetitions. For each rhythmic partial, the corresponding pitch is automatically assigned, using 3300 as the fundamental.

    Arguments:

    - <base-dur>: reference duration in milliseconds.
    - <subdivisions>: values by which to either multiply or divide <base-dur>, based on <mode>.
    - <times>: number of cycles.
    - <mode>: rhythmicon mode, either multiplicative or divisive.
    "
    (let*
        (
            (onsets nil)
            (pitches nil)
            (durations nil)
            (vels (om-scale subdivisions 120 30))
            (f0 (mc->f 3300))
            (output nil)
            (velocities nil)
            (max-onset (* base-dur times))
            (period nil)
            (pitches nil))
        (if (equal mode 'mult)
            (progn 
                (setf subdivisions (om-round subdivisions))
                (setf period (list-lcm subdivisions))))
        (loop for m in subdivisions and v in vels do
            (let*
                (
                    (numnotes nil)
                    (durs nil)
                    (pre-onsets nil)
                )
                (cond
                    (
                        (equal mode 'div)
                        (progn
                            (setf numnotes (nth-value 0 (ceiling (* m times))))
                            (setf durs (repeat-n (/ base-dur m) numnotes))
                            (setf pre-onsets (om-clip (dx->x 0 durs) 0 max-onset))
                            (setf onsets (append onsets (list pre-onsets)))
                            (setf durations (append durations (list (x->dx pre-onsets))))
                            (setf pitches (append pitches (list (repeat-n (f->mc (* f0 m)) numnotes))))
                            (setf velocities (append velocities (list (repeat-n v numnotes))))))
                    (
                        (equal mode 'mult)
                        (progn
                            (setf numnotes (* (/ period m) times))
                            (setf durs (repeat-n (* base-dur m) numnotes))
                            (setf pre-onsets (dx->x 0 durs))
                            (setf onsets (append onsets (list pre-onsets)))
                            (setf durations (append durations (list (x->dx pre-onsets))))
                            (setf pitches (append pitches (list (repeat-n (f->mc (* f0 m)) numnotes))))
                            (setf velocities (append velocities (list (repeat-n v numnotes)))))))))
        (make-instance 'multi-seq :chord-seqs (reverse (loop for o in onsets and p in pitches and d in durations and v in velocities collect 
            (make-instance 'chord-seq 
                :lmidic p 
                :lonset o 
                :ldur d 
                :lvel v))))))

; --------------- Risset rhythm---------------
(defun riss-onsets (onsets totdur period v mode)
    (let*
        (
            (totdur (list-max onsets))
            (v-exp (expt 2 v))
            (tl nil)
            (te))
        (loop for x from 0 to (- v-exp 1) do
            (setf tl (append tl (om+ (butlast onsets) (om* totdur x)))))
        (setf tl (om/ (append tl (list (* v-exp totdur))) (* v-exp totdur)))
        (cond 
            (
                (equal mode 0)
                (setf te (om-log (om+ 1 tl) 2)))
            (
                (equal mode 1)
                (setf te (om- (om^ 2 tl) 1))))
        (om* te period)))

(defun riss-rates (te totdur period v mode)
    (let*
        (
            (normalized (om/ te period))
            (minval nil)
            (maxval nil))
        (if (equal mode 0)
            (setf normalized (om^ 2 (om+ v normalized)))
            (progn 
                (setf minval (expt 2 v))
                (setf maxval (expt 2 (+ v 1)))
                (setf normalized (om-scale (om* (om+ (om-log (om+ normalized 1) 2) 1) minval) maxval minval))))
        (om* (/ (* totdur (log 2)) period) normalized)))

(defun riss-amps (te period voices v mode)
    (let* 
        (
            (index-list (om/ te period)))
        (if (equal mode 1)
            (setf index-list (om- 1 index-list)))
        (om+ 0.5 (om* -0.5 (om-cos (om* (* pi 2) (om+ (/ v voices) (om/ index-list voices))))))))

;--------------- Risset-rhyhtm ---------------
(defmethod! Risset-rhythm ((self chord-seq) (speed number) (voices integer) (rep integer) &optional (onset-mode '0) (mc-mode '2))
    :initvals '((make-instance 'chord-seq :lmidic '(6000 6200 6400 6700) :lonset '(0 250) :ldur '(250)) 2.5 4 4 '0 '2)
    :indoc '("chord-seq" "number" "integer" "integer" "menu")
    :icon 000
    :menuins '((4 (("accelerando" '0) ("ritardando" '1))) (5 (("preserve pitch" '0) ("period-wise transp." '1) ("chord-wise transp." '2))))
    :doc "Builds a seemingly eternal accelerando/ritardando, based on an input CHORD-SEQ, a speed factor, number of voices, number of repetitions.

    Arguments:

    - <self>: a CHORD-SEQ.
    - <speed>: speed factor of central voice.
    - <voices>: number of voices in resulting rhythm.
    - <rep>: number of cycles.

    &optional:
    - <onset-mode>: agogic mode, either accelerando or ritardando.
    - <mc-mode>: mode for midicents treatment, either preserve original values in <self>, transpose pitch period-wise, or chord-wise.
    "
    (let*
        (
            (out nil)
            (mc-list (lmidic self))
            (onsets (lonset self))
            (durs (ldur self))
            (vels (lvel self))
            (offsets (loffset self))
            (chans (lchan self))
            (leg (legato self))
            (totdur (list-max onsets))
            (period (* totdur speed))
            (fq-list (mc->f mc-list)))
        (loop for v from 0 to (- voices 1) do
            (let*
                (
                    (v-exp (expt 2 v))
                    (new-onsets (riss-onsets onsets totdur period v onset-mode))
                    (rates (riss-rates new-onsets totdur period v onset-mode))
                    (amps (riss-amps new-onsets period voices v onset-mode))
                    (new-mc nil)
                    (new-vels nil)
                    (new-durs nil)
                    (te nil)
                    (seq nil))
                (loop for f in (flat (repeat-n fq-list v-exp) 1) and v in (flat (repeat-n vels v-exp) 1) and d in (flat (repeat-n durs v-exp) 1) and a in (butlast amps) and r in (butlast rates) and n from 0 to (- (length rates) 1) do
                    (let*
                        (
                            (r-index nil))
                        (setf new-durs (append new-durs (list (om/ d r))))
                        (if (equal mc-mode 1)
                            (progn
                                (setf r-index (* (nth-value 0 (floor (/ n (length durs)))) (length durs)))
                                (setf r (nth r-index rates))))
                        (if (> mc-mode 0)
                            (setf f (om* r f)))
                        (setf new-mc (append new-mc (list (f->mc f))))
                        (setf new-vels (append new-vels (list (om* v a))))))
                (setf te (flat 
                    (loop for x from 0 to (- rep 1) collect
                        (om+ (butlast new-onsets) (* x period))
                    ) 1))
                (setf te (append te (list (* rep (list-max new-onsets)))))
                (setf seq (make-instance 'chord-seq 
                    :lmidic (flat (repeat-n new-mc rep) 1) 
                    :lonset te 
                    :ldur (flat (repeat-n new-durs rep) 1)  
                    :lvel (flat (repeat-n new-vels rep) 1)
                    :loffset (flat (repeat-n offsets (* rep v-exp)) 1)
                    :lchan (flat (repeat-n chans (* rep v-exp)) 1)
                    :legato leg))
                (setf out (append out (list seq)))))
        (reverse out)))
