#|--------------------------------- * 
|     OM-Data library functions     |
|   [www.felipe-tovar-henao.com]    |
|               2021                |
* --------------------------------- *
|#

(in-package :om)

;--------------- Segment-seq ---------------
(defmethod! Segment-seq ((self chord-seq) (time-pt-list list) (samp-dur number) &optional (detection-mode '0) (clip-mode '0))
    :initvals '((make-instance 'chord-seq :lmidic '((6000 7200) (6400 7000) 6700) :lchan '((1 1) (2 1) 3)) (0 1000) 250 '0 '0)
    :indoc '("chord-seq or multi-seq" "list or number" "number" "menu" "menu")
    :icon 000
    :menuins '((3 (("detect onsets" '0) ("detect onsets and durations" '1))) (4 (("no clipping" '0) ("clip onsets" '1) ("clip durations" '2) ("clip onsets and durations" '3))))
    :doc "Extracts a segment from a CHORD-SEQ, given a list of time points and duration for all segments.
    "
    (let*
        (
            (cents nil) 
            (onsets nil) 
            (durations nil) 
            (velocities nil) 
            (offsets nil) 
            (chans nil)
            (seq-onsets (lonset self)))
        (loop for time-pt in time-pt-list do
            (let*
                (
                    (end-time-pt (+ time-pt samp-dur))
                    (out-onset nil))
                (loop for s-onset in seq-onsets and s-dur in (ldur self) and i from 0 to (- (length seq-onsets) 1) do
                    (let*
                        (
                            (end-pts (om+ s-onset s-dur))
                            (st-pts (om- end-pts s-dur))
                            (time-points (mat-trans (list st-pts end-pts)))
                            (current_chord nil)
                            (current_dur nil)
                            (current_vel nil)
                            (current_offset nil)
                            (current_chan nil))
                        (loop for tp in time-points and j from 0 to (- (length time-points) 1) do
                            (if 
                                (if (eq detection-mode '1)
                                    (or
                                        (and ; detect onset within range
                                            (>= (first tp) time-pt) ; onset after st
                                            (< (first tp) end-time-pt)) ; onset before end
                                        (and ; detect ending within range
                                            (> (second tp) time-pt) ; outset after st
                                            (<= (second tp) end-time-pt)) ; and outset before end
                                        (and  ; detect note during range
                                            (<= (first tp) time-pt) ; onset before st
                                            (>= (second tp) end-time-pt) ; and onset after end
                                        )
                                    )
                                    (and ; detect onset within range only
                                        (>= (first tp) time-pt)
                                        (< (first tp) end-time-pt)))
                                (progn 
                                    (let*
                                        (
                                            (maxdur (- end-time-pt (first tp)))
                                            (out-dur (nth j (nth i (ldur self)))))
                                        (if (or (eq clip-mode 2) (eq clip-mode 3))
                                            (setf out-dur (min out-dur maxdur)))
                                        (setf current_chord (append current_chord (list (nth j (nth i (lmidic self))))))
                                        (setf current_dur (append current_dur (list out-dur)))
                                        (setf current_vel (append current_vel (list (nth j (nth i (lvel self))))))
                                        (setf current_chan (append current_chan (list (nth j (nth i (lchan self))))))
                                        (setf current_offset (append current_offset (list (nth j (nth i (loffset self))))))))))
                        (if (not (eq current_chord nil))
                            (progn 
                                (setf cents (append cents (list current_chord)))
                                (setf out-onset (nth i (lonset self)))
                                (if (or (eq clip-mode '1) (eq clip-mode '3))
                                    (setf out-onset (max out-onset time-pt)))
                                (setf onsets (append onsets (list out-onset)))
                                (setf durations (append durations (list current_dur)))
                                (setf velocities (append velocities (list current_vel)))
                                (setf chans (append chans (list current_chan)))
                                (setf offsets (append offsets (list current_offset)))))))))
        (make-instance 'chord-seq :lmidic cents :lonset (om- onsets (list-min (flat onsets))) :ldur durations :lvel velocities :loffset offsets :lchan chans)))

(defmethod! Segment-seq ((self chord-seq) (time-pt number) (samp-dur number) &optional (detection-mode '0) (clip-mode '0))
    (Segment-seq self (list time-pt) samp-dur detection-mode clip-mode))

(defmethod! Segment-seq ((self multi-seq) (time-pt number) (samp-dur number) &optional (detection-mode '0) (clip-mode '0))
    (make-instance 'multi-seq :chord-seqs 
        (loop for seq in (inside self) collect
            (Segment-seq seq (list time-pt) samp-dur detection-mode clip-mode))))

(defmethod! Segment-seq ((self multi-seq) (time-pt list) (samp-dur number) &optional (detection-mode '0) (clip-mode '0))
    (make-instance 'multi-seq :chord-seqs 
        (loop for seq in (inside self) collect
            (Segment-seq seq time-pt samp-dur detection-mode clip-mode))))

;--------------- Extract-channel ---------------
(defmethod! Extract-channel ((self chord-seq) (midi-chan integer))
    :initvals '((make-instance 'chord-seq :lmidic '((6000 7200) (6400 7000) 6700) :lchan '((1 1) (2 1) 3)) 1)
    :indoc '("chord-seq or multi-seq" "integer")
    :icon 000
    :doc "
        Extracts the notes with the specified midi-channel from a CHORD-SEQ or MULTI-SEQ 
    "
    (let*
        (
            (channels (lchan self))
            (cents nil) 
            (onsets nil) 
            (durations nil) 
            (velocities nil) 
            (offsets nil))
        (loop for chans in channels and i from 0 to (- (length channels) 1) do
            (let*
                (
                    (current_chord nil)
                    (current_dur nil)
                    (current_vel nil)
                    (current_offset nil))
                (loop for ch in chans and j from 0 to (- (length chans) 1) do
                    (if (eq ch midi-chan)
                        (progn
                            (setf current_chord (append current_chord (list (nth j (nth i (lmidic self))))))
                            (setf current_dur (append current_dur (list (nth j (nth i (ldur self))))))
                            (setf current_vel (append current_vel (list (nth j (nth i (lvel self))))))
                            (setf current_offset (append current_offset (list (nth j (nth i (loffset self)))))))))
                (if (not (eq current_chord nil))
                    (progn 
                        (setf cents (append cents (list current_chord)))
                        (setf onsets (append onsets (list (nth i (lonset self)))))
                        (setf durations (append durations (list current_dur)))
                        (setf velocities (append velocities (list current_vel)))
                        (setf offsets (append offsets (list current_offset)))))))
        (make-instance 'chord-seq :lmidic cents :lonset onsets :ldur durations :lvel velocities :loffset offsets :lchan midi-chan)))

(defmethod! Extract-channel ((self multi-seq) (midi-chan integer))
    (let*
        (
            (chord-seq-list (inside self)))
        (make-instance 'multi-seq :chord-seqs (mapcar #'(lambda (input) (Extract-channel input midi-chan)) chord-seq-list))))

;--------------- Multi-join---------------
(defmethod! Multi-join ((seqs list) &optional (mode '0) (concat-offset nil))
    :initvals '(nil '0 nil)
    :indoc '("multi-seq or poly" "join mode" "list")
    :icon 000
    :doc "Joins a list of score objects, either by merging or concatenating them." 
    :menuins '((1 (("concat" '0) ("merge" '1))))
    (let*
        (
            (out nil)
            (seq-list nil))
        (if (and (equal mode '0) (equal concat-offset nil))
            (setf concat-offset 
                (cdr (dx->x 0 (loop for s in seqs collect
                    (list-max (lonset s)))))))
        (setf out (car seqs))
        (setf seq-list (cdr seqs))
        (loop for s in seq-list and i from 0 to (- (length seq-list) 1) do
            (if (equal mode '0)
                (setf out (concat out s (nth i concat-offset)))
                (setf out (merger out s))))
        out))

(defmethod! Multi-join ((self multi-seq) &optional (mode '0) (concat-offset nil))
    (Multi-join (chord-seqs self) mode concat-offset))

(defmethod! Multi-join ((self poly) &optional (mode '0) (concat-offset nil))
    (Multi-join (voices self) mode concat-offset))

;--------------- Get-transients ---------------
(defmethod! Get-transients ((self chord-seq) (threshold number))
    :initvals '(nil 0.05)
    :indoc '("sequence" "number")
    :icon 000
    :doc "Outputs a list of onsets corresponding to detected transients in a score object." 
    (let*
        (
            (out nil)
            (energy (list 0))
            (velocities (lvel self))
            (min-vel (list-min velocities))
            (onsets (lonset self))
            (pvel nil))
        (setf velocities (append (list (list min-vel) (list min-vel)) velocities))
        (loop for v in velocities and i from 0 to (- (length velocities) 1) do
            (let* 
                (
                    (vel (log (reduce #'+ v)))
                    (val nil))
                (if (> i 0)
                    (progn 
                        (setf val (max 0 (- vel pvel)))
                        (if (> val threshold)
                            (setf out (append out (list (nth (- i 2) onsets)))))))
                (setf pvel vel)))
        out))

(defmethod! Get-transients ((self multi-seq) (threshold number))
    (Get-transients (Multi-join self 1) threshold))

;--------------- Score-filter ---------------
(defmethod! Score-filter ((self chord-seq) (ranges list) (filter-params list) (filter-types list) bool)
    :initvals '(
        (make-instance 'chord-seq 
            :lmidic '(4800 6600 7200 8000)
            :onset '(0 300)
            :lvel '(30 50 90 127))
        '((6000 7200) (80 100)) 
        (list 'midic 'vel) (list 'pass 'reject) 'and)
    :indoc '("sequence" "list" "list (midic, pc, onset, dur, vel, offset, and/or chan)" "list (pass or reject)" "menu")
    :menuins '((4 (("and" 'and) ("or" 'or))))
    :icon 000
    :doc "Allows to filter out events within a CHORD-SEQ, based on any number of conditions. Given a list of ranges, filter parameters (see list below), filter types (pass/reject), and a boolean operator (and/or), each event in the CHORD-SEQ evaluated.
    Filter params:
	- onset
	- midic
	- dur
	- vel
	- offset
	- chan
	- dxonset (interonset distance)
	- chordsize
	- pc (pitch class)
" 
    (let*
        (
            (positions (loop for fp in filter-params collect
                (cond 
                    ((equal fp 'midic) 0)
                    ((equal fp 'onset) 1)
                    ((equal fp 'dur) 2)            
                    ((equal fp 'vel) 3)
                    ((equal fp 'offset) 4)
                    ((equal fp 'chan) 5)
                    ((equal fp 'dxonset) 6)
                    ((equal fp 'chordsize) 7)
                    ((equal fp 'pc) 8))))
            (modes (loop for ft in filter-types collect
                (cond 
                    ((equal ft 'pass) 0)
                    ((equal ft 'reject) 1))))
            (bool (position bool (list 'and 'or)))
            (events (chord-seq->event-list self))
            (out nil))
        (setf events (loop for e in events collect
            (append e (list (/ (nth-value 1 (om// (nth 0 e) 1200)) 100)))))
        (loop for e in events do
            (setf out (append out (list (multi-band-filter e positions ranges modes bool)))))   
        (setf out (remove nil out))
        (setf out (mat-trans out))
        (make-instance 'chord-seq
            :lmidic (first out)
            :lonset (second out)
            :ldur (third out)
            :lvel (fourth out)
            :loffset (fifth out)
            :lchan (sixth out)
            :legato (legato self))))

(defun multi-band-filter (in-list n-pos range mode bool)
    (let*
        (
            (out nil))
        (loop for n in n-pos and r in range and m in mode do
            (let*
                (
                    (p nil))
                (if (and (>= (nth n in-list) (list-min r)) (<= (nth n in-list) (list-max r)))
                    (setf p t)
                    (setf p nil))
                (if (eq m 0)
                    nil
                    (setf p (not p)))
                (setf out (append out (list p)))))
        (cond
            (
                (equal bool 0)
                (if (equal (remove-dup out 'equal 1) (list t))
                    in-list))
            (
                (equal bool 1)
                (if (not (equal nil (position t out)))
                    in-list)))))

;--------------- Chord-seq->event-list ---------------
(defmethod! Chord-seq->event-list ((self chord-seq))
    :initvals '(nil)
    :indoc '("sequence")
    :icon 000
    :doc "Extracts event information from a CHORD-SEQ in the form of individual lists per note. Each list is formatted as follows:
       (MIDICENT
        ONSET
        DURATION
        VELOCITY
        OFFSET
        MIDI-CHANNEL
        INTER-ONSET_DISTANCE
        CHORD-SIZE)
    " 
    (let*
        (
            (mc (lmidic self))
            (onsets (lonset self))
            (durs (ldur self))
            (vels (lvel self))
            (offsets (loffset self))
            (chans (lchan self))
            (dxonsets (x->dx onsets))
            (out nil))
        (loop for m in mc and on in onsets and d in durs and v in vels and of in offsets and c in chans and dxo in dxonsets do
            (loop for subm in m and subd in d and subv in v and subof in of and subc in c do
                (setf out (append out (list (list subm on subd subv subof subc dxo (length m)))))))
        out))

;--------------- Granulate ---------------
(defmethod! Granulate ((self chord-seq) (output-dur number) grain-dur write-incr read-incr &optional (st-read-time '0) (dur-var '0) (write-incr-var '0) (read-incr-var '0) (vel-var '0) (read-mode 'linear))
    :initvals '(
        (make-instance 'chord-seq
            :lmidic '((6000 6400 6700) (6000 6500 6900) (5900 6200 6500 6700) (6000 6400 6700))
            :lonset '(0 1000 1500 2000)
            :ldur '(1000 500 500 1000))
        5000
        (500 200 100)
        (250 50)
        (50 200)
        0
        0
        0
        0
        0
        'circular)
    :menuins '((8 (("linear" 'linear) ("circular" 'circular))))
    :indoc '("sequence" "number (ms)" "number or list (ms)" "number or list (ms)" "number or list (ms)" "number (ms)" "number or list (ms)" "number or list (ms)" "number or list (ms)" "menu")
    :icon 000
    :doc "Performs symbolic granulation on a CHORD-SEQ.
        " 
    (let*
        (
            (cents (lmidic self))
            (onsets (mat-trans (list (lonset self))))
            (durs (ldur self))
            (vels (lvel self))
            (chans (lchan self))
            (writehead 0)
            (tblsize 1000)
            (write-incr-var (om-abs write-incr-var))
            (read-incr-var (om-abs write-incr-var))
            (durtbl nil)
            (writetbl nil)
            (readtbl nil)
            (durvartable nil)
            (wincrvartable nil)
            (rincrvartable nil)
            (velvartable nil)
            (event-posn nil)
            (gr-durs nil)
            (gr-onsets nil)
            (gr-vels nil)
            (readhead st-read-time)
            (input-dur (list-max onsets)))
        (setf vel-var (om-abs vel-var))
        (setf dur-var (om-abs dur-var))
        (if (atom grain-dur)
            (setf durtbl (repeat-n grain-dur tblsize))
            (setf durtbl (nth-value 2 (om-sample grain-dur tblsize))))
        (if (atom write-incr)
            (setf writetbl (repeat-n write-incr tblsize))
            (setf writetbl (nth-value 2 (om-sample write-incr tblsize))))
        (if (atom read-incr)
            (setf readtbl (repeat-n read-incr tblsize))
            (setf readtbl (nth-value 2 (om-sample read-incr tblsize)))) 
        (if (atom dur-var)
            (setf durvartable (repeat-n dur-var tblsize))
            (setf durvartable (nth-value 2 (om-sample dur-var tblsize))))    
        (if (atom write-incr-var)
            (setf wincrvartable (repeat-n write-incr-var tblsize))
            (setf wincrvartable (nth-value 2 (om-sample write-incr-var tblsize))))   
        (if (atom read-incr-var)
            (setf rincrvartable (repeat-n read-incr-var tblsize))
            (setf rincrvartable (nth-value 2 (om-sample read-incr-var tblsize))))  
        (if (atom vel-var)
            (setf velvartable (repeat-n vel-var tblsize))
            (setf velvartable (om-clip (nth-value 2 (om-sample vel-var tblsize)) 0 127)))
        (setf velvartable (om/ velvartable 127.0))
        (while (<= writehead output-dur)
            (let*
                (
                    (index (om-round (* (/ writehead output-dur) (- tblsize 1))))
                    (dur (nth index durtbl))
                    (r-incr (nth index readtbl))
                    (w-incr (nth index writetbl))
                    (wvar (nth index wincrvartable))
                    (rvar (nth index rincrvartable))
                    (vvar (nth index velvartable))
                    (dvar (nth index durvartable))
                    (evpos (car (nth-value 1 (NNS (list readhead) onsets nil)))))
                (setf event-posn (append event-posn (list evpos)))
                (setf wvar (om-random (* -1 wvar) wvar))
                (setf rvar (om-random (* -1 rvar) rvar))
                (setf dvar (om-random (* -1 dvar) dvar))
                (setf vvar (om-random (- 1 vvar) 1))
                (setf gr-onsets (append gr-onsets (list writehead)))
                (setf gr-durs (append gr-durs (list (om-clip (nth evpos durs) 10 (max 10 (+ dur dvar))))))
                (setf gr-vels (append gr-vels (list (om-abs (om* (nth evpos vels) vvar)))))
                (setf readhead (+ readhead r-incr rvar))
                (if (and (> readhead input-dur) (equal read-mode 'circular))
                    (setf readhead (nth-value 1 (om// readhead input-dur))))
                (setf writehead (om-clip (+ writehead w-incr wvar) 0 nil))))
        (make-instance 'chord-seq
            :lmidic (posn-match cents event-posn)
            :lonset (om-clip gr-onsets 0 nil)
            :ldur gr-durs
            :lvel gr-vels
            :lchan (posn-match chans event-posn))))

;--------------- Chord-seq->event-list ---------------
(defmethod! Harmonic-mapping ((self chord-seq) (chords list) (time-pts list) &optional (accuracy 1.0) (octave-eq? t))
    :initvals '(
        (make-instance 'chord-seq 
            :lmidic '((6550 6900 7275) (6300) (6100) (6000 6500 7300) (6550 6900 7275) (6300) (6000 6500 7300))
            :lonset '(0 250 500 1000 1500 2000 2500 3000)
            :ldur '(250))
        ((6000 6400 6700) (5900 6200 6500 6700) (6000 6400 6800)) 
        (0 1000 2000 3000)
        1.0
        t)
    :indoc '("sequence" "list" "list" "number" "menu")
    :icon 000
    :menuins '((4 (("nil" nil) ("octave equivalence" t))))
    :doc "Given a list of time points/markers (ms) and a list of target pitch collections/chords corresponding to those markers, pitches in the input chord-seq between those time points are snapped to the corresponding chords.
    " 
    (let*
        (
            (time-pts (list-frames time-pts 2))
            (low-bound (- (list-min (lmidic self)) 600))
            (hi-bound (+ (list-max (lmidic self)) 600))
            (new-chords nil))
        (if octave-eq?
            (setf chords (mapcar #'(lambda (input) (fill-range input low-bound hi-bound)) chords)))
        (loop for mc in (lmidic self) and onset in (lonset self) do
            (let*
                (
                    (new-mc (copy-tree mc)))
                (loop for tp in time-pts and ch in chords do
                    (if (and (>= onset (first tp)) (< onset (second tp)))
                        (setf new-mc (list-quantize mc ch accuracy))))
                (setf new-chords (append new-chords (list new-mc)))))
        (make-instance 'chord-seq 
            :lmidic new-chords
            :lonset (lonset self)
            :ldur (ldur self)
            :lvel (lvel self)
            :loffset (loffset self)
            :lchan (lchan self)
            :legato (legato self))))