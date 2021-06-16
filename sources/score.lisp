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
    (setq cents nil) (setq onsets nil) (setq durations nil) (setq velocities nil) (setq offsets nil) (setq chans nil)
    (setq seq-onsets (lonset self))
    (loop for time-pt in time-pt-list do
        (setq end-time-pt (+ time-pt samp-dur))

        (loop for s-onset in seq-onsets and s-dur in (ldur self) and i from 0 to (- (length seq-onsets) 1) do

            (setq end-pts (om+ s-onset s-dur))
            (setq st-pts (om- end-pts s-dur))
            (setq time-points (mat-trans (list st-pts end-pts)))

            (setq current_chord nil)
            (setq current_dur nil)
            (setq current_vel nil)
            (setq current_offset nil)
            (setq current_chan nil)

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
                        (setq maxdur (- end-time-pt (first tp)))
                        (setq out-dur (nth j (nth i (ldur self))))
                        (if (or (eq clip-mode 2) (eq clip-mode 3))
                            (setq out-dur (min out-dur maxdur)))

                        (setq current_chord (append current_chord (list (nth j (nth i (lmidic self))))))
                        (setq current_dur (append current_dur (list out-dur)))
                        (setq current_vel (append current_vel (list (nth j (nth i (lvel self))))))
                        (setq current_chan (append current_chan (list (nth j (nth i (lchan self))))))
                        (setq current_offset (append current_offset (list (nth j (nth i (loffset self)))))))))
            (if (not (eq current_chord nil))
                (progn 
                    (setq cents (append cents (list current_chord)))
                    (setq out-onset (nth i (lonset self)))
                    (if (or (eq clip-mode '1) (eq clip-mode '3))
                        (setq out-onset (max out-onset time-pt)))
                    (setq onsets (append onsets (list out-onset)))
                    (setq durations (append durations (list current_dur)))
                    (setq velocities (append velocities (list current_vel)))
                    (setq chans (append chans (list current_chan)))
                    (setq offsets (append offsets (list current_offset)))))))
    (setq min-onset (list-min (flat onsets)))
    (make-instance 'chord-seq :lmidic cents :lonset (om- onsets min-onset) :ldur durations :lvel velocities :loffset offsets :lchan chans))

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
    (setq channels (lchan self))
    (setq cents nil) (setq onsets nil) (setq durations nil) (setq velocities nil) (setq offsets nil)
    (loop for chans in channels and i from 0 to (- (length channels) 1) do
        (setq current_chord nil)
        (setq current_dur nil)
        (setq current_vel nil)
        (setq current_offset nil)
        (loop for ch in chans and j from 0 to (- (length chans) 1) do
            (if (eq ch midi-chan)
                (progn
                    (setq current_chord (append current_chord (list (nth j (nth i (lmidic self))))))
                    (setq current_dur (append current_dur (list (nth j (nth i (ldur self))))))
                    (setq current_vel (append current_vel (list (nth j (nth i (lvel self))))))
                    (setq current_offset (append current_offset (list (nth j (nth i (loffset self)))))))))
        (if (not (eq current_chord nil))
            (progn 
                (setq cents (append cents (list current_chord)))
                (setq onsets (append onsets (list (nth i (lonset self)))))
                (setq durations (append durations (list current_dur)))
                (setq velocities (append velocities (list current_vel)))
                (setq offsets (append offsets (list current_offset))))))
    (make-instance 'chord-seq :lmidic cents :lonset onsets :ldur durations :lvel velocities :loffset offsets :lchan midi-chan))

(defmethod! Extract-channel ((self multi-seq) (midi-chan integer))
    (setq chord-seq-list (inside self))
    (make-instance 'multi-seq :chord-seqs (mapcar #'(lambda (input) (Extract-channel input midi-chan)) chord-seq-list)))

;--------------- Multi-join---------------
(defmethod! Multi-join ((seqs list) &optional (mode '0) (concat-offset nil))
    :initvals '(nil '0 nil)
    :indoc '("multi-seq or poly" "join mode" "list")
    :icon 000
    :doc "Joins a list of score objects, either by merging or concatenating them." 
    :menuins '((1 (("concat" '0) ("merge" '1))))
    (if (and (equal mode '0) (equal concat-offset nil))
        (setq concat-offset 
            (cdr (dx->x 0 (loop for s in seqs collect
                (list-max (lonset s)))))))
    (setq out (car seqs))
    (setq seq-list (cdr seqs))
    (loop for s in seq-list and i from 0 to (- (length seq-list) 1) do
        (if (equal mode '0)
            (setq out (concat out s (nth i concat-offset)))
            (setq out (merger out s))))
    out)

(defmethod! Multi-join ((self multi-seq) &optional (mode '0) (concat-offset nil))
    (setq seqs (chord-seqs self))
    (Multi-join seqs mode concat-offset))

(defmethod! Multi-join ((self poly) &optional (mode '0) (concat-offset nil))
    (setq seqs (voices self))
    (Multi-join seqs mode concat-offset))

;--------------- Get-transients ---------------
(defmethod! Get-transients ((self chord-seq) (threshold number))
    :initvals '(nil 0.05)
    :indoc '("sequence" "number")
    :icon 000
    :doc "Outputs a list of onsets corresponding to detected transients in a score object." 
    (setq energy (list 0))
    (setq velocities (lvel self))
    (setq min-vel (list-min velocities))
    (setq velocities (append (list (list min-vel) (list min-vel)) velocities))
    (setq onsets (lonset self))
    (setq out nil)
    (loop for v in velocities and i from 0 to (- (length velocities) 1) do
        (setq vel (log (reduce #'+ v)))
        (if (> i 0)
            (progn 
                (setq val (max 0 (- vel pvel)))
                (if (> val threshold)
                    (setq out (append out (list (nth (- i 2) onsets)))))))
        (setq pvel vel))
    out)

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
    (setq positions (loop for fp in filter-params collect
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
    (setq modes (loop for ft in filter-types collect
        (cond 
            ((equal ft 'pass) 0)
            ((equal ft 'reject) 1))))
    (setq bool (position bool (list 'and 'or)))
    (setq events (chord-seq->event-list self))
    (setq events (loop for e in events collect
        (append e (list (/ (nth-value 1 (om// (nth 0 e) 1200)) 100)))))
    (setq out nil)
    (loop for e in events do
        (setq out (append out (list (multi-band-filter e positions ranges modes bool)))))   
    (setq out (remove nil out))
    (setq out (mat-trans out))
    (make-instance 'chord-seq
        :lmidic (first out)
        :lonset (second out)
        :ldur (third out)
        :lvel (fourth out)
        :loffset (fifth out)
        :lchan (sixth out)
        :legato (legato self)))

(defun multi-band-filter (in-list n-pos range mode bool)
    (setq out nil)
    (loop for n in n-pos and r in range and m in mode do
        (if (and (>= (nth n in-list) (list-min r)) (<= (nth n in-list) (list-max r)))
            (setq p t)
            (setq p nil))
        (if (eq m 0)
            nil
            (setq p (not p)))
        (setq out (append out (list p))))
    (cond
        (
            (equal bool 0)
            (if (equal (remove-dup out 'equal 1) (list t))
                in-list))
        (
            (equal bool 1)
            (if (not (equal nil (position t out)))
                in-list))))

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
    (setq mc (lmidic self))
    (setq onsets (lonset self))
    (setq durs (ldur self))
    (setq vels (lvel self))
    (setq offsets (loffset self))
    (setq chans (lchan self))
    (setq dxonsets (x->dx onsets))
    (setq out nil)
    (loop for m in mc and on in onsets and d in durs and v in vels and of in offsets and c in chans and dxo in dxonsets do
        (loop for subm in m and subd in d and subv in v and subof in of and subc in c do
            (setq out (append out (list (list subm on subd subv subof subc dxo (length m)))))))
    out)

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
    (setq cents (lmidic self))
    (setq onsets (mat-trans (list (lonset self))))
    (setq durs (ldur self))
    (setq vels (lvel self))
    (setq chans (lchan self))
    (setq writehead 0)
    (setq tblsize 1000)

    (setq write-incr-var (om-abs write-incr-var))
    (setq read-incr-var (om-abs write-incr-var))
    (setq vel-var (om-abs vel-var))
    (setq dur-var (om-abs dur-var))

    (if (atom grain-dur)
        (setq durtable (repeat-n grain-dur tblsize))
        (setq durtable (nth-value 2 (om-sample grain-dur tblsize))))
    (if (atom write-incr)
        (setq writetbl (repeat-n write-incr tblsize))
        (setq writetbl (nth-value 2 (om-sample write-incr tblsize))))
    (if (atom read-incr)
        (setq readtbl (repeat-n read-incr tblsize))
        (setq readtbl (nth-value 2 (om-sample read-incr tblsize)))) 
    (if (atom dur-var)
        (setq durvartable (repeat-n dur-var tblsize))
        (setq durvartable (nth-value 2 (om-sample dur-var tblsize))))    
    (if (atom write-incr-var)
        (setq wincrvartable (repeat-n write-incr-var tblsize))
        (setq wincrvartable (nth-value 2 (om-sample write-incr-var tblsize))))   
    (if (atom read-incr-var)
        (setq rincrvartable (repeat-n read-incr-var tblsize))
        (setq rincrvartable (nth-value 2 (om-sample read-incr-var tblsize))))  
    (if (atom vel-var)
        (setq velvartable (repeat-n vel-var tblsize))
        (setq velvartable (om-clip (nth-value 2 (om-sample vel-var tblsize)) 0 127)))
    
    (setq velvartable (om/ velvartable 127.0))
    (setq event-posn nil)
    (setq gr-durs nil)
    (setq gr-onsets nil)
    (setq gr-vels nil)
    (setq readhead st-read-time)

    (setq input-dur (list-max onsets))
    (while (<= writehead output-dur)
        (setq index (om-round (* (/ writehead output-dur) (- tblsize 1))))
        (setq dur (nth index durtable))
        (setq r-incr (nth index readtbl))
        (setq w-incr (nth index writetbl))
        (setq wvar (nth index wincrvartable))
        (setq rvar (nth index rincrvartable))
        (setq vvar (nth index velvartable))
        (setq dvar (nth index durvartable))

        (setq wvar (om-random (* -1 wvar) wvar))
        (setq rvar (om-random (* -1 rvar) rvar))
        (setq dvar (om-random (* -1 dvar) dvar))
        (setq vvar (om-random (- 1 vvar) 1))
        
        (setq evpos (car (nth-value 1 (NNS (list readhead) onsets nil))))
        (setq event-posn (append event-posn (list evpos)))
        (setq gr-onsets (append gr-onsets (list writehead)))
        (setq gr-durs (append gr-durs (list (om-clip (nth evpos durs) 10 (max 10 (+ dur dvar))))))
        (setq gr-vels (append gr-vels (list (om-abs (om* (nth evpos vels) vvar)))))

        (setq readhead (+ readhead r-incr rvar))
        (if (and (> readhead input-dur) (equal read-mode 'circular))
            (setq readhead (nth-value 1 (om// readhead input-dur))))

        (setq writehead (om-clip (+ writehead w-incr wvar) 0 nil)))

    (make-instance 'chord-seq
        :lmidic (posn-match cents event-posn)
        :lonset (om-clip gr-onsets 0 nil)
        :ldur gr-durs
        :lvel gr-vels
        :lchan (posn-match chans event-posn)))

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
    (setq time-pts (list-frames time-pts 2))
    (setq low-bound (- (list-min (lmidic self)) 600))
    (setq hi-bound (+ (list-max (lmidic self)) 600))
    (if octave-eq?
        (setq chords (mapcar #'(lambda (input) (fill-range input low-bound hi-bound)) chords)))
    (setq new-chords nil)
    (loop for mc in (lmidic self) and onset in (lonset self) do
        (setq new-mc (copy-tree mc))
        (loop for tp in time-pts and ch in chords do
            (if (and (>= onset (first tp)) (< onset (second tp)))
                (setq new-mc (list-quantize mc ch accuracy))))
        (setq new-chords (append new-chords (list new-mc))))
    (make-instance 'chord-seq 
        :lmidic new-chords
        :lonset (lonset self)
        :ldur (ldur self)
        :lvel (lvel self)
        :loffset (loffset self)
        :lchan (lchan self)
        :legato (legato self))
)