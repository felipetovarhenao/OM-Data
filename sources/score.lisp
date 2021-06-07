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

;--------------- Get-transients ---------------
(defmethod! Score-filter ((self chord-seq) (ranges list) (filter-params list) (filter-types list))
    :initvals '(
        (make-instance 'chord-seq 
            :lmidic '(4800 6600 7200 8000)
            :onset '(0 300)
            :lvel '(30 50 90 127))
        '((6000 7200) (80 100)) 
        (list 'midic 'vel) (list 'pass 'reject))
    :indoc '("sequence" "list" "list (midic, pc, onset, dur, vel, offset, and/or chan)" "list (pass or reject)")
    :icon 000
    :doc "Outputs a list of onsets corresponding to detected transients in a score object." 
    (setq positions (loop for fp in filter-params collect
        (cond 
            ((equal fp 'midic) 0)
            ((equal fp 'onset) 1)
            ((equal fp 'dur) 2)            
            ((equal fp 'vel) 3)
            ((equal fp 'offset) 4)
            ((equal fp 'chan) 5)
            ((equal fp 'pc) 6))))
    (setq modes (loop for ft in filter-types collect
        (cond 
            ((equal ft 'pass) 0)
            ((equal ft 'reject) 1))))
    (setq events (get-notes self))
    (setq events (loop for e in events collect
        (append e (list (/ (nth-value 1 (om// (nth 0 e) 1200)) 100)))))
    (setq out nil)
    (loop for e in events do
        (setq out (append out (list (multi-band-filter e positions ranges modes)))))   
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

(defun multi-band-filter (in-list n-pos range mode)
    (setq out nil)
    (loop for n in n-pos and r in range and m in mode do
        (if (and (>= (nth n in-list) (list-min r)) (<= (nth n in-list) (list-max r)))
            (setq p t)
            (setq p nil))
        (if (eq m 0)
            nil
            (setq p (not p)))
        (setq out (append out (list p))))
    (if (equal (remove-dup out 'equal 1) (list t))
        in-list))

;--------------- Get-transients ---------------
(defmethod! Get-notes ((self chord-seq))
    :initvals '(nil)
    :indoc '("sequence")
    :icon 000
    :doc "" 
    (setq mc (lmidic self))
    (setq onsets (lonset self))
    (setq durs (ldur self))
    (setq vels (lvel self))
    (setq offsets (loffset self))
    (setq chans (lchan self))
    (setq out nil)
    (loop for m in mc and on in onsets and d in durs and v in vels and of in offsets and c in chans do
        (loop for subm in m and subd in d and subv in v and subof in of and subc in c do
            (setq out (append out (list (list subm on subd subv subof subc))))))
    out)
    

