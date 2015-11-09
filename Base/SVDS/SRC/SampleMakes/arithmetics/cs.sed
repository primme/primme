s/@(pre)/s/g
s/@(rpre)/s/g
s/@(rtpre)/s/g
s/@(trpre)/s/g
s/@(type)/float/g
s/@(rtype)/float/g
s/@(tpone)/+1.0e+00/g
s/@(tmone)/-1.0e+00/g
s/@(tzero)/+0.0e+00/g
s/@(rtpone)/+1.0e+00/g
s/@(rtmone)/-1.0e+00/g
s/@(rtzero)/+0.0e+00/g
s/@(TransConj)/"T"/g
/#ifdefarithm L_DEFCPLX/,/#endifarithm/d
/#ifdefarithm L_DEFREAL/d
/#endifarithm/d
