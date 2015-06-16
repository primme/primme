s/@(pre)/d/g
s/@(rpre)/d/g
s/@(rtpre)/d/g
s/@(trpre)/d/g
s/@(type)/double/g
s/@(rtype)/double/g
s/@(abs)/fabs/g
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
