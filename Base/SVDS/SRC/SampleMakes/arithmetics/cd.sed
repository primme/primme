s/@(pre)/d/g
s/@(rpre)/d/g
s/@(rtpre)/d/g
s/@(trpre)/d/g
s/@(type)/double/g
s/@(rtype)/double/g
/#ifdefarithm L_DEFCPLX/,/#endifarithm/d
/#ifdefarithm L_DEFREAL/d
/#endifarithm/d
