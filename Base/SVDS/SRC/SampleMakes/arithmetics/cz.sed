s/@(pre)/z/g
s/@(rpre)/d/g
s/@(rtpre)/dz/g
s/@(trpre)/zd/g
s/@(type)/complex double/g
s/@(rtype)/double/g
/#ifdefarithm L_DEFREAL/,/#endifarithm/d
/#ifdefarithm L_DEFCPLX/d
/#endifarithm/d
