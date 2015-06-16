s/@(pre)/z/g
s/@(rpre)/d/g
s/@(rtpre)/dz/g
s/@(trpre)/zd/g
s/@(type)/Complex_Z/g
s/@(rtype)/double/g
s/@(abs)/z_abs_primme/g
s/@(tpone)/{+1.0e+00,+0.0e00}/g
s/@(tmone)/{-1.0e+00,+0.0e00}/g
s/@(tzero)/{+0.0e+00,+0.0e00}/g
s/@(rtpone)/+1.0e+00/g
s/@(rtmone)/-1.0e+00/g
s/@(rtzero)/+0.0e+00/g
s/@(TransConj)/"C"/g
/#ifdefarithm L_DEFREAL/,/#endifarithm/d
/#ifdefarithm L_DEFCPLX/d
/#endifarithm/d
