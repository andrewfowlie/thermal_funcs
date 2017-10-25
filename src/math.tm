:Begin:
:Function:      J_F
:Pattern:       JF[ysq_?NumericQ, derivative_Integer: 0]
:Arguments:     {ysq, derivative}
:ArgumentTypes: {Real, Integer}
:ReturnType:    Real
:End:

:Begin:
:Function:      J_B
:Pattern:       JB[ysq_?NumericQ, derivative_Integer: 0]
:Arguments:     {ysq, derivative}
:ArgumentTypes: {Real, Integer}
:ReturnType:    Real
:End:

:Evaluate:      JF::usage = "JF[ysq, derivative: 0] returns the fermionic thermal function, or its
                             first or second derivative for derivative = 1 and 2, respectively."
                             
:Evaluate:      JB::usage = "JB[ysq, derivative: 0] returns the bosonic thermal function, or its
                             first or second derivative for derivative = 1 and 2, respectively."
