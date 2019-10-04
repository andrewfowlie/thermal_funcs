:Begin:
:Function:      J_F
:Pattern:       JF[ysq_?NumericQ, OptionsPattern[{derivative -> 0,
                absError -> 1.*^-7, relError -> 1.*^-7, maxN -> 10000,
                fast -> 0}]]
:Arguments:     {ysq, OptionValue[derivative], N[OptionValue[absError]],
                N[OptionValue[relError]], OptionValue[maxN], OptionValue[fast]}
:ArgumentTypes: {Real, Integer, Real, Real, Integer, Integer}
:ReturnType:    Real
:End:

:Begin:
:Function:      J_B
:Pattern:       JB[ysq_?NumericQ, OptionsPattern[{derivative -> 0,
                absError -> 1.*^-7, relError -> 1.*^-7, maxN -> 10000,
                fast -> 0}]]
:Arguments:     {ysq, OptionValue[derivative], N[OptionValue[absError]],
                N[OptionValue[relError]], OptionValue[maxN], OptionValue[fast]}
:ArgumentTypes: {Real, Integer, Real, Real, Integer, Integer}
:ReturnType:    Real
:End:

:Evaluate:      JF::usage = "JF[ysq] returns the bosonic thermal function, or its
                             first or second derivative for derivative -> 1 and 2, respectively."

:Evaluate:      JB::usage = "JB[ysq] returns the bosonic thermal function, or its
                             first or second derivative for derivative -> 1 and 2, respectively."

:Evaluate:      Derivative[1][JB][ysq_] := JB[ysq, derivative -> 1]
:Evaluate:      Derivative[1][JF][ysq_] := JF[ysq, derivative -> 1]
:Evaluate:      Derivative[2][JB][ysq_] := JB[ysq, derivative -> 2]
:Evaluate:      Derivative[2][JF][ysq_] := JF[ysq, derivative -> 2]
:Evaluate:      Derivative[1][JF][ysq_, derivative->1] := JF[ysq, derivative->2]
:Evaluate:      Derivative[1][JB][ysq_, derivative->1] := JB[ysq, derivative->2]
