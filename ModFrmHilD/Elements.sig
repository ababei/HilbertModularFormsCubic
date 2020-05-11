174,0
T,ModFrmHilDElt,-,0
A,ModFrmHilDElt,4,Parent,Weight,Level,Coefficients
S,PercentM,returns a string to produce f in a magma session,0,1,0,0,0,0,0,0,0,ModFrmHilDElt,,298,-38,-38,-38,-38,-38
S,Print,,0,2,0,0,1,0,0,0,0,298,,0,0,ModFrmHilDElt,,-38,-38,-38,-38,-38,-38
S,in,,0,2,0,0,0,0,0,0,0,ModFrmHilDElt,,0,0,-1,,36,-38,-38,-38,-38,-38
S,IsCoercible,,0,2,0,0,0,0,0,0,0,-1,,0,0,ModFrmHilDElt,,36,-1,-38,-38,-38,-38
S,eq,"compares Parent, Weight, and Coefficients",0,2,0,0,0,0,0,0,0,ModFrmHilDElt,,0,0,ModFrmHilDElt,,36,-38,-38,-38,-38,-38
S,Parent,returns ModFrmHilD space where f lives,0,1,0,0,0,0,0,0,0,ModFrmHilDElt,,ModFrmHilD,-38,-38,-38,-38,-38
S,Weight,returns weight of f,0,1,0,0,0,0,0,0,0,ModFrmHilDElt,,82,-38,-38,-38,-38,-38
S,Coefficients,returns coefficients of f as SeqEnum,0,1,0,0,0,0,0,0,0,ModFrmHilDElt,,82,-38,-38,-38,-38,-38
S,CoefficientsParent,returns coefficients of f as SeqEnum,0,1,0,0,0,0,0,0,0,ModFrmHilDElt,,82,-38,-38,-38,-38,-38
S,BaseField,returns base field of parent of f,0,1,0,0,0,0,0,0,0,ModFrmHilDElt,,27,-38,-38,-38,-38,-38
S,Level,returns level of parent of f,0,1,0,0,0,0,0,0,0,ModFrmHilDElt,,217,-38,-38,-38,-38,-38
S,Precision,returns precision of parent of f,0,1,0,0,0,0,0,0,0,ModFrmHilDElt,,148,-38,-38,-38,-38,-38
S,Ideals,returns ideals indexing f up to bound on the norm,0,1,0,0,0,0,0,0,0,ModFrmHilDElt,,82,-38,-38,-38,-38,-38
S,Dictionary,returns dictionary of (parent of) f,0,1,0,0,0,0,0,0,0,ModFrmHilDElt,,457,-38,-38,-38,-38,-38
S,ModFrmHilDEltInitialize,Create an empty ModFrmHilDElt object,0,0,0,0,0,0,0,ModFrmHilDElt,-38,-38,-38,-38,-38
S,ModFrmHilDEltCopy,new instance of ModFrmHilDElt,0,1,0,0,0,0,0,0,0,ModFrmHilDElt,,ModFrmHilDElt,-38,-38,-38,-38,-38
S,HMF,"Given M and a SeqEnum of coefficients, return ModFrmHilDElt with parent M. WARNING: user is responsible for coefficients and order...proceed with caution",1,2,1,82,0,148,4,0,0,0,0,0,0,0,82,,0,0,82,,0,0,217,,0,0,ModFrmHilD,,ModFrmHilDElt,-38,-38,-38,-38,-38
S,HMF,"given a ModFrmHilD, a weight k, and an associative array of coefficients on the ideals, return ModFrmHilDElt",1,2,1,82,0,148,4,0,0,0,0,0,0,0,457,,0,0,82,,0,0,217,,0,0,ModFrmHilD,,ModFrmHilDElt,-38,-38,-38,-38,-38
S,HMFZero,create zero ModHilFrmDElt of weight k,1,2,1,82,0,148,3,0,0,0,0,0,0,0,82,,0,0,217,,0,0,ModFrmHilD,,ModFrmHilDElt,-38,-38,-38,-38,-38
S,GetCoefficient,returns a_I,0,2,0,0,0,0,0,0,0,217,,0,0,ModFrmHilDElt,,-45,-38,-38,-38,-38,-38
S,SetCoefficient,sets a_I to c,0,3,0,0,1,0,0,0,0,-45,,0,0,217,,0,0,ModFrmHilDElt,,-38,-38,-38,-38,-38,-38
S,GetCoefficient,returns a_nu,0,2,0,0,0,0,0,0,0,216,,0,0,ModFrmHilDElt,,-45,-38,-38,-38,-38,-38
S,SetCoefficient,sets a_nu to c,0,3,0,0,1,0,0,0,0,-45,,0,0,216,,0,0,ModFrmHilDElt,,-38,-38,-38,-38,-38,-38
S,GetCoefficientsIdeal,"given a ModFrmHilDElt, return associative array of coefficients on the ideals",0,1,0,0,0,0,0,0,0,ModFrmHilDElt,,457,-38,-38,-38,-38,-38
S,GetCoefficientsRepresentatives,"given a ModFrmHilDElt, return associative array of coefficients on the representatives",0,1,0,0,0,0,0,0,0,ModFrmHilDElt,,457,-38,-38,-38,-38,-38
S,EigenformToHMF,Construct the ModFrmHilDElt in M determined (on prime ideals up to norm prec) by hecke_eigenvalues,1,2,1,82,0,148,4,0,0,0,0,0,0,0,457,,0,0,82,,0,0,217,,0,0,ModFrmHilD,,ModFrmHilDElt,-38,-38,-38,-38,-38
S,NewformsToHMF,returns Hilbert newforms,1,2,1,82,0,148,3,0,0,0,0,0,0,0,82,,0,0,217,,0,0,ModFrmHilD,,82,-38,-38,-38,-38,-38
S,GaloisOrbit,returns the full Galois orbit of a modular form,0,1,0,0,0,0,0,0,0,ModFrmHilDElt,,82,-38,-38,-38,-38,-38
S,GaloisOrbitDescent,returns the full Galois orbit of a modular form over Q,0,1,0,0,0,0,0,0,0,ModFrmHilDElt,,82,-38,-38,-38,-38,-38
S,CuspFormBasis,returns a basis for M of weight k,1,2,1,82,0,148,3,0,0,0,0,0,0,0,82,,0,0,217,,0,0,ModFrmHilD,,82,148,-38,-38,-38,-38
S,HeckeOperator,returns T(n)(f) for the trivial character,0,2,0,0,0,0,0,0,0,217,,0,0,ModFrmHilDElt,,ModFrmHilDElt,-38,-38,-38,-38,-38
S,HeckeOperator,"returns T(n)(f) with loss of precision. If a basis for the space is provided, returns T(n)(f) without loss of precision",0,3,0,0,0,0,0,0,0,GrpHeckeElt,,0,0,217,,0,0,ModFrmHilDElt,,ModFrmHilDElt,-38,-38,-38,-38,-38
S,EisensteinSeries,"Let aa*bb be the modulus of psi*eta^-1. Return the Eisenstein series E_k(eta,psi) in M_k(aa*bb,eta*psi)",1,4,1,82,0,148,5,0,0,0,0,0,0,0,82,,0,0,GrpHeckeElt,,0,0,GrpHeckeElt,,0,0,217,,0,0,ModFrmHilD,,ModFrmHilDElt,-38,-38,-38,-38,-38
S,LValue_Recognized,This is a toolbox function to compute L values in the right space,1,3,1,82,0,148,4,0,0,0,0,0,0,0,82,,0,0,GrpHeckeElt,,0,0,217,,0,0,ModFrmHilD,,28,-38,-38,-38,-38,-38
S,Basis,returns a Basis for the space,1,2,1,82,0,148,3,0,0,0,0,0,0,0,82,,0,0,217,,0,0,ModFrmHilD,,82,148,-38,-38,-38,-38
S,*,scale f by integer c,0,2,0,0,0,0,0,0,0,ModFrmHilDElt,,0,0,148,,ModFrmHilDElt,-38,-38,-38,-38,-38
S,*,scale f by some scalar c,0,2,0,0,0,0,0,0,0,ModFrmHilDElt,,0,0,-1,,ModFrmHilDElt,-38,-38,-38,-38,-38
S,+,return f+g,0,2,0,0,0,0,0,0,0,ModFrmHilDElt,,0,0,ModFrmHilDElt,,ModFrmHilDElt,-38,-38,-38,-38,-38
S,-,return f-g,0,2,0,0,0,0,0,0,0,ModFrmHilDElt,,0,0,ModFrmHilDElt,,ModFrmHilDElt,-38,-38,-38,-38,-38
S,*,return f*g,0,2,0,0,0,0,0,0,0,ModFrmHilDElt,,0,0,ModFrmHilDElt,,ModFrmHilDElt,-38,-38,-38,-38,-38
S,/,return f/g,0,2,0,0,0,0,0,0,0,ModFrmHilDElt,,0,0,ModFrmHilDElt,,ModFrmHilDElt,-38,-38,-38,-38,-38
S,Inverse,return 1/f,0,1,0,0,0,0,0,0,0,ModFrmHilDElt,,ModFrmHilDElt,-38,-38,-38,-38,-38
S,^,return f^e,0,2,0,0,0,0,0,0,0,148,,0,0,ModFrmHilDElt,,ModFrmHilDElt,-38,-38,-38,-38,-38
S,!,returns f such that a_I := R!a_I,0,2,0,0,0,0,0,0,0,ModFrmHilDElt,,0,0,-44,,ModFrmHilDElt,-38,-38,-38,-38,-38
S,IsCoercible,,0,2,0,0,0,0,0,0,0,-1,,0,0,ModFrmHilD,,36,-1,-38,-38,-38,-38
S,Swap,"given a hilbert modular form f(z_1, z_2), returns the swapped form f(z_2,z_1)",0,1,0,0,0,0,0,0,0,ModFrmHilDElt,,ModFrmHilDElt,-38,-38,-38,-38,-38
S,GaloisInvariantBasis,returns a basis for the GaLois invariant subspace,1,2,1,82,0,148,3,0,0,0,0,0,0,0,82,,0,0,217,,0,0,ModFrmHilD,,82,-38,-38,-38,-38,-38
S,LSeriesHil,The L-series of the Hilbert modular newform f,0,1,0,0,0,0,0,0,0,ModFrmHilDElt,,LSer,-38,-38,-38,-38,-38
