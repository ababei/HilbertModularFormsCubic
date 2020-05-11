/*****
ModFrmHilDElt
*****/

////////// ModFrmHilDElt attributes //////////

declare type ModFrmHilDElt;
declare attributes ModFrmHilDElt:
  Parent, // M
  Weight, // SeqEnum[RngIntElt] : a sequence of [ModFrmHilDBaseField : QQ] integers
  Level, // RngOrdIdl : ideal of Integers(Field)
  Coefficients; // SeqEnum : this also keeps track of coefficient ring

////////// ModFrmHilDElt fundamental intrinsics //////////

intrinsic PercentM(f::ModFrmHilDElt) -> MonStgElt
  {returns a string to produce f in a magma session.}
  return Sprintf("HMF(%m, %m, %m)", Parent(f), Weight(f), Coefficients(f));
end intrinsic;

intrinsic Print(f::ModFrmHilDElt, level::MonStgElt : num_coeffs := 20)
  {}
  if level in ["Default", "Minimal", "Maximal"] then
    M := Parent(f);
    ideals := Ideals(M);
    k := Weight(f);
    prec := Precision(f);
    coeffs := Coefficients(f);
    N := Level(f);
    printf "Hilbert modular form expansion with precision %o.\n", prec;
    printf "Level: (Norm, Ideal) = (%o, %o)\n", Norm(N),  Generators(N);
    printf "Weight: %o\n", k;
    printf "Parent: %o\n", M;
    printf "Coefficients:\n";
    printf "\n\t(Norm, nn)  |--->   a_nn";
    printf "\n\t(%o, %o)  |--->   %o", 0,  Generators(ideals[1]), coeffs[1];
    for i:= 2 to Min(num_coeffs, #coeffs) do
      printf "\n\t(%o, %o)  |--->   %o", Norm(ideals[i]),  Generators(ideals[i]), coeffs[i];
    end for;
  elif level eq "Magma" then
    printf "%o", PercentM(f);
  else
    error "not a valid printing level.";
  end if;
end intrinsic;

intrinsic 'in'(x::., y::ModFrmHilDElt) -> BoolElt
  {}
  return false;
end intrinsic;

intrinsic IsCoercible(x::ModFrmHilDElt, y::.) -> BoolElt, .
  {}
  return false;
end intrinsic;

intrinsic 'eq'(f::ModFrmHilDElt, g::ModFrmHilDElt) -> BoolElt
  {compares Parent, Weight, and Coefficients.}
  if Parent(f) ne Parent(g) then
    return false;
  end if;
  if Weight(f) ne Weight(g) then
    return false;
  end if;
  if Coefficients(f) ne Coefficients(g) then
    return false;
  end if;
  return true;
end intrinsic;

////////// ModFrmHilDElt access to attributes //////////

intrinsic Parent(f::ModFrmHilDElt) -> ModFrmHilD
  {returns ModFrmHilD space where f lives.}
  return f`Parent;
end intrinsic;

intrinsic Weight(f::ModFrmHilDElt) -> SeqEnum[RngIntElt]
  {returns weight of f.}
  return f`Weight;
end intrinsic;

intrinsic Coefficients(f::ModFrmHilDElt) -> SeqEnum
  {returns coefficients of f as SeqEnum.}
  return f`Coefficients;
end intrinsic;

intrinsic CoefficientsParent(f::ModFrmHilDElt) -> SeqEnum
  {returns coefficients of f as SeqEnum.}
  return Parent(f`Coefficients[1]);
end intrinsic;

intrinsic BaseField(f::ModFrmHilDElt) -> FldNum
  {returns base field of parent of f.}
  return BaseField(Parent(f));
end intrinsic;

intrinsic Level(f::ModFrmHilDElt) -> RngOrdIdl
  {returns level of parent of f.}
  return f`Level;
end intrinsic;

intrinsic Precision(f::ModFrmHilDElt) -> RngIntElt
  {returns precision of parent of f.}
  return Precision(Parent(f));
end intrinsic;

intrinsic Ideals(f::ModFrmHilDElt) -> SeqEnum[RngOrdIdl]
  {returns ideals indexing f up to bound on the norm.}
  return Ideals(Parent(f));
end intrinsic;

intrinsic Dictionary(f::ModFrmHilDElt) -> Assoc
  {returns dictionary of (parent of) f.}
  return Dictionary(Parent(f));
end intrinsic;



////////// ModFrmHilDElt elementary creation functions //////////

intrinsic ModFrmHilDEltInitialize() -> ModFrmHilDElt
  {Create an empty ModFrmHilDElt object.}
  M := New(ModFrmHilDElt);
  return M;
end intrinsic;

intrinsic ModFrmHilDEltCopy(f::ModFrmHilDElt) -> ModFrmHilDElt
  {new instance of ModFrmHilDElt.}
  g := ModFrmHilDEltInitialize();
  for attr in GetAttributes(Type(f)) do
    if assigned f``attr then
      g``attr := f``attr;
    end if;
  end for;
  return g;
end intrinsic;


intrinsic HMF(M::ModFrmHilD, N::RngOrdIdl, k::SeqEnum[RngIntElt], coeffs::SeqEnum) -> ModFrmHilDElt
  {Given M and a SeqEnum of coefficients, return ModFrmHilDElt with parent M.
  WARNING: user is responsible for coefficients and order...proceed with caution.
  }
  // assertions
  ideals := Ideals(M);
  if #coeffs ne #ideals then
    error "not the right amount of coefficients to match precision of M.";
  end if;
  assert #k eq Degree(BaseField(M));
  // make form f
  f := ModFrmHilDEltInitialize();
  f`Parent := M;
  f`Level := N;
  f`Weight := k;
  f`Coefficients := coeffs;
  return f;
end intrinsic;

intrinsic HMF(M::ModFrmHilD, N::RngOrdIdl, k::SeqEnum[RngIntElt], cmap::Assoc) -> ModFrmHilDElt
  {given a ModFrmHilD, a weight k, and an associative array of coefficients on the ideals, return ModFrmHilDElt.}
  ideals := Ideals(M);
  if #Keys(cmap) ne #ideals then
    error "not the right amount of coefficients to match precision of M.";
  end if;
  dict_ideal := Dictionary(M);
  coeffs_seqenum := [ cmap[ideals[i]]  : i in [1..#Keys(cmap)] ];
  return HMF(M, N, k, coeffs_seqenum);
end intrinsic;

intrinsic HMFZero(M::ModFrmHilD, N::RngOrdIdl, k::SeqEnum[RngIntElt]) -> ModFrmHilDElt
  {create zero ModHilFrmDElt of weight k.}
  coeffs := [ 0 : i in [1..#Ideals(M)]];
  return HMF(M, N, k, coeffs);
end intrinsic;




////////// ModFrmHilDElt user convenience functions //////////


intrinsic GetCoefficient(f::ModFrmHilDElt, I::RngOrdIdl) -> RngElt
  {returns a_I}
  return Coefficients(f)[Dictionary(f)[I]];
end intrinsic;

intrinsic SetCoefficient(f::ModFrmHilDElt, I::RngOrdIdl, c::RngElt)
  {sets a_I to c}
  f`Coefficients[Dictionary(f)[I]] := c;
end intrinsic;

intrinsic GetCoefficient(f::ModFrmHilDElt, nu::RngOrdElt) -> RngElt
  {returns a_nu}
  if not assigned Parent(f)`Representatives then
    error "You must first equip the space with a multiplication table";
  end if;
  return Coefficients(f)[DictionaryRepresentatives(f)[nu]];
end intrinsic;

intrinsic SetCoefficient(f::ModFrmHilDElt, nu::RngOrdElt, c::RngElt)
  {sets a_nu to c}
  if not assigned Parent(f)`Representatives then
    error "You must first equip the space with a multiplication table";
  end if;
  f`Coefficients[DictionaryRepresentatives(f)[nu]] := c;
end intrinsic;

intrinsic GetCoefficientsIdeal(f::ModFrmHilDElt) -> Assoc
  {given a ModFrmHilDElt, return associative array of coefficients on the ideals.}
  coeffs := Coefficients(f);
  ideals := Ideals(f);
  result := AssociativeArray();
  for i := 1 to #ideals do
    result[ideals[i]] := coeffs[i];
  end for;
  return result;
end intrinsic;

intrinsic GetCoefficientsRepresentatives(f::ModFrmHilDElt) -> Assoc
  {given a ModFrmHilDElt, return associative array of coefficients on the representatives.}
  if not assigned Parent(f)`Representatives then
    error "You must first equip the space with a multiplication table";
  end if;
  coeffs := Coefficients(f);
  reps := Representatives(f);
  result := AssociativeArray();
  for i := 1 to #reps do
    result[reps[i]] := coeffs[i];
  end for;
  return result;
end intrinsic;

////////// ModFrmHilDElt creation functions //////////

//TODO: change hecke_eigenvalues to a list
intrinsic EigenformToHMF(M::ModFrmHilD, N::RngOrdIdl, k::SeqEnum[RngIntElt], hecke_eigenvalues::Assoc) -> ModFrmHilDElt
  {Construct the ModFrmHilDElt in M determined (on prime ideals up to norm prec) by hecke_eigenvalues.}
  // pull info from M
  F := BaseField(M);
  prec := Precision(M);
  // a prime
  pp := Random(Keys(hecke_eigenvalues));
  // assertions
  if not Set(PrimesUpTo(prec, F)) subset Keys(hecke_eigenvalues) then
    print "Not enough primes";
    assert false;
  end if;

  k0 := Max(k);

  // power series ring
  log_prec := Floor(Log(prec)/Log(2)); // prec < 2^(log_prec+1)
  ZK := Parent(hecke_eigenvalues[pp]);
  ZKX<X, Y> := PolynomialRing(ZK, 2);
  R<T> := PowerSeriesRing(ZKX : Precision := log_prec + 1);
  // If good, then 1/(1 - a_p T + Norm(p) T^2) = 1 + a_p T + a_{p^2} T^2 + ...
  // If bad, then 1/(1 - a_p T) = 1 + a_p T + a_{p^2} T^2 + ...
  recursion := Coefficients(1/(1 - X*T + Y*T^2));
  ideals := Ideals(M);
  coeffs := [ZK!0: i in ideals];
  set := [false : c in coeffs];
  coeffs[1] := 0; //a_0
  coeffs[2] := 1; //a_1
  set[1] := true;
  set[2] := true;
  dict := Dictionary(M);
  for i := 1 to #coeffs do
    if not set[i] then
      // is[i] is a prime
      pp := ideals[i];
      assert IsPrime(pp);
      coeffs[i] := hecke_eigenvalues[pp];
      set[i] := true;
      Np := Norm(pp)^(k0-1);
      // if pp is bad
      if N subset pp then
        Np := 0;
      end if;

      r := 2;
      pp_power := pp * pp;
      //deals with powers of p
      while pp_power in Keys(dict) do
        ipower := dict[pp_power];
        coeffs[ipower] := Evaluate(recursion[r + 1], [coeffs[i], Np]);
        set[ipower] := true;
        pp_power *:= pp;
        r +:= 1;
      end while;

      //deal with multiples of its powers by smaller numbers
      for j := 3 to #coeffs do
        if set[j] and not (ideals[j] subset pp) then
          mm := ideals[j];
          pp_power := pp;
          mmpp_power := mm * pp_power;
          while mmpp_power in Keys(dict) do
            l := dict[mmpp_power];
            assert set[l] eq false;
            ipower := dict[pp_power];
            // a_{m * pp_power} := a_{m} * a_{pp_power}
            coeffs[l] := coeffs[j] * coeffs[ipower];
            set[l] := true;
            mmpp_power *:= pp;
            pp_power *:= pp;
          end while;
        end if; //check if it's set
      end for; //loop in j

    end if; // check if it's set
  end for; // loop in i
  for i := 1 to #coeffs do
    assert set[i];
  end for;
  return HMF(M, N, k, coeffs);
end intrinsic;

intrinsic NewformsToHMF(M::ModFrmHilD, N::RngOrdIdl, k::SeqEnum[RngIntElt]) -> SeqEnum[ModFrmHilDElt]
  {returns Hilbert newforms}
  F := BaseField(M);
  prec := Precision(M);
  MF := HilbertCuspForms(F, N, k);
  S := NewSubspace(MF);
  newspaces  := NewformDecomposition(S);
  newforms := [* Eigenform(U) : U in newspaces *];
  eigenvalues := AssociativeArray();
  primes := PrimesUpTo(prec, F);
  HMFnewforms := [];
  for newform in newforms do
    for pp in primes do
        eigenvalues[pp] := HeckeEigenvalue(newform, pp);
    end for;
    ef := EigenformToHMF(M, N, k, eigenvalues);
    Append(~HMFnewforms, ef);
  end for;
  return HMFnewforms;
end intrinsic;


/*
intrinsic NewformsToHMF2(M::ModFrmHilD, k::SeqEnum[RngIntElt]) -> SeqEnum[ModFrmHilDElt]
  {returns Hilbert newforms}
  F := BaseField(M);
  N := Level(M); //input
  prec := Precision(M);
  HeckeEigenvalue := HeckeEigenvalues(M);
  key :=  [* N, k *];
  if not IsDefined(M, key) then
    MF := HilbertCuspForms(F, N, k);
    S := NewSubspace(MF);
    newspaces  := NewformDecomposition(S);
    newforms := [* Eigenform(U) : U in newspaces *];
    primes := Primes(M);
    EVnewforms := [];
    for newform in newforms do
      eigenvalues := [];
      for i in [1..#primes] do
          eigenvalues[i] := HeckeEigenvalue(newform, primes[i]);
      end for;
      Append(~EVnewforms, eigenvalues);
    end for;
    HeckeEigenvalue[key] := EVnewforms;
  else
    EVnewforms := HeckeEigenvalue[key];
  end if;
  HMFnewforms := [];
  for eigenvalues in EVnewforms do
      ef := EigenformToHMF(M, k, eigenvalues); //FIXME, this is not correct
      Append(~HMFnewforms, ef);
    end for;
  return HMFnewforms;
end intrinsic;
*/

intrinsic GaloisOrbit(f::ModFrmHilDElt) -> SeqEnum[ModFrmHilDElt]
  {returns the full Galois orbit of a modular form}
  M := Parent(f);
  k := Weight(f);
  coeff := Coefficients(f);
  K := NumberField(CoefficientsParent(f));
  G, Pmap, Gmap := AutomorphismGroup(K);
  result := [];
  for g in G do
    Append(~result, HMF(M, Level(f), k, [Gmap(g)(elt) : elt in coeff]) );
  end for;
  return result;
end intrinsic;

intrinsic GaloisOrbitDescent(f::ModFrmHilDElt) -> SeqEnum[ModFrmHilDElt]
  {returns the full Galois orbit of a modular form over Q}
  M := Parent(f);
  k := Weight(f);
  coeff := Coefficients(f);
  K := NumberField(CoefficientsParent(f));
  result := [];
  for b in Basis(K) do
    Append(~result, HMF(M, Level(f), k, [Trace(b * elt) : elt in coeff]) );
  end for;
  return result;
end intrinsic;

intrinsic CuspFormBasis(M::ModFrmHilD, N::RngOrdIdl, k::SeqEnum[RngIntElt]) -> SeqEnum[ModFrmHilDElt], RngIntElt
  {returns a basis for M of weight k}
  prec := Precision(M);
  F := BaseField(M);
  ideals := Ideals(M);
  dict := Dictionary(M);
  set_ideals := Keys(dict);
  zero_coeffs := [0 : i in [1..#ideals]];
  basis := [];
  newforms_dimension := 0;
  //TODO should we sieve?
  for dd in Divisors(N) do
    basis := [];
    orbit_representatives := NewformsToHMF(M, N, k);
    orbits := [GaloisOrbitDescent(elt) : elt in orbit_representatives];
    old_space_basis := &cat orbits;
    old_space_coeffs := [Coefficients(elt) : elt in old_space_basis];
    for ee in Divisors(N/dd) do
      new_coeffs_ee := [ zero_coeffs : i in [1..#old_space_basis]];
      for i := 1 to #ideals do
        Iee := ideals[i]*ee;
        if Iee in set_ideals then
          Iee_cord := dict[Iee];
          for j := 1 to #old_space_basis do
            new_coeffs_ee[j][Iee_cord] := old_space_coeffs[j][i];
          end for;
        else
          //this assumes ideals are sorted by norm!!
          break i;
        end if;
      end for;
      basis := basis cat [HMF(M, N, k, elt) : elt in new_coeffs_ee];
    end for;
    if dd eq N then
      if #orbits eq 0 then
        newforms_dimension := 0;
      else
        newforms_dimension := &+[ #orbit : orbit in orbits];
      end if;
    end if;
  end for;
  return basis, newforms_dimension;
end intrinsic;

intrinsic HeckeOperator(f::ModFrmHilDElt, nn::RngOrdIdl : Basis:=[]) -> ModFrmHilDElt
  {returns T(n)(f) for the trivial character}
  return HeckeOperator(f, nn, HeckeCharacterGroup(Level(f))! 1 : Basis:=Basis);
end intrinsic;


//TODO:
//Tests:
// - Apply Hecke on a Galois Orbit, and see that it doesn't move
// - Apply Hecke to a Eisensten series, and check that is a multiple
// - Apply Hecke to a Theta series, and see if we get the whole space




intrinsic HeckeOperator(f::ModFrmHilDElt, nn::RngOrdIdl, chi::GrpHeckeElt : Basis:= []) -> ModFrmHilDElt
  {returns T(n)(f) with loss of precision. If a basis for the space is provided, returns T(n)(f) without loss of precision}
  M := Parent(f);
  fcoeffs := Coefficients(f);
  ideals := Ideals(f);
  dict := Dictionary(f);
  ZC := Parent(fcoeffs[1]);
  k0 := Max(Weight(f));
  prec:=Precision(M);
  //We work in smaller precision and obtain a function the space of precision Prec/Norm(nn)
  idealssmallprec:=[ideals[1]];
  for x:=2 to #ideals do
    if Norm(ideals[x]) le prec/Norm(nn) then
      idealssmallprec:= idealssmallprec cat [ideals[x]];
    end if;
  end for;
  coeffssmallprec:=AssociativeArray();
  coeffssmallpreclist := [ZC!0 : i in [1..#idealssmallprec]];
  for i:=1 to #idealssmallprec do
    c := 0;
    // loop over divisors
    // Formula 2.23 in Shimura - The Special Values of the zeta functions associated with Hilbert Modular Forms
    for aa in Divisors(ideals[i] + nn) do
      c +:= chi(aa) * Norm(aa)^(k0 - 1) * fcoeffs[ dict[ aa^(-2) * (ideals[i] * nn)]];
      end for;
    coeffssmallprec[ideals[i]] := c;
    coeffssmallpreclist[i]:=c;
    end for;
  if #Basis ne 0 then //if provided a basis, we reconstruct the form to the same precision
    smallprecbasis:=[];
    for g in Basis do
      Append(~smallprecbasis, [Coefficients(g)[i] : i in [1..#idealssmallprec]]);
    end for;
    LinComb:=LinearDependence([coeffssmallpreclist] cat smallprecbasis);
    if ISA(Type(LinComb), RngIntElt) or #LinComb ne 1 then
      return "Error: either precision is too small, or the sequence provided is not a basis";
    end if;
    g:=0*f;
    for i:=2 to #Basis+1 do
      g:=g-LinComb[1][i]*Basis[i-1];
    end for;
    return g;
  end if;
  if #Basis eq 0 then
    M1:=HMFSpace(BaseField(M), Floor(prec/Norm(nn)));
    print "Warning: the Hecke operator calculation decreases precision to ",
    Floor(prec/Norm(nn));
    return HMF(M1, Level(f), Weight(f), coeffssmallprec);
  end if;
end intrinsic;


// TODO needs testing
intrinsic EisensteinSeries(M::ModFrmHilD, N::RngOrdIdl, eta::GrpHeckeElt, psi::GrpHeckeElt, k::SeqEnum[RngIntElt]) -> ModFrmHilDElt
  {Let aa*bb be the modulus of psi*eta^-1. Return the Eisenstein series E_k(eta,psi) in M_k(aa*bb,eta*psi).}
  Cl := NarrowClassGroup(M);
  mp := NarrowClassGroupMap(M);
  X := Parent(eta);
  assert X eq Parent(psi);
  K := X`TargetRing; // where the character values live
  if #Cl eq 1 then
    tt := mp(Cl.0);
  else
    error "not implemented for narrow class number > 1.";
  end if;
  n := Degree(BaseField(M));
  assert #SequenceToSet(k) eq 1; // only parallel weight for now
  nn := N;
  // aa := Conductor(eta);
  aa := Modulus(eta);
  // bb := Conductor(psi);
  bb := Modulus(psi);
  assert nn subset aa;
  assert nn subset bb;
  Haa := HeckeCharacterGroup(aa);
  Hbb := HeckeCharacterGroup(bb);
  ideals := Ideals(M);
  coeffs := [ K!0 : i in [1..#ideals]];
  if k[1] ge 2 then
    // constant term
    if aa eq ideal<Order(aa)|1> then
      prim := AssociatedPrimitiveCharacter(psi*eta^(-1));
      coeffs[1] := 2^(-n)*(eta^(-1))(tt)*LValue_Recognized(M, N, prim, k);
    else
      coeffs[1] := 0;
    end if;
  elif k[1] eq 1 then // wt 1 case
    if aa eq ideal<Order(aa)|1> and bb ne ideal<Order(bb)|1> then
      prim := AssociatedPrimitiveCharacter(psi*eta^(-1));
      coeffs[1] := 2^(-n)*(eta^(-1))(tt)*LValue_Recognized(M, N, prim, k);
    elif aa ne ideal<Order(aa)|1> and bb eq ideal<Order(bb)|1> then
      prim := AssociatedPrimitiveCharacter(psi^(-1)*eta);
      coeffs[1] := 2^(-n)*(psi^(-1))(tt)*LValue_Recognized(M, N, prim, k);
    elif aa eq ideal<Order(aa)|1> and bb eq ideal<Order(bb)|1> then
      prim1 := AssociatedPrimitiveCharacter(psi*eta^(-1));
      prim2 := AssociatedPrimitiveCharacter(psi^(-1)*eta);
      coeffs[1] := 2^(-n)*((eta^(-1))(tt)*LValue_Recognized(M, N, prim1, k)
                          +(psi^(-1))(tt)*LValue_Recognized(M, N, prim2, k));
    elif aa ne ideal<Order(aa)|1> and bb ne ideal<Order(bb)|1> then
      coeffs[1] := 0;
    end if;
  else // nonpositive and half-integral weights...
    error "Not implemented";
  end if;
  // other terms
  for i := 2 to #ideals do // 2 here assumes #Cl = 1 FIXME
    mm := ideals[i];
    sum := 0;
    for rr in Divisors(mm) do
      sum +:= eta(mm/rr)*psi(rr)*Norm(rr^(k[1]-1));
    end for;
    coeffs[i] := sum;
  end for;
  if not (coeffs[1] in [0,1]) then
    factor := 1/coeffs[1];
    coeffs := [factor * elt : elt in coeffs];
  end if;
  if IsIsomorphic(K, RationalsAsNumberField()) then
    coeffs := [ Rationals() ! elt : elt in coeffs ];
  end if;
  return HMF(M, N, k, coeffs);
end intrinsic;

// TODO finish this and use in EisensteinSeries intrinsic

//Toolbox function to use in the Eisenstein series function--gives us an L value
intrinsic LValue_Recognized(M::ModFrmHilD, N::RngOrdIdl, prim::GrpHeckeElt, k::SeqEnum[RngIntElt]) -> FldNumElt
{This is a toolbox function to compute L values in the right space}
// Lf := LSeries(prim : Precision := 50);
// TODO clean up precision
// Maybe a separate function to compute L-values?
K := Parent(prim)`TargetRing; // where the character values live
Lf := LSeries(prim : Precision := 100);
LSetPrecision(Lf, 100); // do we need this?
Lvalue := Evaluate(Lf, 1-k[1]);
// figure out the right place
primes := PrimesUpTo(Precision(M), BaseField(M));
places := InfinitePlaces(K);
i := 1;
while #places gt 1 and i le #primes do
  pp := primes[i];
  app := prim(pp);
  places := [pl : pl in places | Evaluate(app, pl) eq -Coefficients(EulerFactor(Lf, pp : Degree := 1))[2] ];
  // why is this the right way to find the correct place to recognize?
  i +:=1;
end while;
assert #places eq 1;
pl := places[1];
CC<I> := ComplexField(Precision(Lvalue));
return RecognizeOverK(CC!Lvalue, K, pl, false);
end intrinsic;

intrinsic Basis(M::ModFrmHilD, N::RngOrdIdl, k::SeqEnum[RngIntElt]) -> SeqEnum[ModFrmHilDElt], RngIntElt
  { returns a Basis for the space }
  CB, newforms_dimension := CuspFormBasis(M, N, k);
  H := HeckeCharacterGroup(N);
  //FIXME this is wrong for level not 1!
  //print "FIXME this is wrong for level not 1!";
  eta := H ! 1;
  /* psi := H ! 1; */
  psi := eta;
  E := EisensteinSeries(M, N, eta, psi, k);
  return [E] cat CB, newforms_dimension;
end intrinsic;

////////// ModFrmHilDElt arithmetic //////////

intrinsic '*'(c::RngIntElt, f::ModFrmHilDElt) -> ModFrmHilDElt
  {scale f by integer c.}
  coeffs := Coefficients(f);
  ZK := Parent(coeffs[1]);
  assert c in ZK;
  czk := ZK ! c;
  new_coeffs := [ czk * elt : elt in coeffs];
  return HMF(Parent(f), Level(f), Weight(f), new_coeffs);
end intrinsic;

intrinsic '*'(c::Any, f::ModFrmHilDElt) -> ModFrmHilDElt
  {scale f by some scalar c.}
  coeffs := Coefficients(f);
  ZK := Parent(coeffs[1]);
  //if ZK eq Rationals() then
  //  K := Parent(c);
  //  else
  //    K := Compositum(Parent(c), FieldOfFractions(ZK));
  //  end if;
  //ck := K!c;
  new_coeffs := [ c * elt : elt in coeffs];
  return HMF(Parent(f), Level(f), Weight(f), new_coeffs);
end intrinsic;

intrinsic '+'(f::ModFrmHilDElt, g::ModFrmHilDElt) -> ModFrmHilDElt
  {return f+g.}
  // assert Parent(f) eq Parent(g);
  assert Level(f) eq Level(g);
  assert Weight(f) eq Weight(g);
  assert BaseField(Parent(f)) eq BaseField(Parent(g));
  prec_f := Precision(Parent(f));
  prec_g := Precision(Parent(g));
  if prec_f lt prec_g then
    new_coeffs := [ Coefficients(f)[i] + Coefficients(g)[i] : i in [1..#Coefficients(f)] ];
    return HMF(Parent(f), Level(f) meet Level(g), Weight(f), new_coeffs);
  else
    new_coeffs := [ Coefficients(f)[i] + Coefficients(g)[i] : i in [1..#Coefficients(g)] ];
    return HMF(Parent(g), Level(g) meet Level(f), Weight(g), new_coeffs);
  end if;
end intrinsic;

intrinsic '-'(f::ModFrmHilDElt, g::ModFrmHilDElt) -> ModFrmHilDElt
  {return f-g.}
  // assert Parent(f) eq Parent(g);
  assert Level(f) eq Level(g);
  assert Weight(f) eq Weight(g);
  assert BaseField(Parent(f)) eq BaseField(Parent(g));
  prec_f := Precision(Parent(f));
  prec_g := Precision(Parent(g));
  if prec_f lt prec_g then
    new_coeffs := [ Coefficients(f)[i] - Coefficients(g)[i] : i in [1..#Coefficients(f)] ];
    return HMF(Parent(f), Level(f) meet Level(g), Weight(f), new_coeffs);
  else
    new_coeffs := [ Coefficients(f)[i] - Coefficients(g)[i] : i in [1..#Coefficients(g)] ];
    return HMF(Parent(g), Level(g) meet Level(f), Weight(g), new_coeffs);
  end if;
end intrinsic;

// TODO for varied precision
intrinsic '*'(f::ModFrmHilDElt, g::ModFrmHilDElt) -> ModFrmHilDElt
  {return f*g}
  N := Level(f) meet Level(g);
  M := Parent(f);
  //assert Parent(f) eq Parent(g);
  if not assigned M`MultiplicationTable then
    assert HMFEquipWithMultiplication(M);
  end if;
  prec := Precision(M);
  fcoeffs := Coefficients(f);
  gcoeffs := Coefficients(g);
  ZC := FieldOfFractions(Parent(fcoeffs[1]));
  MTable := MultiplicationTable(M);
  assert ZC eq FieldOfFractions(Parent(gcoeffs[1]));
  coeffs := [ZC!0 :  i in [1..#fcoeffs]];
  for i := 1 to #fcoeffs do
    c := ZC!0;
    for pair in MTable[i] do
      c +:= fcoeffs[ pair[1] ] * gcoeffs[ pair[2] ];
    end for;
    coeffs[i] := c;
  end for;
  kf := Weight(f);
  kg := Weight(g);
  k := [ kf[i] + kg[i] : i in [1..#kf] ];
  return HMF(M, N, k, coeffs);
end intrinsic;

intrinsic '/'(f::ModFrmHilDElt, g::ModFrmHilDElt) -> ModFrmHilDElt
  {return f/g}
  N := Level(f) meet Level(g);
  M := Parent(f);
  assert Parent(f) eq Parent(g);
  if not assigned M`MultiplicationTable then
    assert HMFEquipWithMultiplication(M);
  end if;
  prec := Precision(M);
  fcoeffs := Coefficients(f);
  gcoeffs := Coefficients(g);
  // TODO make it work if Valuation(g) <= Valuation(f)
  require gcoeffs[1] ne 0: "Denominator must have nonzero constant term";
  ib0 := 1/gcoeffs[1];
  ZC := CoefficientsParent(f);
  MTable := MultiplicationTable(M);
  assert ZC eq CoefficientsParent(g);
  coeffs := [ZC!0 :  i in [1..#fcoeffs]];
  for i := 1 to #fcoeffs do
    c := fcoeffs[i];
    for pair in MTable[i] do
      if pair[1] ne 1 then
        c -:= gcoeffs[ pair[1] ] * coeffs[ pair[2] ];
      else
        assert pair[2] eq i;
      end if;
    end for;
    coeffs[i] := c * ib0;
  end for;
  kf := Weight(f);
  kg := Weight(g);
  k := [ kf[i] - kg[i] : i in [1..#kf] ];
  return HMF(M, N, k, coeffs);
end intrinsic;

intrinsic Inverse(f::ModFrmHilDElt) -> ModFrmHilDElt
  {return 1/f}
  return (Parent(f) ! ( CoefficientsParent(f) ! 1) ) / f;
end intrinsic;

intrinsic '^'(f::ModFrmHilDElt, n::RngIntElt) -> ModFrmHilDElt
  {return f^e}
  if n lt 0 then
    f := Inverse(f);
  end if;
  g := Parent(f) ! (CoefficientsParent(f) ! 1);
  if n eq 0 then
    return g;
  end if;
  while n gt 1 do
    if n mod 2 eq 0 then
      f := f * f;
      n := Integers() ! (n/2);
    else
      g := f * g;
      f := f * f;
      n := Integers() ! ((n - 1)/2);
    end if;
  end while;
  return f * g;
end intrinsic;



intrinsic '!'(R::Rng, f::ModFrmHilDElt) -> ModFrmHilDElt
  {returns f such that a_I := R!a_I}
  coeffs := [R!c : c in Coefficients(f)];
  return HMF(Parent(f), Level(f), Weight(f), coeffs);
end intrinsic;


intrinsic IsCoercible(M::ModFrmHilD, f::.) -> BoolElt, .
  {}
  if ISA(Type(f), RngElt) then
    P := Parent(f);
    N := 1*Integers(M); // FIXME only level 1 for now
    coeffs := [P!0 : c in [1..#Ideals(M)]];
    coeffs[1] := f;
    k := [0 : c in [1..Degree(BaseField(M))]];
    return true, HMF(M, N, k, coeffs);
  end if;

  if Type(f) ne ModFrmHilDElt then
    return false;
  end if;

  if Parent(f) eq M then
    return true, f;
  elif (BaseField(M) eq BaseField(f)) and (Precision(M) eq Precision(f)) then
    coeffs := Coefficients(f);
    return true, HMF(M, Level(f), Weight(f) , coeffs);
  else
    return false;
  end if;
end intrinsic;

/*
intrinsic '!'(M::ModFrmHilD, f::ModFrmHilDElt) -> ModFrmHilDElt
  {returns f with parent M}
  nn := Level(M);
  nnf := Level(f);
  assert nn subset nnf; // nnf divides nn
  coeffs := Coefficients(f);
  return HMF(M, Weight(f) , coeffs);
end intrinsic;
*/

/*swap map
*/

intrinsic Swap(f::ModFrmHilDElt) -> ModFrmHilDElt
  {given a hilbert modular form f(z_1, z_2), returns the swapped form f(z_2,z_1)}
  M:=Parent(f);
  g:=M!(1*f);
  F:=BaseField(M);
  ZF<w>:=Integers(F);
  ideals:=Ideals(f);
  for i in ideals do
    x:=GetCoefficient(f, Conjugate(i));
    SetCoefficient(g, i, x);
  end for;
  return g;
  end intrinsic;

  intrinsic GaloisInvariantBasis(M::ModFrmHilD, N::RngOrdIdl, k::SeqEnum[RngIntElt]) -> SeqEnum[ModFrmHilDElt]
  {returns a basis for the GaLois invariant subspace}
  B:=Basis(M,N,k);
  InvariantGenerators:=[];
  for x in B do
    Append(~InvariantGenerators, 1/2*(x+Swap(x)));
  end for;
  InvariantBasis:=[];
  for x in InvariantGenerators do
    if #LinearDependence(InvariantBasis cat [x]) eq 0 then
      Append(~InvariantBasis, x);
    end if; 
  end for;
  return InvariantBasis;
end intrinsic;


intrinsic LSeriesHil(f::ModFrmHilDElt : Precision:=8, Embedding:=1) -> LSer
{The L-series of the Hilbert modular newform f}
 Mf:=Parent(f); K:=BaseField(f); N:=Level(f); W:=Weight(f);
//to do: what if f is a bianchi form?
//require f is a newform
// TO DO: won't work for indefinite? because bad Heckes not availabe
 require &and[IsEven(w) : w in W]: "All weights must be even";
 w:=W[#W];

 E:=Parent(Coefficients(f)[3]);
 A:=AbsoluteField(E);
 r,c:=Signature(A); RF:=c eq 0 select RealField else ComplexField;
 prec:= Precision;
 R1 := PowerSeriesRing(RF(prec),1+2*Degree(K)); T := R1.1;
 R2 := PolynomialRing(RF(prec)); U := R2.1; // hack: are T,U compatible?
 ip:=InfinitePlaces(A)[Embedding]; // only first place...
 scal:=Coefficients(f)[2];

 function cfK(p,d : Precision:=prec)  fp:=Degree(p); P:=Norm(p); // prec?
  if fp gt d then return 1+O(T^(d+1)); end if; // think here?
  if P le  820 then
    cc:=GetCoefficient(f,p)/scal;
  else 
    cc:=0;
  end if;
  if Degree(A) eq 1 then ap:=cc;
  else ap:=Evaluate(A!cc,ip : Precision:=prec); end if;
  eps:=Valuation(N,p) gt 0 select 0 else 1; // need character generally
  return 1-ap*(U^fp)+eps*P^(w-1)*U^(2*fp);  end function;

 function cf(p,d : Precision:=prec) // need prec compatible
  F:=Factorization(p*Integers(K));
  B:=[* cfK(f[1],d : Precision:=Precision) : f in F *];
  if &and[Type(Parent(b)) eq Type(R2) and Parent(b) eq R2 : b in B]
   then B:=[b : b in B]; else B:=[R1!b :  b in B]; end if; return &*B;
  end function;
 name:=<"L-series of Hilbert modular form",f>;
 gamma:=&cat[[0-e,1-e] where e:=(w-W[i]) div 2 : i in [1..Degree(K)]];
 L:=LSeries(w,gamma,Norm(N)*Discriminant(Integers(K))^2,cf :
      Name:=name, Precision:=prec, Sign:=1);
 L`cffuncK:=cfK; L`degreeK:=<2,K>; L`condK:=N;
 IP:=InfinitePlaces(K);
 L`hodgeK:=
   &cat[[<[w-W[i],W[i]-1,2],IP[i]>,<[W[i]-1,w-W[i],2],IP[i]>] : i in [1..#IP]];
 return L; end intrinsic;

