
//Base functions
function Lam(f,s: Precision:=8, Embedding:=1) //the LStar function would probably work just as well
  F:=BaseField(Parent(f));
  Nd:=Abs(Discriminant(F));
  n:=Degree(F);
  if Type(f) eq ModFrmHilDElt then L:=LSeriesHil(f: Precision:=Precision, Embedding:=Embedding);
  	else L:=LSeries(f: Precision:=Precision, Emb:=Embedding);
  end if; 
  C<i>:=ComplexField(Precision);
  return Nd^s*(Gamma(s)^n)*(2*Pi(C))^(-n*s)*Evaluate(L,s);
end function;



//Period polynomial functions

function rfx(f: Precision:=8, Embedding:=1)
	F:=BaseField(Parent(f));
	Nd:=Discriminant(F);
	n:=Degree(F);
	C<i> := ComplexField(Precision);
	if Type(f) eq ModFrmHilDElt then 
		W:=Weight(f); 
		if #Set(W) ne 1 then return "Error: only supported for parallel weight"; end if;
		k:=Random(Set(W));
  	else 
  		W:=Weight(Parent(f));
		if #Set(W) ne 1 then return "Error: only supported for parallel weight"; end if;
		k:=Random(Set(W));
  	end if; 
	pol:=0;
	PolRing<x>:=PolynomialRing(Rationals());
	for m:=0 to k-2 do
  		pol+:= (-1)^m*i^(n*(k-m-1))*Binomial(k-2,m)*Lam(f,k-m-1: Precision:=Precision, Embedding:=Embedding)*(x)^m;
	end for;
	return pol;// 1/2*(Evaluate(pol, x)+Evaluate(pol, -x)), 1/2*(Evaluate(pol, x)-Evaluate(pol, -x)); 
end function;

function L2Dist(L1, L2)
	sum:=0;
	if #L1 ne #L2 then
		print("The sequences don't have the same length");
		return -1;
	else
		for x in [1..#L1] do
			sum+:= (L1[x]-L2[x])^2;
		end for;
		return Real(Sqrt(sum));
	end if;
end function;



function CompEmb(f: Bound:=1000)
  K:=FieldOfFractions(Parent(Coefficients(f)[4]));
  embs:=[**];
  for i:=1 to #Conjugates(K!1) do
  	coeffs := [Conjugates(K!Coefficients(f)[x])[i]: x in [1..Bound]];
    Append(~embs, coeffs);
  end for;
  return embs;
 end function;




function CoeffsUpToPrec(f, n)
  K:=FieldOfFractions(Parent(Coefficients(f)[4]));
  ids:=IdealsUpTo(n, BaseField(Parent(f)));
  dict:=Dictionary(f);
  coe:=AssociativeArray();
  coe[0*Integers(BaseField(Parent(f)))]:=Coefficients(f)[1];
  for x in ids do
  	coe[x]:=Coefficients(f)[dict[x]];
  end for;
  g:=HMF(HMFSpace(BaseField(Parent(f)), n), Level(f), Weight(f), coe);
  return Coefficients(g);
 end function;



 function CheckRoots(f: Prec:=8, Emb:=1)
 	C<i>:=ComplexField();
 	F:=BaseField(Parent(f));
 	n:=Degree(F);
 	D:=Abs(Discriminant(F));
 	print(D);
 	if Type(f) eq ModFrmHilDElt then 
 		L:=LSeriesHil(f: Precision:=Prec, Embedding:=Emb);
		W:=Weight(f); 
		if #Set(W) ne 1 then return "Error: only supported for parallel weight"; end if;
		k:=Random(Set(W));
  	else 
  		W:=Weight(Parent(f));
  		L:=LSeries(f: Precision:=Prec, Emb:=Emb);
		if #Set(W) ne 1 then return "Error: only supported for parallel weight"; end if;
		k:=Random(Set(W));
  	end if; 
  	m:=(k-2)/2;
  	//Lis:=[[Evaluate(L,m+1), Evaluate(L, 2*m+1), Evaluate(L,m+1)/Evaluate(L, 2*m+1)]];
 	s1:=1/2*Gamma(m+1)^(n-2)/Gamma(2*m+1)^(n-1)*(2*Pi(C))^(m*n)/(D^m)*Evaluate(L,m+1)/Evaluate(L, 2*m+1);
 	for j in [1..m-1] do
 		//Append(~Lis, [Evaluate(L, 2*m+1-j), Evaluate(L, 2*m+1), Evaluate(L, 2*m+1-j)/Evaluate(L, 2*m+1)]);
 		s1+:=1/Factorial(j)*(2^n*Pi(C)^n/D)^j*(Gamma(2*m+1-j)/Gamma(2*m+1))^(n-1)*Evaluate(L, 2*m+1-j)/Evaluate(L, 2*m+1);
 	end for;
 	return Abs(s1);//Lis
 end function;




// function CoeffsUpToPrec(f, n)
//   F:=BaseField(Parent(f));
//   K:=FieldOfFractions(Parent(Coefficients(f)[1*ZF][1*ZF]));
//   ids:=IdealsUpTo(n, BaseField(Parent(f)));
//   coe:=[];
//   Append(~coe, Coefficients(f)[1*ZF][0*ZF]);
//   for x in ids do
//   	Append(~coe, Coefficients(f)[1*ZF][x]);
//   end for;
//   return coe;
//  end function;

//quadratic case
QuadDiscs:=[5,8,12,13,17,21,24,29,33];
Weights:=[6,5,4,5,4,4,4,4,4];
Precs:=[10,10,10,10,10,10,10,10,10];
//

//loop through spaces; add period polynomials to list and check that roots of the derivative inside the unit disk
PerPolys:=[**]; //Returns list of all period polynomials
RemainingPerPolys:=[**];   //Returns list of period polynomials that still remain to be checked 


for j in [1..9] do
	printf "Computing period polynomials for field with discriminant  %o\n", QuadDiscs[j];
	Append(~PerPolys, [**]);
	Append(~RemainingPerPolys, [**]);
	F:=QuadraticField(QuadDiscs[j]);
	ZF:=Integers(F);
	prec:=Precs[j];
	MaxWeight:=Weights[j];
	if j  le 9 then                        ///Can modify if the calculation runs out of memory or it crashes (Windows bug)
		for l in [1..MaxWeight] do
			print(2*l);
			Append(~PerPolys[j], [**]);
			Append(~RemainingPerPolys[j], [**]);
			if l eq 5 or l eq 6 then
				H := HilbertCuspForms(F,1*ZF,[2*l, 2*l]);
				S:= NewformDecomposition(NewSubspace(H));
				for subspace in S do
					g:=Eigenform(subspace);
					K:=BaseField(g);
					embs:=Conjugates(K!1);
					for x in [1..#embs] do
						pol:=rfx(g: Precision:=15, Embedding:=x);
						Append(~PerPolys[j][l], pol);
						s:=CheckRoots(g: Prec:=4, Emb:=x);
						if Abs(s) gt 1 then
							Append(~RemainingPerPolys[j][l], pol);							
						end if;
					end for;
				end for;
			end if;
		end for;
	end if;
end for;




//Cubics

PolRing<x>:= PolynomialRing(Integers());

//First cubic field
//D=49;
t:=Cputime();
F := NumberField(x^3 - x^2 - 2*x + 1);
ZF:=Integers(F);
sizes:=[];

ids:=IdealsUpTo(150, F);
prec := 5000;                       //Can change precision as you will
M := HMFSpace(F, prec);
H := HeckeCharacterGroup(1*ZF);
Y2:=EisensteinSeries(M, 1*ZF,  H!1, H!1, [2,2,2]);   // 1 Form of weight 2
Y4:=EisensteinSeries(M, 1*ZF,  H!1, H!1, [4,4,4]);   // 1 form of weight 4
Y6:=EisensteinSeries(M, 1*ZF,  H!1, H!1, [6,6,6]);   // 1 form of weight 6
Y8:=EisensteinSeries(M, 1*ZF,  H!1, H!1, [8,8,8]);   // 1 form of weight 8




/////////Create all the possible weighted monomials from the ones created above

MonBas:=[**];
P<x2, x4, x6, x8> := PolynomialRing(RationalField(), [2,4,6,8]);
for i in [1..4] do
	Append(~MonBas, []);
	Mons:=MonomialsOfWeightedDegree(P, 2*i);
	for j in Mons do
		print(j);
		factor:=Factorization(j);
		mult:=1;
		for l in factor do
			if l[1] eq x2 then
				mult:=mult*Y2^l[2];
			end if;
			if l[1] eq x4 then
				mult:=mult*Y4^l[2];
			end if;
			if l[1] eq x6 then
				mult:=mult*Y6^l[2];
			end if;
			if l[1] eq x8 then
				mult:=mult*Y8^l[2];
			end if;
		end for;
		Append(~MonBas[i], mult);
	end for;
end for;



//////// The monomials obtained above give possible non-cusp forms. Refine it to a cusp form basis.

CuspBas:=[**];
for i:=1 to 4 do
	print(2*i);
	Append(~CuspBas, []);
	for j:=1 to #MonBas[i] do
		print([2*i,j]);
		if i ge 2 then
			for j:=1 to #MonBas[i]-1 do
				//print([2*i,j]);
				for l:=j+1 to #MonBas[i] do
					if #CuspBas[i] eq 0 then
						Append(~CuspBas[i], MonBas[i][j]-MonBas[i][l]);
						print([2*i, j,l,Coefficients(MonBas[i][j]-MonBas[i][l])[4],1]);
					else
						//print(LinearDependence([CoeffsUpToPrec(x, 100): x in CuspBas[i]] cat [CoeffsUpToPrec(MonBas[i][j]-MonBas[i][l], 100)]));
						if  #LinearDependence([CoeffsUpToPrec(x, 50): x in CuspBas[i]] cat [CoeffsUpToPrec(MonBas[i][j]-MonBas[i][l], 50)]) eq 0 then
							Append(~CuspBas[i], MonBas[i][j]-MonBas[i][l]);
							print([2*i, j,l,Coefficients(MonBas[i][j]-MonBas[i][l])[4], 2]);
						end if;
					end if;
				end for;
			end for;
		end if;
	end for;
end for;

// Run Hecke Operators to obtain basis of eigenforms for each weight. For each of them, compute period polynomials


NewBas:=[**];
BasTrack:=[**];
PerPolys:=[**];
for i in [1..#CuspBas] do
	if #CuspBas[i] ne 0 then
		Mk:=HMFSpace(F, 2000);
		Append(~NewBas,[**]);
		Append(~BasTrack,[**]);
		Append(~PerPolys,[**]);
		tempor:=[];
		CoeffLi:=[];
		HeckeOpAct:=[CoeffsUpToPrec(HeckeOperator(p, ids[2]), 50): p in CuspBas[i]];
		LinDeps:=[];
		for a in CuspBas[i] do
			Lin:=LinearDependence(HeckeOpAct cat [CoeffsUpToPrec(a, 50)]);
			Append(~LinDeps, #Lin);
			Append(~CoeffLi, [Lin[1][r]/Lin[1][#Lin[1]] : r in [1..#Lin[1]-1]]);
		end for;
		//print(LinDeps);
		H:=Matrix(CoeffLi);
		//print(H);
		//Determinant(H);
		newroots:=[r[1]: r in Factorization(CharacteristicPolynomial(H))];
		K:=ext< Rationals() | newroots>;
		HeckeMatrix:=Matrix(K, CoeffLi);
		print([2*i, #LinDeps]);
		Eig:=[t[1] : t in Eigenvalues(HeckeMatrix)];
		print "Eigenvalues computed";
		print(#Eig);
		for e in Eig do
			B:= Basis(Kernel(HeckeMatrix-e));
			print "Basis of eigenvectors computed for eigenvalue";
			print(e);
			print(B);
			for j in B do
				eigvec:=HMFZero(Mk, 1*ZF, [2*i, 2*i, 2*i]);
				for l in [1..NumberOfColumns(j)] do	
					eigvec+:=j[l]*CuspBas[i][l];
				end for;
				if Coefficients(eigvec)[2] ne 0 then
					newf:=1/Coefficients(eigvec)[2]*eigvec;
				else
					newf:= eigvec;
				end if;
				if CoefficientField(eigvec) eq Rationals()  then
					print(1);
					Append(~NewBas[i], newf);
					Append(~BasTrack[i], newf);
					Append(~sizes, CheckRoots(newf: Prec:=6, Emb:=1));
				else
					Append(~BasTrack[i], newf);
					C:=CompEmb(newf: Bound:=30);
					print(Degree(Parent(Coefficients(newf)[4])));
					for x in [1..#C] do
						count:=0;
						for y in NewBas[i] do
							if Type(y) eq SeqEnum then
								if L2Dist(C[x],y) le 0.0001 then
									count+:=1;
								end if;
							end if;
						end for;
						if count eq 0 then
							Append(~NewBas[i], C[x]);
							Append(~sizes, CheckRoots(newf: Prec:=6, Emb:=x));
						end if;
					end for;
				end if;
			end for;
		end for;
	else
		Append(~NewBas, []);
		Append(~BasTrack, []);
		Append(~PerPolys, []);
	end if;
end for;


Cputime(t);























//D=81;
//Second cubic field
t:=Cputime();
F := NumberField(x^3-3*x-1);
ZF:=Integers(F);
sizes:=[];


ids:=IdealsUpTo(150, F);
prec := 35000;                        
M := HMFSpace(F, prec);
H := HeckeCharacterGroup(1*ZF);
Y2:=EisensteinSeries(M, 1*ZF,  H!1, H!1, [2,2,2]);
Y4:=EisensteinSeries(M, 1*ZF,  H!1, H!1, [4,4,4]);
Y6:=EisensteinSeries(M, 1*ZF,  H!1, H!1, [6,6,6]);
Y8:=EisensteinSeries(M, 1*ZF,  H!1, H!1, [8,8,8]);
XX:=Y2^2;


prec2:=16464;                         
M2 := HMFSpace(F, prec2);
E2:=EisensteinSeries(M2, 1*ZF,  H!1, H!1, [2,2,2]);
H42:=HeckeOperator(E2^2, ids[3]);
Y42:=1/Coefficients(H42)[1]*H42;
H62:=HeckeOperator(E2^3, ids[3]);
Y62:=1/Coefficients(H62)[1]*H62;




prec1:=2058;                           
M1 := HMFSpace(F, prec1);
Y21:=EisensteinSeries(M1, 1*ZF,  H!1, H!1, [2,2,2]);
Y41:=EisensteinSeries(M1, 1*ZF,  H!1, H!1, [4,4,4]);
Y61:=EisensteinSeries(M1, 1*ZF,  H!1, H!1, [6,6,6]);
Y81:=EisensteinSeries(M1, 1*ZF,  H!1, H!1, [8,8,8]);
H82:=HeckeOperator(XX^2, ids[5]);
Y82:=1/Coefficients(H82)[1]*H82;
H83:=HeckeOperator(XX^2, ids[6]);
Y83:=1/Coefficients(H83)[1]*H83;






///Basis of monomials

MonBas:=[**];
P<x21, x41, x42, x61, x62, x81, x82, x83> := PolynomialRing(RationalField(), [2,4,4,6,6,8,8,8]);
for i in [1..4] do
	Append(~MonBas, []);
	Mons:=MonomialsOfWeightedDegree(P, 2*i);
	for j in Mons do
		print(j);
		factor:=Factorization(j);
		mult:=1;
		for l in factor do
			if l[1] eq x21 then
				mult:=mult*Y21^l[2];
			end if;
			if l[1] eq x41 then
				mult:=mult*Y41^l[2];
			end if;
			if l[1] eq x42 then
				mult:=mult*Y42^l[2];
			end if;
			if l[1] eq x61 then
				mult:=mult*Y61^l[2];
			end if;
			if l[1] eq x62 then
				mult:=mult*Y62^l[2];
			end if;
			if l[1] eq x81 then
			  	mult:=mult*Y81^l[2];
			end if;
			if l[1] eq x82 then
			  	mult:=mult*Y82^l[2];
			end if;
			if l[1] eq x83 then
			  	mult:=mult*Y83^l[2];
			end if;
		end for;
		Append(~MonBas[i], mult);
	end for;
end for;



//Refine to sup basis

CuspBas:=[**];
for i:=1 to 4 do
	print(2*i);
	Append(~CuspBas, []);
	for j:=1 to #MonBas[i] do
		print([2*i,j]);
		if i ge 2 then
			for j:=1 to #MonBas[i]-1 do
				//print([2*i,j]);
				for l:=j+1 to #MonBas[i] do
					if #CuspBas[i] eq 0 then
						Append(~CuspBas[i], MonBas[i][j]-MonBas[i][l]);
						print([2*i, j,l,Coefficients(MonBas[i][j]-MonBas[i][l])[4],1]);
					else
						//print(LinearDependence([CoeffsUpToPrec(x, 100): x in CuspBas[i]] cat [CoeffsUpToPrec(MonBas[i][j]-MonBas[i][l], 100)]));
						if  #LinearDependence([CoeffsUpToPrec(x, 50): x in CuspBas[i]] cat [CoeffsUpToPrec(MonBas[i][j]-MonBas[i][l], 50)]) eq 0 then
							Append(~CuspBas[i], MonBas[i][j]-MonBas[i][l]);
							print([2*i, j,l,Coefficients(MonBas[i][j]-MonBas[i][l])[4], 2]);
						end if;
					end if;
				end for;
			end for;
		end if;
	end for;
end for;



//Forther refine to eigenform basis. Compute period polynomials.
NewBas:=[**];
BasTrack:=[**];
PerPolys:=[**];
for i in [1..#CuspBas] do
	if #CuspBas[i] ne 0 then
		Mk:=HMFSpace(F, 2000);
		Append(~NewBas,[**]);
		Append(~BasTrack,[**]);
		Append(~PerPolys,[**]);
		tempor:=[];
		CoeffLi:=[];
		HeckeOpAct:=[CoeffsUpToPrec(HeckeOperator(p, ids[5]), 50): p in CuspBas[i]];
		LinDeps:=[];
		for a in CuspBas[i] do
			Lin:=LinearDependence(HeckeOpAct cat [CoeffsUpToPrec(a, 50)]);
			Append(~LinDeps, #Lin);
			Append(~CoeffLi, [Lin[1][r]/Lin[1][#Lin[1]] : r in [1..#Lin[1]-1]]);
		end for;
		//print(LinDeps);
		H:=Matrix(CoeffLi);
		//print(H);
		//Determinant(H);
		newroots:=[r[1]: r in Factorization(CharacteristicPolynomial(H))];
		K:=ext< Rationals() | newroots>;
		HeckeMatrix:=Matrix(K, CoeffLi);
		print([2*i, #LinDeps]);
		Eig:=[t[1] : t in Eigenvalues(HeckeMatrix)];
		print "Eigenvalues computed";
		print(#Eig);
		for e in Eig do
			B:= Basis(Kernel(HeckeMatrix-e));
			print "Basis of eigenvectors computed for eigenvalue";
			print(e);
			print(B);
			for j in B do
				eigvec:=HMFZero(Mk, 1*ZF, [2*i, 2*i, 2*i]);
				for l in [1..NumberOfColumns(j)] do	
					eigvec+:=j[l]*CuspBas[i][l];
				end for;
				if Coefficients(eigvec)[2] ne 0 then
					newf:=1/Coefficients(eigvec)[2]*eigvec;
				else
					newf:= eigvec;
				end if;
				if CoefficientField(eigvec) eq Rationals()  then
					print(1);
					Append(~NewBas[i], newf);
					Append(~BasTrack[i], newf);
					Append(~sizes, CheckRoots(newf: Prec:=6, Emb:=1));
				else
					Append(~BasTrack[i], newf);
					C:=CompEmb(newf: Bound:=30);
					print(Degree(Parent(Coefficients(newf)[4])));
					for x in [1..#C] do
						count:=0;
						for y in NewBas[i] do
							if Type(y) eq SeqEnum then
								if L2Dist(C[x],y) le 0.0001 then
									count+:=1;
								end if;
							end if;
						end for;
						if count eq 0 then
							Append(~NewBas[i], C[x]);
							Append(~sizes, CheckRoots(newf: Prec:=6, Emb:=x));
						end if;
					end for;
				end if;
			end for;
		end for;
	else
		Append(~NewBas, []);
		Append(~BasTrack, []);
		Append(~PerPolys, []);
	end if;
end for;


Cputime(t);



