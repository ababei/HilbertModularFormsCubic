SetDebugOnError(true);
load "config.m";

F := QuadraticField(13);
ZF<w> := Integers(F);
N := ideal<ZF | {1}>;
k := [2, 2];
K := Rationals();
prec := 100;
M := HMFSpace(F,  prec);

time g, r := ConstructGeneratorsAndRelations(M, N, 2, 8);
time g1, r1 := Relations(M, g, r, 20);

X := MakeScheme(g1, r1);

print X;
