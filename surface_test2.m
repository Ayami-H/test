(*PRB77 ver. unitcell*)
Clear["Global`*"]
mmat = 32;
dim = 2 mmat;(*numbers of atom*)

ho := SparseArray[
  {
   Band[{1, 2}, {-1, -1}] -> -t,
   Band[{2, 1}, {-1, -1}] -> -t,
   Band[{1, 3}, {-1, -1}, {4, 4}] -> I so,
   Band[{2, 4}, {-1, -1}, {4, 4}] -> -I so,
   Band[{3, 5}, {-1, -1}, {4, 4}] -> -I so,
   Band[{4, 6}, {-1, -1}, {4, 4}] -> I so,
   Band[{3, 1}, {-1, -1}, {4, 4}] -> -I so,
   Band[{4, 2}, {-1, -1}, {4, 4}] -> I so,
   Band[{5, 3}, {-1, -1}, {4, 4}] -> I so,
   Band[{6, 4}, {-1, -1}, {4, 4}] -> -I so
   }
  , {dim, dim}];

ht := SparseArray[
   {
    Band[{3, 4}, {-1, -1}, {4, 4}] -> -t,
    Band[{2, 1}, {-1, -1}, {4, 4}] -> -t,
    Band[{1, 1}, {-1, -1}] -> -2 I so,
    Band[{2, 4}, {-1, -1}, {4, 4}] -> I so,
    Band[{3, 5}, {-1, -1}, {4, 4}] -> I so,
    Band[{3, 1}, {-1, -1}, {4, 4}] -> I so,
    Band[{6, 4}, {-1, -1}, {4, 4}] -> I so
    }
   , {dim, dim}];
t = 2.66; a = 1; so = 0.0001;
z = 0.0002 I;

itmatl := Inverse[ht] // Normal;
ig0l := (z IdentityMatrix[dim] - ho).itmatl // Normal;

genZL[z_] := Module[{},
   mzl = SparseArray[{}, 2 dim {1, 1}];
   mzl[[1 ;; dim, dim + 1 ;; 2 dim]] = itmatl;
   mzl[[dim + 1 ;; 2 dim, 1 ;; dim]] = -ht\[Conjugate]\[Transpose];
   mzl[[dim + 1 ;; 2 dim, dim + 1 ;; 2 dim]] = ig0l;
   ];

itmatr := Inverse[ht\[Conjugate]\[Transpose]] // Normal;
ig0r := (z IdentityMatrix[dim] - ho).itmatr // Normal;

genZR[z_] := Module[{},
   mzr = SparseArray[{}, 2 dim {1, 1}];
   mzr[[1 ;; dim, dim + 1 ;; 2 dim]] = itmatr;
   mzr[[dim + 1 ;; 2 dim, 1 ;; dim]] = -ht;
   mzr[[dim + 1 ;; 2 dim, dim + 1 ;; 2 dim]] = ig0r;
   ];

igd := z IdentityMatrix[dim] - ho // Normal;
(*-------------------------------------------------------------*)
(*DOS*)
ParallelEvaluate[
  Off[Eigenvectors::arh];
  Off[Eigensystem::arh];
  Off[Inverse::luc];
  ];
SetSharedVariable[dat];
dat = {};
ParallelDo[
  genZL[z];
  genZR[z];
  igd[z];
  evecsl = Eigenvectors[mzl, dim];
  Ul12 = Take[evecsl, {1, dim}, {1, dim}]\[Transpose];
  Ul22 = Take[evecsl, {1, dim}, {dim + 1, 2 dim}]\[Transpose];
  gsemil = Ul12.Inverse[Ul22];
  \[Sigma]l = ht\[Conjugate]\[Transpose].gsemil.ht;
  evecsr = Eigenvectors[mzr, dim];
  Ur12 = Take[evecsr, {1, dim}, {1, dim}]\[Transpose];
  Ur22 = Take[evecsr, {1, dim}, {dim + 1, 2 dim}]\[Transpose];
  gsemir = Ur12.Inverse[Ur22];
  \[Sigma]r = ht.gsemir.ht\[Conjugate]\[Transpose];
  hd = igd - \[Sigma]l - \[Sigma]r;
  ghd = Inverse[hd];
  dos = -(1/\[Pi]) 1/dim Im[Sum[ghd[[m, m]], {m, 1, dim}]];
  AppendTo[dat, {Re[z]/t, dos}];
  Export["zig_1_trans_dos_1.dat", dat // Sort];
  , {z, -t + 0.0002 I, t + 0.0002 I, 0.01}];
(*-------------------------------------------------------------*)
(*transport*)
ParallelEvaluate[
  Off[Eigenvectors::arh];
  Off[Eigensystem::arh];
  Off[Inverse::luc];
  ];
SetSharedVariable[dat77];
dat77 = {};
ParallelDo[
  genZL[z];
  genZR[z];
  igd[z];
  evecsl = Eigenvectors[mzl, dim];
  Ul12 = Take[evecsl, {1, dim}, {1, dim}]\[Transpose];
  Ul22 = Take[evecsl, {1, dim}, {dim + 1, 2 dim}]\[Transpose];
  gsemil = Ul12.Inverse[Ul22];
  sil = ht\[Conjugate]\[Transpose].gsemil.ht;
  gaml = I (sil - sil\[Conjugate]\[Transpose]);
  evecsr = Eigenvectors[mzr, dim];
  Ur12 = Take[evecsr, {1, dim}, {1, dim}]\[Transpose];
  Ur22 = Take[evecsr, {1, dim}, {dim + 1, 2 dim}]\[Transpose];
  gsemir = Ur12.Inverse[Ur22];
  sir = ht.gsemir.ht\[Conjugate]\[Transpose];
  gamr = I (sir - sir\[Conjugate]\[Transpose]);
  hd = igd - sil - sir;
  ghd = Inverse[hd];
  tra = gaml.ghd.gamr.ghd\[Conjugate]\[Transpose];
  trans = Sum[tra[[m, m]], {m, 1, dim}] // Chop;
  AppendTo[dat77, {Re[z]/t, trans}];
  Export["zig_1_trans_conduct_1.dat", dat77 // Sort];
  , {z, -t + 0.0002 I, t + 0.0002 I, 0.01}];
