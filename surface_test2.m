(*PRB 77, 085408 (2008) 'Quantum conductance of graphene nanoribbons with edge defects' *)
Clear["Global`*"]
mmat = 32;(*nanoribbon width*)
dim = 2 mmat;(*numbers of atom*)

ho := SparseArray[                              (*疎な配列*)(*On-site matrix H00*)
  {
   Band[{1, 2}, {-1, -1}] -> -t,                (*帯対角行列*)
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

ht := SparseArray[                                (*Transfer matrix H01*)
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
t = 2.66; a = 1; so = 0.0001;                     (*so:second nearest neighbot hopping term, do not make zero matrix*)
z = 0.0002 I;                                     (*z+in : n=0.0002 infinitesimal*)

itmatl := Inverse[ht] // Normal;
ig0l := (z IdentityMatrix[dim] - ho).itmatl // Normal;

genZL[z_] := Module[{},                           (*Left lead*)
   mzl = SparseArray[{}, 2 dim {1, 1}];
   mzl[[1 ;; dim, dim + 1 ;; 2 dim]] = itmatl;
   mzl[[dim + 1 ;; 2 dim, 1 ;; dim]] = -Transpose[Conjugate[ht]];
   mzl[[dim + 1 ;; 2 dim, dim + 1 ;; 2 dim]] = ig0l;
   ];

itmatr := Inverse[Transpose[Conjugate[ht]]] // Normal;
ig0r := (z IdentityMatrix[dim] - ho).itmatr // Normal;

genZR[z_] := Module[{},                           (*Right lead*)
   mzr = SparseArray[{}, 2 dim {1, 1}];
   mzr[[1 ;; dim, dim + 1 ;; 2 dim]] = itmatr;
   mzr[[dim + 1 ;; 2 dim, 1 ;; dim]] = -ht;
   mzr[[dim + 1 ;; 2 dim, dim + 1 ;; 2 dim]] = ig0r;
   ];

igd := z IdentityMatrix[dim] - ho // Normal;      (*Central lead*)
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
  evecsl = Eigenvectors[mzl, dim];                              (*固有ベクトル、固有値の大きい順からdim個*)
  Ul12 = Take[evecsl, {1, dim}, {1, dim}]\[Transpose];          (*eveclからdim個取ってきて、{{a1,a2,...,a_dim},{b1,b2,...,b_dim},...}を*)
  Ul22 = Take[evecsl, {1, dim}, {dim + 1, 2 dim}]\[Transpose];  (*{{a1,b1,...,[dim]_1},{a2,b2,...,[dim]2},...}と並べ替える。*)
  gsemil = Ul12.Inverse[Ul22];                                  (*表面グリーン関数(Left lead)*)
  \[Sigma]l = Transpose[Conjugate[ht]].gsemil.ht;               (*自己エネルギー*)
  evecsr = Eigenvectors[mzr, dim];
  Ur12 = Take[evecsr, {1, dim}, {1, dim}]\[Transpose];
  Ur22 = Take[evecsr, {1, dim}, {dim + 1, 2 dim}]\[Transpose];
  gsemir = Ur12.Inverse[Ur22];                                  (*表面グリーン関数(Right lead)*)
  \[Sigma]r = ht.gsemir.Transpose[Conjugate[ht]];
  hd = igd - \[Sigma]l - \[Sigma]r;
  ghd = Inverse[hd];                                            (*デバイスのグリーン関数(central lead)*)
  dos = -(1/\[Pi]) 1/dim Im[Sum[ghd[[m, m]], {m, 1, dim}]];     (*状態密度*)
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
  gaml = I (sil - sil\[Conjugate]\[Transpose]);                     (*広がり行列*)
  evecsr = Eigenvectors[mzr, dim];
  Ur12 = Take[evecsr, {1, dim}, {1, dim}]\[Transpose];
  Ur22 = Take[evecsr, {1, dim}, {dim + 1, 2 dim}]\[Transpose];
  gsemir = Ur12.Inverse[Ur22];
  sir = ht.gsemir.ht\[Conjugate]\[Transpose];
  gamr = I (sir - sir\[Conjugate]\[Transpose]);                     (*広がり行列*)
  hd = igd - sil - sir;
  ghd = Inverse[hd];
  tra = gaml.ghd.gamr.ghd\[Conjugate]\[Transpose];
  trans = Sum[tra[[m, m]], {m, 1, dim}] // Chop;                    (*トランスミッション*)
  AppendTo[dat77, {Re[z]/t, trans}];
  Export["zig_1_trans_conduct_1.dat", dat77 // Sort];
  , {z, -t + 0.0002 I, t + 0.0002 I, 0.01}];
