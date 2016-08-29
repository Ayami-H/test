(*m file*)
Clear["Global`*"]
t = 2.66; ho = 0; h1 = -t;
tau1 = -h1;
tau2 = tau1\[Conjugate];

wmat[z_] :=
  Module[{},
   ws = z - ho;
   wb = z - ho;
   tau1 = -h1;
   tau2 = tau1\[Conjugate];
   wbi = 1/wb;
   For[k = 1, k < 20, k++,
    ws = ws - tau1 wbi tau2;
    wb = wb - tau1 wbi tau2 - tau2 wbi tau1;
    tau1 = -tau1 wbi tau1;
    tau2 = -tau2 wbi tau2;
    wbi = 1/wb;
    tval = Max[Abs[tau1], Abs[tau2]];
    If[tval < 1 10^-8, Break[]];
    ];
   ];

dat = {};
Do[
  wmat[z];
  inws = 1/ws;
  imgs = -(1/\[Pi]) Im[inws];
  AppendTo[dat, {Re[z], imgs}];
  , {z, -6 + 0.0001 I, 6 + 0.0001 I, 0.01}];
