//Tb=T2;
fvScalarMatrix TEqnb
(
   -fvm::div(-phi,Tb) 
   +fvm::laplacian(DT, Tb)
 ==
   RHS_Tb /*Q*/  // 于是 Q
);
TEqnb.solve();
