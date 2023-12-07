fvScalarMatrix TEqnc
(
   -fvm::div(-phi,Tc) 
   +fvm::laplacian(DT, Tc)
 ==
    (Talpha*(3*lambda+2*mu))*fvc::div(Da)  // 如果量纲有问题可以直接改在字典中 Tc 的量纲
);
TEqnc.solve();