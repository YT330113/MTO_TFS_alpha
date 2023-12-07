T=T2;
fvScalarMatrix TEqn
(
    fvm::div(phi,T)
   -fvm::laplacian(DT, T)
 ==
    Q
);
TEqn.solve();

// while (simple.correctNonOrthogonal())
// {
//          fvScalarMatrix TEqn
//          (
//              fvm::div(phi,T)
//             -fvm::laplacian(DT, T)
//           ==
//              Q
//          );
//          TEqn.relax();
//          TEqn.solve();
//          fvOptions.correct(T);
// }

// 输出温度场
if (runTime.writeTime())
{
  T.write();
}