Info<< "\nslove AdjD_Complicande\n" << endl;

initialResidual = 0;    // 收敛准则暂时设置和D一样
iCorr=0;
do
{
  {
      fvVectorMatrix DaEqn
      (
          fvm::d2dt2(Da)
       ==
          fvm::laplacian(2*mu + lambda, Da, "laplacian(DD,Da)")
        + divSigmaExp_Da
      );

      initialResidual = DaEqn.solve().max().initialResidual();

      if (!compactNormalStress)
      {
          divSigmaExp_Da = fvc::div(DaEqn.flux());
      }
  }

  {
      volTensorField gradDa(fvc::grad(Da));
      sigmaDa = mu*twoSymm(gradDa) + (lambda*I)*tr(gradDa);

      if (compactNormalStress)
      {
          divSigmaExp_Da = fvc::div
          (
              sigmaDa - (2*mu + lambda)*gradDa,
              "div(sigmaDa)"
          );
      }
      else
      {
          divSigmaExp_Da += fvc::div(sigmaDa);
      }
  }

} while (initialResidual > convergenceTolerance && ++iCorr < nCorr);
