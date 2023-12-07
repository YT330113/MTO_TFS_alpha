Info<< "\nslove linear elasticity equation\n" << endl;

initialResidual = 0;
iCorr=0;
do
{
  {
      fvVectorMatrix DEqn
      (
          fvm::d2dt2(D)
       ==
          fvm::laplacian(2*mu + lambda, D, "laplacian(DD,D)")
        + divSigmaExp
        - fvc::grad(threeKalpha*T)
      );

      initialResidual = DEqn.solve().max().initialResidual();

      if (!compactNormalStress)
      {
          divSigmaExp = fvc::div(DEqn.flux());
      }
  }

  {
      volTensorField gradD(fvc::grad(D));
      sigmaD = mu*twoSymm(gradD) + (lambda*I)*tr(gradD);

      if (compactNormalStress)
      {
          divSigmaExp = fvc::div
          (
              sigmaD - (2*mu + lambda)*gradD,
              "div(sigmaD)"
          );
      }
      else
      {
          divSigmaExp += fvc::div(sigmaD);
      }
  }

} while (initialResidual > convergenceTolerance && ++iCorr < nCorr);

// 计算应力应变
gradD = fvc::grad(D);
sigmaD = mu * twoSymm(gradD) + (lambda * I) * tr(gradD);

// 输出位移及等效应力
if (runTime.writeTime())
{
   D.write();

   volSymmTensorField sigma
    (
        IOobject
        (
            "sigma",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        rho*sigmaD
    );

    if (thermalStress)
    {
        sigma = sigma - I*(rho*threeKalpha*T);
    }

    volScalarField sigmaEq
    (
        IOobject
        (
            "sigmaEq",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sqrt((3.0/2.0)*magSqr(dev(sigma)))
    );
    sigmaEq.write();
}

