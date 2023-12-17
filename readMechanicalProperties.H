Info<< "Reading mechanical properties\n" << endl;
IOdictionary mechanicalProperties
(
    IOobject
    (
        "mechanicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);
int iCorr = 0, nCorr = 200;
scalar convergenceTolerance(5e-6), initialResidual(0);

Switch planeStress(mechanicalProperties.lookup("planeStress"));
dimensionedScalar rhoE("rhoE", dimMass / dimLength / sqr(dimTime), readScalar(mechanicalProperties.lookup("rhoE")));
dimensionedScalar rho("rho", dimMass / dimVolume, readScalar(mechanicalProperties.lookup("rho")));
scalar Po(readScalar(mechanicalProperties.lookup("Po")));
dimensionedScalar E(rhoE / rho);
dimensionedScalar Emin(E * 1e-9);
double qd = 0.1;
volScalarField mu(qd * (1 - xh) / (qd + xh) * (E - Emin) / (2.0 * (1.0 + Po)) + Emin / (2.0 * (1.0 + Po)));
volScalarField lambda(qd * (1 - xh) / (qd + xh) * Po * (E - Emin) / ((1.0 + Po) * (1.0 - 2.0 * Po)) + Po * Emin / ((1.0 + Po) * (1.0 - 2.0 * Po)));

volScalarField threeK(3 * lambda + 2 * mu);   // 热膨胀相关

Info<< "Reading field D\n" << endl;
volVectorField D
(
    IOobject
    (
        "D",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


// gradD 和 sigmaD 用于计算柔度
volTensorField gradD
(
    IOobject
    (
        "gradD",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
      fvc::grad(D)
);
volSymmTensorField sigmaD // 注意，sigmaD不是真的应力，差个密度
(
    IOobject
    (
        "sigmaD",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mu*twoSymm(gradD) + (lambda*I)*tr(gradD)
);

// explicit part of div(sigma)
volVectorField divSigmaExp
(
    IOobject
    (
        "divSigmaExp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::div(sigmaD)
);
Switch compactNormalStress("yes");
if (compactNormalStress)
{
    divSigmaExp -= fvc::laplacian(2*mu + lambda, D, "laplacian(DD,D)");
}
else
{
    divSigmaExp -= fvc::div((2*mu + lambda)*fvc::grad(D), "div(sigmaD)");
}
mesh.setFluxRequired(D.name());

// 柔度的灵敏度
volScalarField gsensCompliance
(
    IOobject
    (
        "gsensCompliance",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -gradD&&(rho*(3*xh*xh*(E-Emin)/(2.0*(1.0 + Po))*twoSymm(gradD))),  // 这里量纲不对，要给成 0 2 -3 0
    zeroGradientFvPatchScalarField::typeName
);
volScalarField gsenshCompliance(gsensCompliance);


// 热应力伴随
volVectorField Da  // 柔度的伴随位移
(
    IOobject
    (
        "Da",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
volTensorField gradDa
(
    IOobject
    (
        "gradDa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
      fvc::grad(Da)
);
volSymmTensorField sigmaDa 
(
    IOobject
    (
        "sigmaDa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mu*twoSymm(gradDa) + (lambda*I)*tr(gradDa)
);

// explicit part of div(sigmaDa)
volVectorField divSigmaExp_Da
(
    IOobject
    (
        "divSigmaExp_Da",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::div(sigmaDa)
);
if (compactNormalStress)
{
    divSigmaExp_Da -= fvc::laplacian(2*mu + lambda, Da, "laplacian(DD,Da)");
}
else
{
    divSigmaExp_Da -= fvc::div((2*mu + lambda)*fvc::grad(Da), "div(sigmaDa)");
}
mesh.setFluxRequired(Da.name());
