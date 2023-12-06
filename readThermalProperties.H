Info<< "Reading thermal properties\n" << endl;
IOdictionary thermalProperties
(
    IOobject
    (
        "thermalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar rhoc(thermalProperties.lookup("rhoc"));//rho_fluid*heat_capacity
dimensionedScalar kf(thermalProperties.lookup("kf"));//Thermal Conductivity of fluid
dimensionedScalar ks(thermalProperties.lookup("ks"));//Thermal Conductivity of solid
dimensionedScalar Q("Q", dimensionSet(0,0,-1,1,0,0,0),10000000.0/4180000);

volScalarField DT
(
    IOobject
    (
        "DT",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    (ks+(kf-ks)*xh*(1+qu)/(qu+xh))/rhoc
);
Info<< "Reading field T\n" << endl;
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
Info<< "Reading field Tb\n" << endl;
   volScalarField Tb
   (
       IOobject
       (
           "Tb",
           runTime.timeName(),
           mesh,
           IOobject::MUST_READ,
           IOobject::AUTO_WRITE
       ),
       mesh
   ); 
volScalarField T2(T);
volScalarField fsensMeanT
(
    IOobject
    (
        "fsensMeanT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -alphaMax*(1+qu)*qu/((qu+xh)*(qu+xh))*(U&Ub),
    zeroGradientFvPatchScalarField::typeName
);
volScalarField fsenshMeanT(fsensMeanT);

// 热应力相关
dimensionedScalar Talpha(thermalProperties.lookup("Taplha"));
volScalarField threeKalpha
(
    IOobject
    (
        "threeKalpha",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimensionSet(0, 2, -2 , -1, 0), 0.0)
);
threeKalpha = threeK*Talpha;