Info<<"cost function"<<endl;

// 应变能约束，energy = sigma && epsilon
volScalarField Energy(0.25 * (gradD + gradD.T()) && sigmaD);
C = fvc::domainIntegrate(Energy).value() / C0; // domainIntegrate是全局积分

// 平均温度目标
MeanT = fvc::domainIntegrate(T).value() / area;

// 体积约束
V = fvc::domainIntegrate(xh).value() / area - (voluse /*+ 0* Foam::max(0, volinit - opt * 0.01)*/);

// 能量耗散约束
PowerDiss = 0;
for (i=0; i<nObjPatch; i++)
{
     PowerDiss = PowerDiss - sum(
     phi.boundaryField()[conPatchList[i]]*(p.boundaryField()[conPatchList[i]]
     + 0.5*magSqr(U.boundaryField()[conPatchList[i]]))); 
}
PowerDiss = PowerDiss / PowerDiss0;
if (GeoDim == 2)
{
   PowerDiss = PowerDiss / std::pow(mesh.V()[0], 1.0 / 3);
}
reduce(PowerDiss, sumOp<scalar>()); 

// 目标及约束值输出
if (Pstream::master())
{
   ofstream outfile1("Voluse.txt", std::ios::app);
   outfile1 << V << "\n";
   outfile1.close();

   ofstream outfile2("PowerDiss.txt", std::ios::app);
   outfile2 << PowerDiss << "\n";
   outfile2.close();

   ofstream outfile3("MeanT.txt", std::ios::app);
   outfile3 << MeanT << "\n";
   outfile3.close();

   ofstream outfile4("Time.txt", std::ios::app);
   outfile4 << runTime.elapsedClockTime() << "\n";
   outfile4.close();

   ofstream outfile5("Compliance.txt", std::ios::app);
   outfile5 << C << "\n";
   outfile5.close();
}

