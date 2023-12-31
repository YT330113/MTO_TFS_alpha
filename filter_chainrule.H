
//不可设计域的灵敏度置零
if(solid_area)
{
   forAll(cells_solid, i)
   {
      fsenshMeanT[cells_solid[i]]=0;
      gsenshPowerDiss[cells_solid[i]]=0;
      gsenshCompliance[cells_solid[i]]=0;
      gsenshVol[cells_solid[i]]=0;
   }
}
if(fluid_area)
{
   forAll(cells_fluid, i)
   {
      fsenshMeanT[cells_fluid[i]]=0;
      gsenshPowerDiss[cells_fluid[i]]=0;
      gsenshCompliance[cells_fluid[i]]=0;
      gsenshVol[cells_fluid[i]]=0;
   }
}

// //灵敏度归一化
// if (opt<300)
// {
// fsenshMeanT=fsenshMeanT/fvc::domainIntegrate(mag(fsenshMeanT)).value()*area/100000;
// gsenshPowerDiss=gsenshPowerDiss/fvc::domainIntegrate(mag(gsenshPowerDiss)).value()*area/100000;
// gsenshCompliance=gsenshCompliance/fvc::domainIntegrate(mag(gsenshCompliance)).value()*area/100000;
// gsenshVol=gsenshVol/100000;
// }
// else
// {
gsenshPowerDiss=gsenshPowerDiss/gMax(mag(gsenshPowerDiss.primitiveField()));
gsenshCompliance=gsenshCompliance/gMax(mag(gsenshCompliance.primitiveField()))/100;
gsenshVol=gsenshVol/gMax(mag(gsenshVol.primitiveField()))/10000;
fsenshMeanT=fsenshMeanT/gMax(mag(fsenshMeanT.primitiveField()));
// }

if (Pstream::master())
{
   ofstream outfile1("Sensitivity.txt", std::ios::app);
   outfile1 << opt << '\t' << fsenshMeanT[222]
                   << '\t' << gsenshPowerDiss[222]
                   << '\t' << gsenshCompliance[222]
                   << '\t' << gsenshVol[222]
                   << "\n";
   outfile1.close();
}

//Heaviside过滤的链式求导
if(opt>1)
{
  for(i=0;i<n;i++)
  {
     if(xp[i]<=eta5)
     {
       drho[i]=del*Foam::exp(-del*(1-xp[i]/eta5))+Foam::exp(-del); 
     }
     else
     {
       drho[i]=del*Foam::exp(-del*(xp[i]-eta5)/(1-eta5))+Foam::exp(-del); 
     } 
  }
  fsenshMeanT=fsenshMeanT*drho;
  gsenshPowerDiss=gsenshPowerDiss*drho;
  gsenshCompliance=gsenshCompliance*drho;
  gsenshVol=gsenshVol*drho;
}

//PDE 过滤的链式求导
solve(fvm::laplacian(fsensMeanT)-fvm::Sp(b,fsensMeanT)+fsenshMeanT*b);

solve(fvm::laplacian(gsensPowerDiss)-fvm::Sp(b,gsensPowerDiss)+gsenshPowerDiss*b);

solve(fvm::laplacian(gsensVol)-fvm::Sp(b,gsensVol)+gsenshVol*b);

solve(fvm::laplacian(gsensCompliance)-fvm::Sp(b,gsensCompliance)+gsenshCompliance*b);

fsensMeanT.primitiveFieldRef()=fsensMeanT.primitiveFieldRef()*mesh.V()/gMax(mesh.V());
gsensPowerDiss.primitiveFieldRef()=gsensPowerDiss.primitiveFieldRef()*mesh.V()/gMax(mesh.V());
gsensCompliance.primitiveFieldRef()=gsensCompliance.primitiveFieldRef()*mesh.V()/gMax(mesh.V());
gsensVol.primitiveFieldRef()=gsensVol.primitiveFieldRef()*mesh.V()/gMax(mesh.V());