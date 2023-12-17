
Info<<"sensitivity analysis"<<endl;
//计算灵敏度
gsenshPowerDiss=-alphaMax*(1+qu)*qu/((qu+xh)*(qu+xh))*(U&Ua);  // 这里少乘一个单元体积
fsenshMeanT=-alphaMax*(1+qu)*qu/((qu+xh)*(qu+xh))*(U&Ub)+b1*(kf-ks)*(1+qu)*qu/((qu+xh)*(qu+xh))*(fvc::grad(T)&fvc::grad(Tb))/rhoc; // 也没有乘单元体积
mu=-(1+qd)*qd/((qd+xh)*(qd+xh))*(E-Emin)/(2.0*(1.0+Po));
lambda=-(1+qd)*qd/((qd+xh)*(qd+xh))*Po*(E-Emin)/((1.0+Po)*(1.0-2.0*Po));
gsenshCompliance=-alphaMax*(1+qu)*qu/((qu+xh)*(qu+xh))*(U&Uc)
   +b1*(kf-ks)*(1+qu)*qu/((qu+xh)*(qu+xh))*(fvc::grad(T)&fvc::grad(Tc))/rhoc
   -b3_c*(gradDa && (rho*(mu*twoSymm(gradD)+ (lambda*I)*tr(gradD))));
gsenshVol.primitiveFieldRef() = 1.0;

//链式求导计算对物理密度的灵敏度
#include "filter_chainrule.H"

//往MMA中传入数据
gx[0]=1*(PowerDiss/Foam::max(PowerDissRelax-opt*0.4,PowerDissMax)-1.0);
gx[1]=1*(C/Foam::max(CRelax-opt*0.2,CMax)-1.0);
gx[2]=1*V;

for(i=0;i<n;i++)
{
 xmma[i] = x[i];
 dfdx[i]=fsensMeanT[i]/1e3;//sensitivities of objective function
 dgdx[0][i]=gsensPowerDiss[i]/1e3;//sensitivities of constraint function
 dgdx[1][i]=gsensCompliance[i]/1e3;//sensitivities of constraint function
 dgdx[2][i]=gsensVol[i]/1e3;//sensitivities of constraint function
}

Info<< "\n run MMA \n" << endl;
// mma.movelimit = move;
mma.MMAsolver(xmma, dfdx, gx, dgdx);
for(i=0;i<n;i++)
{
   x[i]=xmma[i];
}
//PDE过滤和Heaviside过滤计算物理密度
#include "filter_x.H"
