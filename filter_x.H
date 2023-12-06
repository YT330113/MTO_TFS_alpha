
//************************* PDE filter ********************************//

fvScalarMatrix Eqn4(fvm::laplacian(xp)-fvm::Sp(b,xp)+x*b);
Eqn4.solve();

//************************ Heaviside_filter ***************************//

// del就是heaviside过滤中的控制参数beta，beta随着迭代次数增大而增大。
del=Foam::min(0.3*opt,100); 

// 二分法求解heaviside过滤截止参数eta
eta0=0.0001;  // eta左边界
eta1=1.0;     // eta右边界
y0=diff(xp,mesh.V(),del,eta0,n); // 计算eta=0左边界时的函数值
reduce(y0, sumOp<scalar>());
// 迭代求解保体积参数eta
do
{
  eta5=(eta0+eta1)/2.0;            // // 计算eta区间中点值
  y5=diff(xp,mesh.V(),del,eta5,n); // 计算eta区间中点的函数值
  reduce(y5, sumOp<scalar>());
  if(y0*y5<0)
  { 
     eta1=eta5;
  }
  else
  { 
     eta0=eta5;
     y0=y5;
  }
} while ((eta1-eta0)>0.0001);
// eta5=0.5;

d=0; // 用来统计算例收敛的情况

// 根据保体积heaviside过滤论文式(19)计算xh
for(i=0;i<n;i++)
{
   if(xp[i]<=eta5)
   {
     xh[i]=eta5*(Foam::exp(-del*(1-xp[i]/eta5))-(1-xp[i]/eta5)*Foam::exp(-del)); 
   }
   else
   {
     xh[i]=eta5+(1-eta5)*(1-Foam::exp(-del*(xp[i]-eta5)/(1-eta5))+(xp[i]-eta5)*Foam::exp(-del)/(1-eta5));
   }
   if (xh[i] > 0.95 || xh[i] < 0.05)
     d++;
}

// 统计密度值收敛情况
MPI_Allreduce(&d, &dg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); // 计算全局的d
dg=dg/N; // 百分占比，有dg%的单元密度达到了xh[i] > 0.95 || xh[i] < 0.05
if (Pstream::master())
{ 
   ofstream outfile1("lisan.txt",std::ios::app);
   outfile1<<dg<<"\n";
   outfile1.close();
}