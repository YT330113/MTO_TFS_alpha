// 计算二分法待求函数值
scalar diff(volScalarField &gamma,const scalarField &V,double del,double eta,int n) // gamma就是xp，
{
     int i;
     scalar z=0;
     double *x =new double[n];
     
     // 先计算出截断值eta下的xh
     for(i=0;i<n;i++)
     {
        if(gamma[i]<=eta)
        {
          x[i]=eta*(Foam::exp(-del*(1-gamma[i]/eta))-(1-gamma[i]/eta)*Foam::exp(-del));
        }
        else
        {
          x[i]=eta+(1-eta)*(1-Foam::exp(-del*(gamma[i]-eta)/(1-eta))+(gamma[i]-eta)*Foam::exp(-del)/(1-eta));
        }
     }
     // 计算二分法待求函数值，函数是宝体积heaviside过滤论文的式(21)
     for(i=0;i<n;i++)
     {
        z=z+(gamma[i]-x[i])*V[i];
     }
     delete x;
     return {z};
}

