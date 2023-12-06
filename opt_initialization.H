Info<< "Initialization\n" << endl;
std::remove("Voluse.txt");
std::remove("Time.txt");
std::remove("MeanT.txt"); 
std::remove("PowerDiss.txt"); 
std::remove("Compliance.txt"); 
std::remove("lisan.txt"); 

int i=0,opt=0,rex;
double del,d,dg;// del就是heaviside过滤中的控制参数beta，d，dg用来统计xh的分布情况
scalar n(0),N(0),eta0,eta5,eta1,y0,y5;
n=x.size();//n:length of x in the i-th core
N=x.size();//N:design freedom
reduce(N, sumOp<scalar>()); 

scalar C(0),V(0),area=gSum(mesh.V()),MeanT(0),PowerDiss(0);
double len=std::pow(area/N, 1.0/3);
b.value()=1.0/(filterR*len/3.464)/(filterR*len/3.464);
double move=movlim;
//***********************************//MMA
int m=3;//number of constraint
MMA mma(n,m);
mma.asyminit = mma_init;
mma.asymdec = mma_dec;
mma.asyminc = mma_inc;
mma.raa0 = raa0;
mma.movelimit = movlim;
std::vector<double> xmma(n);
std::vector<double> dfdx(n);
std::vector<double> gx(m);
std::vector<std::vector<double>> dgdx(m);
for (int i = 0; i < m; i++)
{
  dgdx[i].resize(n);
}