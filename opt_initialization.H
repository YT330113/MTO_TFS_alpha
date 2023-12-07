Info<< "Initialization\n" << endl;
std::remove("Voluse.txt");
std::remove("Time.txt");
std::remove("MeanT.txt"); 
std::remove("PowerDiss.txt"); 
std::remove("Compliance.txt"); 
std::remove("lisan.txt");

// 一些常量
int i = 0, opt = 0, rex;         // 迭代计数
double del, d, dg;               // del就是heaviside过滤中的控制参数beta，d，dg用来统计xh的分布情况
scalar eta0, eta5, eta1, y0, y5; // heaviside 过滤二分法相关
scalar n = x.size();             // n:length of x in the i-th core
scalar N = x.size();
reduce(N, sumOp<scalar>());      // N:design freedom
scalar area = gSum(mesh.V());                                        // area 计算域总体积
scalar C(0), V(0), MeanT(0), PowerDiss(0);                           // 目标及约束值
dimensionedScalar b("b", dimensionSet(0, -2, 0, 0, 0, 0, 0), 1.0);   // filter radius of the PDE filter
dimensionedScalar b1("b1", dimensionSet(0, 2, -2, -2, 0, 0, 0), 1.0);
double len = std::pow(area / N, 1.0 / 3); // 估算单元平均尺寸
b.value() = 1.0 / (filterR * len / 3.464) / (filterR * len / 3.464); // PDE过滤参数

dimensionedScalar RHS_Tb("RHS_Tb", dimensionSet(0,0,-1,1,0,0,0),1/area/rhoc.value()); // 平均温度的伴随温度方程的右端项，为了匹配量纲

// MMA算法相关
double move = movlim; // move值用于放松，随迭代变化，movlim是固定的
int m = 3;            // number of constraint
MMA mma(n, m);
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
