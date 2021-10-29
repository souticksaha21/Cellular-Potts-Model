double R = 10; // um, cell radius (half side length -- area = 4*R^2)
double a = 2; // um, lattice size
double tau = 100; // s, MC time step
double d = a; // um, pixel depth
double T = 9*60*60+1; // s, total time
double deltat = 15*60; // s, measurement interval
double alpha = 2; // 1/um, cell-ECM line energy
double lambdaA = 1e-2; // 1/um^4, area deviation penalty
int Nrecep = 5000;
//double c0 = 2.5*.6; // 1/um^-3, concentration at x = 0
double epsilon = 56.23; // 1/um, bias strength
int Z = 1000; // number of trials
double r = 1/tau; //1/deltat; // 1/s, inverse memory time
double eta = 107.15; // 1/um, persistence strength

char base[] = "../../dat/cpm5g";
