#include "simulation_x.hpp"
#include "random_x.hpp"
#include <chrono>
#include <thread>

XNum XSimulation::Distance(XParticle *p1, XParticle *p2) {
  XNum res = 0;
  int i;
  for (i = 0; i < D; ++i) {
    XNum a = MinimumImage1(p1->coor[i]-p2->coor[i]);
    res += a*a;
  }
  return res;
}

void XSimulation::RelativeDistance(XParticle *p1, XParticle *p2, XNum *displace) {
  int i;
  for (i = 0; i < D; ++i)
    displace[i] = MinimumImage2(p1->coor[i]-p2->coor[i]);
}

void XSimulation::Initial(std::istream &in) {
  std::string parameter;
  while (in >> parameter) {
    if (parameter == "Na")
      in >> Na;
    else if (parameter == "Nb")
      in >> Nb;
    else if (parameter == "P")
      in >> P;
    else if (parameter == "T") {
      in >> T;
      beta = 1/(kB*T);
    }
    else if (parameter == "L")
      in >> L;
    else if (parameter == "seed") {
      int s;
      in >> s;
      XSetRandSeed(s);
    }
    else if (parameter == "vi")
      in >> vi;
    else if (parameter == "omg0")
      in >> omg0;
    else if (parameter == "omgz")
      in >> omgz;
    else if (parameter == "g")
      in >> g;
    else if (parameter == "s")
      in >> s;
    else if (parameter == "h")
      in >> h;
    else if (parameter == "skip")
      in >> skip;
    else if (parameter == "step")
      in >> step;
    else if (parameter == "den_num") {
      int num;
      in >> num;
      stat.Init(num);
    }
    else if (parameter == "thread_count")
      in >> thread_count;
  }
  N = Na + Nb;
  particles = new XParticle[N*P];
  int l, j;
  XNum vs[D];
  int i;
  for (i = 0; i < D; ++i)
    vs[i] = 0;
  for (l = 1; l <= N; ++l)
    for (j = 1; j <= P; ++j) {
      int index = Index(l, j);
      particles[index].m = 1;
      particles[index].n = l;
      particles[index].p = j;
      for (i = 0; i < D; ++i) {
        XNum tmp = XRandGauss() * std::sqrt(1 / (particles[index].m * beta));
        vs[i] += tmp;
        particles[index].vel[i] = tmp;
      }
      for (i = 0; i < D; ++i)
        particles[index].coor[i] = L / 2 + 1.0 * (XRandFloat()-0.5);
    }
  for (i = 0; i < D; ++i)
    vs[i] /= N*P;
  for (l = 1; l <= N; ++l)
    for (j = 1; j <= P; ++j) {
      int index = Index(l, j);
      for (i = 0; i < D; ++i)
        particles[index].vel[i] -= vs[i];
      int j;
      for (i = 0; i < D; ++i)
        for (j = 0; j < M; ++j) {
          particles[index].theta[i][j] = 1;
          particles[index].vtheta[i][j] = 1;
          particles[index].Q[i][j] = 0.1;
        }
    }
  VelocityRescale();
  omgP = std::sqrt(P)/(beta*hBar);
  ENkCache = new XNum*[N];
  for (l = 1; l <= N; ++l) {
    if (l <= Na)
      ENkCache[l-1] = new XNum[l];
    else
      ENkCache[l-1] = new XNum[l-Na];
  }
  VBCache = new XNum[N+2];
  ForceVBCache = new XNum[D*N*P];
  ENkCache2 = new XNum[N];
  for (i = 0; i < D; ++i)
    center.coor[i] = L/2;
  t = 0;
  ForceCache = NULL;
  count = 0;
  incre = L/2/stat.num;
  ok = true;
}

void XSimulation::Dump(std::ostream &out) {
  int l, j;
  for (l = 1; l <= N; ++l)
    for (j = 1; j <= P; ++j) {
      int index = Index(l, j);
      int i;
      out << particles[index].n << " " << particles[index].p << " " << particles[index].m << std::endl;
      for (i = 0; i < D; ++i)
        out << particles[index].coor[i] << " ";
      for (i = 0; i < D; ++i)
        out << particles[index].vel[i] << " ";
      out << std::endl;
      int j;
      for (i = 0; i < D; ++i)
        for (j = 0; j < M; ++j)
          out << particles[index].theta[i][j] << " " << particles[index].vtheta[i][j] << " " << particles[index].Q[i][j] << " ";
      out << std::endl << std::endl;
    }
}

XSimulation::~XSimulation() {
  delete[] particles;
  int l;
  for (l = 1; l <= N; ++l)
    delete[] ENkCache[l-1];
  delete[] ENkCache;
  delete[] VBCache;
  delete[] ForceVBCache;
  if (ForceCache)
    delete[] ForceCache;
  delete[] ENkCache2;
}

void XSimulation::VelocityRescale() {
  XNum sum = 0;
  int l, j;
  for (l = 1; l <= N; ++l)
    for (j = 1; j <= P; ++j) {
      int index = Index(l, j);
      int i;
      for (i = 0; i < D; ++i)
        sum += particles[index].m*particles[index].vel[i]*particles[index].vel[i];
    }
  XNum lam = std::sqrt((N*P)*D*T/sum);
  for (l = 1; l <= N; ++l)
    for (j = 1; j <= P; ++j) {
      int index = Index(l, j);
      int i;
      for (i = 0; i < D; ++i)
        particles[index].vel[i] *= lam;
    }
}

XNum XSimulation::ENk(int N2, int k) {
  XNum res = 0;
  int l, j;
  for (l = N2 - k + 1; l <= N2; ++l)
    for (j = 1; j <= P; ++j) {
      int index = Index(l, j);
      int index2 = NextIndex(l, j, N2, k);
      res += 0.5*particles[index].m*omgP*omgP*Distance(&particles[index], &particles[index2]);
    }
  return res;
}

int XSimulation::NextIndex(int l, int j, int N2, int k) {
  int res = Index(l,j+1);
  if (j == P) {
    if (l == N2)
      res = Index(N2-k+1,1);
    else
      res = Index(l+1,1);
  }
  return res;
}

int XSimulation::PrevIndex(int l, int j, int N2, int k) {
  int res = Index(l,j-1);
  if (j == 1) {
    if (l == N2 - k + 1)
      res = Index(N2,P);
    else
      res = Index(l-1,P);
  }
  return res;
}

XNum XSimulation::XExp(XNum k, XNum E, XNum EE) {
  if (vi == 0)
    return std::exp(-beta*E+EE);
  else
    return std::exp((k-1)*std::log(vi)-beta*E+EE);
}

XNum XSimulation::XMinE(int N2) {
  if (vi == 0)
    return beta*(ENkCache[N2-1][0]+VBCache[N2-1]);
  int k;
  XNum res = 100000000;
  for (k = 1; k <= N2; ++k) {
    XNum tmp = -(k-1)*std::log(vi)+beta*(ENkCache[N2-1][k-1]+VBCache[N2-k]);
    if (tmp < res)
      res = tmp;
  }
  return res;
}

void XSimulation::FillVB() {
  int N2, k;
  VBCache[0] = 0;
  for (N2 = 1; N2 <= N; ++N2) {
    XNum sum = 0;
    XNum tmp = XMinE(N2);
    for (k = 1; k <= N2; ++k) {
      if (vi == 0 && k-1 != 0)
        continue;
      sum += XExp(k,ENkCache[N2-1][k-1]+VBCache[N2-k],tmp);
    }
    VBCache[N2] = (tmp-std::log(sum)+std::log(N2))/beta;
  }
}

void XSimulation::dENk(int N2, int k, int l, int j, XNum *res) {
  int i;
  for (i = 0; i < D; ++i)
    res[i] = 0;
  if (l >= N2 - k + 1 && l <= N2) {
    int index = Index(l, j);
    int index2 = NextIndex(l, j, N2, k);
    int index3 = PrevIndex(l, j, N2, k);
    XNum displace[D];
    RelativeDistance(&particles[index], &particles[index2], displace);
    for (i = 0; i < D; ++i)
      res[i] += particles[index].m*omgP*omgP*displace[i];
    RelativeDistance(&particles[index], &particles[index3], displace);
    for (i = 0; i < D; ++i)
      res[i] += particles[index].m*omgP*omgP*displace[i];
  }
}

void XSimulation::FillForceVB() {
  int N2, k, l, j;
  XNum *res = new XNum[D*(N+1)];
  XNum grad[D];
  for (l = 1; l <= N; ++l)
    for (j = 1; j <= P; ++j) {
      int i;
      for (i = 0; i < D; ++i)
        res[i] = 0;
      for (N2 = 1; N2 <= N; ++N2) {
        XNum sum2 = 0;
        XNum tmp = XMinE(N2);
        for (k = 1; k <= N2; ++k) {
          if (vi == 0 && k-1 != 0)
            continue;
          sum2 += XExp(k,ENkCache[N2-1][k-1]+VBCache[N2-k],tmp);
        }
        for (i = 0; i < D; ++i) {
          XNum sum = 0;
          for (k = 1; k <= N2; ++k) {
            if (vi == 0 && k-1 != 0)
              continue;
            dENk(N2, k, l, j, grad);
            sum += (grad[i]+res[D*(N2-k)+i])*XExp(k,ENkCache[N2-1][k-1]+VBCache[N2-k],tmp);
          }
          res[D*N2+i] = sum / sum2;
        }
      }
      int index = Index(l, j);
      for (i = 0; i < D; ++i)
        ForceVBCache[D*index+i] = res[D*N+i];
    }
  delete[] res;
}

void XSimulation::FillENk() {
  int l, j;
  for (l = 1; l <= N; ++l) {
    for (j = 1; j <= l; ++j)
      ENkCache[l-1][j-1] = ENk(l, j);
  }
}

XNum XSimulation::ENk2(int N2) {
  int j;
  XNum res = 0;
  for (j = 1; j < P; ++j) {
    int index = Index(N2, j);
    int index2 = Index(N2, j+1);
    res += 0.5*particles[index].m*omgP*omgP*Distance(&particles[index], &particles[index2]);
  }
  return res;
}

void XSimulation::FillENk2() {
  int l;
  for (l = 1; l <= N; ++l)
    ENkCache2[l-1] = ENk2(l);
  int N2, k;
  for (N2 = 1; N2 <= N; ++N2)
    for (k = 1; k <= N2; ++k) {
      XNum res = 0;
      for (l = N2-k+1; l <= N2; ++l) {
        res += ENkCache2[l-1];
        int index = Index(l, P);
        int index2 = NextIndex(l, P, N2, k);
        res += 0.5*particles[index].m*omgP*omgP*Distance(&particles[index], &particles[index2]);
      }
      ENkCache[N2-1][k-1] = res;
    }
}

void XSimulation::TrapForce(int index, XNum *res) {
  int i;
  XNum displace[D];
  RelativeDistance(&particles[index], &center, displace);
  for (i = 0; i < D; ++i)
    res[i] = -particles[index].m*omg0*omg0*displace[i];
}

void XSimulation::PairForce(int index, int index2, XNum *res) {
  int i;
  XNum displace[D];
  XNum inter = (g/(M_PI*s*s))*std::exp(-Distance(&particles[index], &particles[index2])/(s*s));
  RelativeDistance(&particles[index], &particles[index2], displace);
  for (i = 0; i < D; ++i)
    res[i] = ((2*displace[i])/(s*s))*inter;
}

XNum *XSimulation::Force() {
  XNum *res = new XNum[D*N*P];
  //FillENk();
  FillENk2();
  FillVB();
  FillForceVB();
  int l, j, k;
  for (l = 1; l <= N; ++l)
    for (j = 1; j <= P; ++j) {
      int index = Index(l, j);
      int i;
      for (i = 0; i < D; ++i)
        res[D*index+i] = -ForceVBCache[D*index+i]/particles[index].m;
    }
  XNum trap[D];
  for (l = 1; l <= N; ++l)
    for (j = 1; j <= P; ++j) {
      int index = Index(l, j);
      TrapForce(index, trap);
      int i;
      for (i = 0; i < D; ++i)
        res[D*index+i] += trap[i]/particles[index].m/P;
    }
  XNum inter[D];
  for (l = 1; l <= N; ++l)
    for (j = 1; j <= P; ++j) {
      int index = Index(l, j);
      for (k = 1; k <= N; ++k) {
        if (l == k) continue;
        int index2 = Index(k, j);
        PairForce(index, index2, inter);
        int i;
        for (i = 0; i < D; ++i)
          res[D*index+i] += inter[i]/particles[index].m/P;
      }
    }
  return res;
}

void XSimulation::DumpForce(std::ostream &out) {
  int l, j;
  for (l = 1; l <= N; ++l) {
    for (j = 1; j <= l; ++j)
      out << ENkCache[l-1][j-1] << " ";
    out << std::endl;
  }
  out << std::endl;
  for (l = 0; l <= N; ++l)
    out << VBCache[l] << " ";
  out << std::endl;
  out << std::endl;
  for (l = 1; l <= N; ++l)
    for (j = 1; j <= P; ++j) {
      int index = Index(l, j);
      int i;
      for (i = 0; i < D; ++i)
        out << ForceVBCache[D*index+i] << " ";
    }
  out << std::endl;
  out << std::endl;
  XNum *f = Force();
  for (l = 1; l <= N; ++l)
    for (j = 1; j <= P; ++j) {
      int index = Index(l, j);
      int i;
      for (i = 0; i < D; ++i)
        out << f[D*index+i] << " ";
    }
  out << std::endl;
  out << std::endl;
  delete[] f;
}

XNum XSimulation::VBEnergy() {
  XNum *res = new XNum[N+1];
  res[0] = 0;
  int N2, k;
  for (N2 = 1; N2 <= N; ++N2) {
    XNum tmp = XMinE(N2);
    XNum sum2 = 0;
    for (k = 1; k <= N2; ++k) {
      if (vi == 0 && k-1 != 0)
        continue;
      sum2 += XExp(k,ENkCache[N2-1][k-1]+VBCache[N2-k],tmp);
    }
    XNum sum = 0;
    for (k = 1; k <= N2; ++k) {
      if (vi == 0 && k-1 != 0)
        continue;
      sum += (res[N2-k]-ENkCache[N2-1][k-1])*XExp(k,ENkCache[N2-1][k-1]+VBCache[N2-k],tmp);
    }
    res[N2] = sum / sum2;
  }
  XNum e = res[N];
  delete[] res;
  return e;
}

XNum XSimulation::TrapEnergy(int index) {
  return 0.5*particles[index].m*omg0*omg0*Distance(&particles[index], &center);
}

XNum XSimulation::PairEnergy(int index, int index2) {
  XNum inter = (g/(M_PI*s*s))*std::exp(-Distance(&particles[index], &particles[index2])/(s*s));
  return 0.5*inter;
}

XNum XSimulation::TotalEnergy() {
  XNum e = 0;
  e += P*D*N/(2*beta);
  e += VBEnergy();
  int l, j, k;
  for (l = 1; l <= N; ++l)
    for (j = 1; j <= P; ++j) {
      int index = Index(l, j);
      e += TrapEnergy(index)/P;
    }
  for (l = 1; l <= N; ++l)
    for (j = 1; j <= P; ++j) {
      int index = Index(l, j);
      for (k = 1; k <= N; ++k) {
        if (l == k) continue;
        int index2 = Index(k, j);
        e += PairEnergy(index, index2)/P;
      }
    }
  return e;
}

XNum XSimulation::TotalEnergyMT() {
  XNum e = 0;
  int tmpN = N;
  N = Na;
  e += VBEnergy();
  N = Nb;
  ENkCache = &ENkCache[Na];
  VBCache = &VBCache[Na+1];
  e += VBEnergy();
  N = tmpN;
  ENkCache = &ENkCache[-Na];
  VBCache = &VBCache[-Na-1];
  int i;
  for (i = 0; i < 2*thread_count; ++i)
    e += XWorkers[i]->energy;
  return e;
}

void XSimulation::DumpEnergy(std::ostream &out) {
  out << VBEnergy() << " " << TotalEnergy() << std::endl;
  out << std::endl;
}

void XSimulation::NHForce(XNum m, int f, XNum *v, XNum *Q, XNum *vtheta, XNum *res) {
  XNum sum = 0;
  int i;
  for (i = 0; i < f; ++i)
    sum += m*v[i]*v[i];
  res[0] = (sum-f*(1/beta))/Q[0];
  for (i = 1; i < M; ++i)
    res[i] = (Q[i-1]*vtheta[i-1]*vtheta[i-1]-(1/beta))/Q[i];
}

void XSimulation::UpdateMNHC_VV3() {
  if (ForceCache == NULL)
    ForceCache = Force();
  int j;
  XNum NHF[M];
  XNum v[1];
  for (j = 0; j < N*P; ++j) {
    int i;
    for (i = 0; i < D; ++i) {
      if (IsNan(ForceCache[D*j+i]) || IsInf(ForceCache[D*j+i])) {
        ok = false;
        return;
      }
      v[0] = particles[j].vel[i];
      NHForce(particles[j].m, 1, v, particles[j].Q[i], particles[j].vtheta[i], NHF);
      particles[j].vel[i] = particles[j].vel[i]*std::exp(-0.5*h*particles[j].vtheta[i][0])+0.5*h*ForceCache[D*j+i]*std::exp(-0.25*h*particles[j].vtheta[i][0]);
      int M2 = M/2;
      int k;
      for (k = 1; k <= M2; ++k)
        particles[j].theta[i][2*k-2] = particles[j].theta[i][2*k-2]+h*particles[j].vtheta[i][2*k-2]/2;
      for (k = 1; k <= M2; ++k)
        particles[j].vtheta[i][2*k-1] = particles[j].vtheta[i][2*k-1]*std::exp(-0.5*h*((k==M2)?0:particles[j].vtheta[i][2*k]))+0.5*h*NHF[2*k-1]*std::exp(-0.25*h*((k==M2)?0:particles[j].vtheta[i][2*k]));
      particles[j].coor[i] = particles[j].coor[i]+h*particles[j].vel[i];
      for (k = 1; k <= M2; ++k)
        particles[j].theta[i][2*k-1] = particles[j].theta[i][2*k-1]+h*particles[j].vtheta[i][2*k-1];
      v[0] = particles[j].vel[i];
      NHForce(particles[j].m, 1, v, particles[j].Q[i], particles[j].vtheta[i], NHF);
      for (k = 1; k <= M2; ++k)
        particles[j].vtheta[i][2*k-2] = particles[j].vtheta[i][2*k-2]*std::exp(-h*particles[j].vtheta[i][2*k-1])+h*NHF[2*k-2]*std::exp(-0.5*h*particles[j].vtheta[i][2*k-1]);
    }
  }
  t += h;
  delete[] ForceCache;
  ForceCache = Force();
  for (j = 0; j < N*P; ++j) {
    int i;
    for (i = 0; i < D; ++i) {
      if (IsNan(ForceCache[D*j+i]) || IsInf(ForceCache[D*j+i])) {
        ok = false;
        return;
      }
      int M2 = M/2;
      particles[j].vel[i] = particles[j].vel[i]*std::exp(-0.5*h*particles[j].vtheta[i][0])+0.5*h*ForceCache[D*j+i]*std::exp(-0.25*h*particles[j].vtheta[i][0]);
      int k;
      for (k = 1; k <= M2; ++k)
        particles[j].theta[i][2*k-2] = particles[j].theta[i][2*k-2]+h*particles[j].vtheta[i][2*k-2]/2;
      v[0] = particles[j].vel[i];
      NHForce(particles[j].m, 1, v, particles[j].Q[i], particles[j].vtheta[i], NHF);
      for (k = 1; k <= M2; ++k)
        particles[j].vtheta[i][2*k-1] = particles[j].vtheta[i][2*k-1]*std::exp(-0.5*h*((k==M2)?0:particles[j].vtheta[i][2*k]))+0.5*h*NHF[2*k-1]*std::exp(-0.25*h*((k==M2)?0:particles[j].vtheta[i][2*k]));
    }
  }
}

void XSimulation::PeriodBoundary() {
  int j;
  for (j = 0; j < N*P; ++j) {
    int i;
    for (i = 0; i < D; ++i) {
      if (particles[j].coor[i] < 0) {
        int n = (int)(std::abs(particles[j].coor[i]) / L);
        particles[j].coor[i] += (n + 1) * L;
      }
      else {
        int n = (int)(std::abs(particles[j].coor[i]) / L);
        particles[j].coor[i] -= n * L;
      }
    }
  }
}

XNum XSimulation::Temperature() {
  XNum sum = 0;
  int j;
  for (j = 0; j < N*P; ++j) {
    int i;
    for (i = 0; i < D; ++i)
      sum += particles[j].m*particles[j].vel[i]*particles[j].vel[i];
  }
  return sum/(D*N*P);
}

void XSimulation::MakeStat() {
  ++count;
  //stat.e += TotalEnergy();
  stat.e += TotalEnergyMT();
  XNum temp = Temperature();
  stat.T += temp;
  stat.dT += (temp-T)*(temp-T);
  /*int j;
  for (j = 0; j < N*P; ++j) {
    //XNum dis = std::sqrt(Distance(&particles[j], &center));
    XNum dis = std::abs(MinimumImage2(particles[j].coor[0]-L/2));
    int index = int(dis/incre);
    if (index < 0)
      index = 0;
    if (index >= stat.num)
      continue;
    stat.den[index] += 1;
  }*/
  int j, l;
  for (j = 1; j <= P; ++j) {
    XNum c = 0;
    for (l = 1; l <= N; ++l) {
      int index = Index(l, j);
      c += particles[index].coor[D-1];
    }
    c /= N;
    XNum dis = std::abs(MinimumImage2(c-L/2));
    /*if (dis > incre/3) continue;
    for (l = 1; l <= N; ++l) {
      int i = Index(l, j);
      XNum dis = std::abs(MinimumImage2(particles[i].coor[D-1]-c));
      int index = int(dis/incre);
      if (index < 0)
        index = 0;
      if (index >= stat.num)
        continue;
      stat.den[index] += 1;
    }*/
    int index = int(dis/incre);
    if (index < 0)
      index = 0;
    if (index >= stat.num)
      continue;
    stat.den[index] += 1;
  }
}

void XSimulation::DumpStat(std::ostream &out) {
  out << "E " << stat.e << std::endl;
  out << "T " << stat.T << std::endl;
  out << "dT " << stat.dT << " " << "ideal " << 2.0/(D*N*P) << std::endl;
  int i;
  for (i = 0; i < stat.num; ++i)
    out << "{" << i*incre << "," << stat.den[i] << "},";
  out << std::endl;
  out << std::endl;
}

void XSimulation::NormalizeStat() {
  stat.e /= count;
  stat.T /= count;
  stat.dT /= count;
  stat.dT /= T*T;
  int i;
  //for (i = 0; i < stat.num; ++i)
    //stat.den[i] /= std::pow((i+0.5)*incre,D-1);
  XNum norm = 0;
  for (i = 0; i < stat.num; ++i)
    norm += stat.den[i]*incre;
  for (i = 0; i < stat.num; ++i)
    stat.den[i] /= norm;
}

XSimulationWorker **XWorkers;

void XCreateWorker(int index, int N, int P, XNum omgP, XNum L, XParticle *particles, XNum *ENkCache2, XNum *ENkCache3, XNum **ENkCache, XParticle *particles2, int NN, XNum *VBCache, XNum *ForceVBCache, int NT, XNum vi, XNum beta, XParticle center, XNum omg0, XNum g, XNum s, XNum h, int Nb, XParticle *particlesb, XNum omgz) {
  XSimulationWorker *worker = new XSimulationWorker();
  XWorkers[index] = worker;
  worker->N = N;
  worker->P = P;
  worker->omgP = omgP;
  worker->L = L;
  worker->particles = particles;
  worker->ENkCache2 = ENkCache2;
  worker->ENkCache3 = ENkCache3;
  worker->ENkCache = ENkCache;
  worker->particles2 = particles2;
  worker->NN = NN;
  worker->VBCache = VBCache;
  worker->ForceVBCache = ForceVBCache;
  worker->NT = NT;
  worker->vi = vi;
  worker->beta = beta;
  worker->center = center;
  worker->omg0 = omg0;
  worker->g = g;
  worker->s = s;
  worker->h = h;
  worker->Nb = Nb;
  worker->particlesb = particlesb;
  worker->omgz = omgz;
  worker->StartWorker();
}

void XSimulation::SpamWorkers() {
  XRun.store(true);
  XGo1.store(false);
  XGo2.store(false);
  XGo3.store(false);
  XGo4.store(false);
  XGo5.store(false);
  XGo6.store(false);
  XIni.store(true);
  XCount.store(2*thread_count);
  int tmpN = N;
  N = Na;
  int N2 = N/thread_count;
  if (N % thread_count != 0)
    ++N2;
  XWorkers = new XSimulationWorker*[2*thread_count];
  int i;
  for (i = 0; i < thread_count; ++i) {
    int N3 = N2;
    if (i == thread_count-1)
      N3 = N-N2*i;
    std::thread t(XCreateWorker, i, N3, P, omgP, L, &particles[i*N2*P], &ENkCache2[i*N2], ENkCache2, ENkCache, particles, i*N2, VBCache, &ForceVBCache[D*P*i*N2], N, vi, beta, center, omg0, g, s, h, Nb, &particles[Na*P], omgz);
    t.detach();
  }
  ENkCache = &ENkCache[Na];
  VBCache = &VBCache[Na+1];
  ForceVBCache = &ForceVBCache[D*Na*P];
  particles = &particles[Na*P];
  ENkCache2 = &ENkCache2[Na];
  N = Nb;
  N2 = N/thread_count;
  if (N % thread_count != 0)
    ++N2;
  for (i = 0; i < thread_count; ++i) {
    int N3 = N2;
    if (i == thread_count-1)
      N3 = N-N2*i;
    std::thread t(XCreateWorker, Na+i, N3, P, omgP, L, &particles[i*N2*P], &ENkCache2[i*N2], ENkCache2, ENkCache, particles, i*N2, VBCache, &ForceVBCache[D*P*i*N2], N, vi, beta, center, omg0, g, s, h, Na, &particles[-Na*P], omgz);
    t.detach();
  }
  N = tmpN;
  ENkCache = &ENkCache[-Na];
  VBCache = &VBCache[-Na-1];
  ForceVBCache = &ForceVBCache[-D*Na*P];
  particles = &particles[-Na*P];
  ENkCache2 = &ENkCache2[-Na];
}

XNum *XSimulation::Force2() {
  XCount.store(2*thread_count);
  XGo1.store(true);
  while (XCount.load() != 0) {
    //std::this_thread::sleep_for(std::chrono::milliseconds(10));
  }
  XGo1.store(false);
  XCount.store(2*thread_count);
  XGo3.store(true);
  while (XCount.load() != 0) {
    //std::this_thread::sleep_for(std::chrono::milliseconds(10));
  }
  XGo3.store(false);
  int tmpN = N;
  N = Na;
  FillVB();
  ENkCache = &ENkCache[Na];
  VBCache = &VBCache[Na+1];
  ForceVBCache = &ForceVBCache[D*Na*P];
  particles = &particles[Na*P];
  N = Nb;
  FillVB();
  N = tmpN;
  ENkCache = &ENkCache[-Na];
  VBCache = &VBCache[-Na-1];
  ForceVBCache = &ForceVBCache[-D*Na*P];
  particles = &particles[-Na*P];
  XCount.store(2*thread_count);
  XGo4.store(true);
  while (XCount.load() != 0) {
    //std::this_thread::sleep_for(std::chrono::milliseconds(10));
  }
  XGo4.store(false);
  XNum *res = new XNum[D*N*P];
  int i;
  int N2 = Na/thread_count;
  if (Na % thread_count != 0)
    ++N2;
  for (i = 0; i < thread_count; ++i)
    XWorkers[i]->ForceCache = &res[D*P*i*N2];
  N2 = Nb/thread_count;
  if (Nb % thread_count != 0)
    ++N2;
  for (i = 0; i < thread_count; ++i)
    XWorkers[Na+i]->ForceCache = &res[D*P*(i*N2+Na)];
  XCount.store(2*thread_count);
  XGo5.store(true);
  while (XCount.load() != 0) {
    //std::this_thread::sleep_for(std::chrono::milliseconds(10));
  }
  XGo5.store(false);
  return res;
}

void XSimulation::UpdateMNHC_VV3_MT() {
  if (ForceCache == NULL) {
    ForceCache = Force2();
    XIni.store(true);
  }
  else {
    delete[] ForceCache;
    ForceCache = Force2();
    XIni.store(false);
  }
  XCount.store(2*thread_count);
  XGo6.store(true);
  while (XCount.load() != 0) {
    //std::this_thread::sleep_for(std::chrono::milliseconds(10));
  }
  XGo6.store(false);
  XCount.store(2*thread_count);
  XGo2.store(true);
  while (XCount.load() != 0) {
    //std::this_thread::sleep_for(std::chrono::milliseconds(10));
  }
  XGo2.store(false);
  int i;
  for (i = 0; i < 2*thread_count; ++i)
    if (!XWorkers[i]->ok)
      ok = false;
}