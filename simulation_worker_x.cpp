#include "simulation_worker_x.hpp"
#include <chrono>
#include <thread>

std::atomic_int XCount;
std::atomic_bool XGo1;
std::atomic_bool XGo2;
std::atomic_bool XGo3;
std::atomic_bool XGo4;
std::atomic_bool XGo5;
std::atomic_bool XGo6;
std::atomic_bool XIni;
std::atomic_bool XRun;

//std::shared_mutex XReadMutex;

XNum XSimulationWorker::ENk2(int N2) {
  int j;
  XNum res = 0;
  for (j = 1; j < P; ++j) {
    int index = Index(N2, j);
    int index2 = Index(N2, j+1);
    res += 0.5*particles[index].m*omgP*omgP*Distance(&particles[index], &particles[index2]);
  }
  return res;
}

void XSimulationWorker::FillENk2() {
  int l;
  for (l = 1; l <= N; ++l)
    ENkCache2[l-1] = ENk2(l);
}

XNum XSimulationWorker::Distance(XParticle *p1, XParticle *p2) {
  XNum res = 0;
  int i;
  for (i = 0; i < D; ++i) {
    XNum a = MinimumImage1(p1->coor[i]-p2->coor[i]);
    res += a*a;
  }
  return res;
}

void XSimulationWorker::StartWorker() {
  while (XRun.load()) {
    while (!XGo1.load()) {
      //std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
    FillENk2();
    XCount.fetch_sub(1, std::memory_order_relaxed);
    while (!XGo3.load()) {
      //std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
    FillENk();
    XCount.fetch_sub(1, std::memory_order_relaxed);
    while (!XGo4.load()) {
      //std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
    FillForceVB();
    XCount.fetch_sub(1, std::memory_order_relaxed);
    while (!XGo5.load()) {
      //std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
    FillForce();
    XCount.fetch_sub(1, std::memory_order_relaxed);
    while (!XGo6.load()) {
      //std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
    if (XIni.load())
      Update1();
    else
      Update2();
    XCount.fetch_sub(1, std::memory_order_relaxed);
    while (!XGo2.load()) {
      //std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
    XCount.fetch_sub(1, std::memory_order_relaxed);
  }
}

void XSimulationWorker::FillENk() {
  int N = NN+XSimulationWorker::N;
  int N2, k, l;
  for (N2 = NN; N2 <= N; ++N2)
    for (k = 1; k <= N2; ++k) {
      XNum res = 0;
      //XReadMutex.lock_shared();
      for (l = N2-k+1; l <= N2; ++l) {
        res += ENkCache3[l-1];
        int index = Index(l, P);
        int index2 = NextIndex(l, P, N2, k);
        res += 0.5*particles2[index].m*omgP*omgP*Distance(&particles2[index], &particles2[index2]);
      }
      ENkCache[N2-1][k-1] = res;
      //XReadMutex.unlock_shared();
    }
}

int XSimulationWorker::NextIndex(int l, int j, int N2, int k) {
  int res = Index(l,j+1);
  if (j == P) {
    if (l == N2)
      res = Index(N2-k+1,1);
    else
      res = Index(l+1,1);
  }
  return res;
}

void XSimulationWorker::FillForceVB() {
  int N2, k, l, j;
  XNum *res = new XNum[D*(NT+1)];
  XNum grad[D];
  for (l = 1; l <= N; ++l)
    for (j = 1; j <= P; ++j) {
      int i;
      for (i = 0; i < D; ++i)
        res[i] = 0;
      for (N2 = 1; N2 <= NT; ++N2) {
        //XReadMutex.lock_shared();
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
            dENk(N2, k, NN+l, j, grad);
            sum += (grad[i]+res[D*(N2-k)+i])*XExp(k,ENkCache[N2-1][k-1]+VBCache[N2-k],tmp);
          }
          res[D*N2+i] = sum / sum2;
        }
        //XReadMutex.unlock_shared();
      }
      int index = Index(l, j);
      for (i = 0; i < D; ++i)
        ForceVBCache[D*index+i] = res[D*NT+i];
    }
  delete[] res;
}

XNum XSimulationWorker::XExp(XNum k, XNum E, XNum EE) {
  if (vi == 0)
    return std::exp(-beta*E+EE);
  else
    return std::exp((k-1)*std::log(vi)-beta*E+EE);
}

XNum XSimulationWorker::XMinE(int N2) {
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

void XSimulationWorker::dENk(int N2, int k, int l, int j, XNum *res) {
  int i;
  for (i = 0; i < D; ++i)
    res[i] = 0;
  if (l >= N2 - k + 1 && l <= N2) {
    int index = Index(l, j);
    int index2 = NextIndex(l, j, N2, k);
    int index3 = PrevIndex(l, j, N2, k);
    XNum displace[D];
    RelativeDistance(&particles2[index], &particles2[index2], displace);
    for (i = 0; i < D; ++i)
      res[i] += particles2[index].m*omgP*omgP*displace[i];
    RelativeDistance(&particles2[index], &particles2[index3], displace);
    for (i = 0; i < D; ++i)
      res[i] += particles2[index].m*omgP*omgP*displace[i];
  }
}

void XSimulationWorker::RelativeDistance(XParticle *p1, XParticle *p2, XNum *displace) {
  int i;
  for (i = 0; i < D; ++i)
    displace[i] = MinimumImage2(p1->coor[i]-p2->coor[i]);
}

int XSimulationWorker::PrevIndex(int l, int j, int N2, int k) {
  int res = Index(l,j-1);
  if (j == 1) {
    if (l == N2 - k + 1)
      res = Index(N2,P);
    else
      res = Index(l-1,P);
  }
  return res;
}

void XSimulationWorker::FillForce() {
  energy = P*D*N/(2*beta);
  int j, l, k;
  for (l = 1; l <= N; ++l)
    for (j = 1; j <= P; ++j) {
      int index = Index(l, j);
      int i;
      for (i = 0; i < D; ++i)
        ForceCache[D*index+i] = -ForceVBCache[D*index+i]/particles[index].m;
    }
  XNum trap[D];
  for (l = 1; l <= N; ++l)
    for (j = 1; j <= P; ++j) {
      int index = Index(l, j);
      TrapForce(index, trap);
      energy += TrapEnergy(index)/P;
      int i;
      for (i = 0; i < D; ++i)
        ForceCache[D*index+i] += trap[i]/particles[index].m/P;
    }
  XNum inter[D];
  /*for (l = 1; l <= N; ++l)
    for (j = 1; j <= P; ++j) {
      int index = Index(l, j);
      //XReadMutex.lock_shared();
      for (k = 1; k <= NT; ++k) {
        if (NN+l == k) continue;
        int index2 = Index(k, j);
        int index3 = Index(NN+l, j);
        PairForce(index3, index2, inter);
        energy += PairEnergy(index3, index2)/P;
        int i;
        for (i = 0; i < D; ++i)
          ForceCache[D*index+i] += inter[i]/particles[index].m/P;
      }
      //XReadMutex.unlock_shared();
    }*/
  for (l = 1; l <= N; ++l)
    for (j = 1; j <= P; ++j) {
      int index = Index(l, j);
      //XReadMutex.lock_shared();
      for (k = 1; k <= Nb; ++k) {
        int index2 = Index(NT+k, j);
        int index3 = Index(NN+l, j);
        PairForce(index3, index2, inter);
        energy += PairEnergy(index3, index2)/P;
        int i;
        for (i = 0; i < D; ++i)
          ForceCache[D*index+i] += inter[i]/particles[index].m/P;
      }
      //XReadMutex.unlock_shared();
    }
}

void XSimulationWorker::TrapForce(int index, XNum *res) {
  int i;
  XNum displace[D];
  RelativeDistance(&particles[index], &center, displace);
  for (i = 0; i < D; ++i) {
    if (i == D-1)
      res[i] = -particles[index].m*omgz*omgz*displace[i];
    else
      res[i] = -particles[index].m*omg0*omg0*displace[i];
  }
}

void XSimulationWorker::PairForce(int index, int index2, XNum *res) {
  int i;
  XNum displace[D];
  XNum inter;
  if (index2 < P*NT) {
    inter = (g/(M_PI*s*s))*std::exp(-Distance(&particles2[index], &particles2[index2])/(s*s));
    RelativeDistance(&particles2[index], &particles2[index2], displace);
  }
  else {
    inter = (g/(M_PI*s*s))*std::exp(-Distance(&particles2[index], &particlesb[index2-P*NT])/(s*s));
    RelativeDistance(&particles2[index], &particlesb[index2-P*NT], displace);
  }
  for (i = 0; i < D; ++i)
    res[i] = ((2*displace[i])/(s*s))*inter;
}

void XSimulationWorker::Update1() {
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
  PeriodBoundary();
}

void XSimulationWorker::Update2() {
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
  Update1();
}

void XSimulationWorker::NHForce(XNum m, int f, XNum *v, XNum *Q, XNum *vtheta, XNum *res) {
  XNum sum = 0;
  int i;
  for (i = 0; i < f; ++i)
    sum += m*v[i]*v[i];
  res[0] = (sum-f*(1/beta))/Q[0];
  for (i = 1; i < M; ++i)
    res[i] = (Q[i-1]*vtheta[i-1]*vtheta[i-1]-(1/beta))/Q[i];
}

void XSimulationWorker::PeriodBoundary() {
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

XNum XSimulationWorker::TrapEnergy(int index) {
  int i;
  XNum displace[D];
  XNum res = 0;
  RelativeDistance(&particles[index], &center, displace);
  for (i = 0; i < D; ++i) {
    if (i == D-1)
      res += 0.5*particles[index].m*omgz*omgz*displace[i]*displace[i];
    else
      res += 0.5*particles[index].m*omg0*omg0*displace[i]*displace[i];
  }
  return res;
  //return 0.5*particles[index].m*omg0*omg0*Distance(&particles[index], &center);
}

XNum XSimulationWorker::PairEnergy(int index, int index2) {
  XNum inter;
  if (index2 < P*NT)
    inter = (g/(M_PI*s*s))*std::exp(-Distance(&particles2[index], &particles2[index2])/(s*s));
  else
    inter = (g/(M_PI*s*s))*std::exp(-Distance(&particles2[index], &particlesb[index2-P*NT])/(s*s));
  return 0.5*inter;
}