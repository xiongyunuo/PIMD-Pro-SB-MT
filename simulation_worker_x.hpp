#ifndef SIMULATION_WORKER_X_HPP
#define SIMULATION_WORKER_X_HPP

#include "particle_x.hpp"
#include <cmath>
#include <atomic>
//#include <shared_mutex>

extern std::atomic_int XCount;
extern std::atomic_bool XGo1;
extern std::atomic_bool XGo2;
extern std::atomic_bool XGo3;
extern std::atomic_bool XGo4;
extern std::atomic_bool XGo5;
extern std::atomic_bool XGo6;
extern std::atomic_bool XIni;
extern std::atomic_bool XRun;

//extern std::shared_mutex XReadMutex;

class XSimulationWorker {
public:
  XParticle *particles;
  XParticle *particles2;
  int N, P;
  XNum omgP;
  XNum L;
  XNum *ENkCache2;
  int NN;
  XNum *ENkCache3;
  XNum **ENkCache;
  XNum *VBCache;
  XNum *ForceVBCache;
  XNum *ForceCache;
  int Nb;
  XParticle *particlesb;
  XParticle center;
  XNum omg0;
  XNum omgz;
  XNum g, s;
  int NT;
  XNum vi;
  XNum beta;
  XNum h;
  XNum energy;
  bool ok;
  XSimulationWorker() { ok = true; }
  void Update1();
  void Update2();
  void NHForce(XNum m, int f, XNum *v, XNum *Q, XNum *vtheta, XNum *res);
  int PrevIndex(int l, int j, int N2, int k);
  void RelativeDistance(XParticle *p1, XParticle *p2, XNum *displace);
  XNum XExp(XNum k, XNum E, XNum EE);
  XNum XMinE(int N2);
  XNum ENk2(int N2);
  void FillENk2();
  void FillENk();
  void StartWorker();
  inline int Index(int l, int j) { return (l-1)*P+j-1; }
  XNum Distance(XParticle *p1, XParticle *p2);
  inline XNum MinimumImage1(XNum a);
  inline XNum MinimumImage2(XNum a);
  int NextIndex(int l, int j, int N2, int k);
  void dENk(int N2, int k, int l, int j, XNum *res);
  void FillForceVB();
  void FillForce();
  XNum TrapEnergy(int index);
  XNum PairEnergy(int index, int index2);
  void TrapForce(int index, XNum *res);
  void PairForce(int index, int index2, XNum *res);
  void PeriodBoundary();
};

XNum XSimulationWorker::MinimumImage1(XNum a) {
  if (std::abs(a) > L / 2)
    return L - std::abs(a);
  return a;
}

XNum XSimulationWorker::MinimumImage2(XNum a) {
  if (std::abs(a) > L / 2) {
    if (a < 0)
      return L - std::abs(a);
    else
      return -(L - std::abs(a));
  }
  return a;
}

#endif