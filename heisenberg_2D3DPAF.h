/*****************************************************************************
* Copyright (c) 2018, Hiroaki Kadowaki (kadowaki@tmu.ac.jp)
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without 
* modification, are permitted provided that the following conditions are met:
*
* 1. Redistributions of source code must retain the above copyright notice, 
* this list of conditions and the following disclaimer.
*
* 2. Redistributions in binary form must reproduce the above copyright notice, 
* this list of conditions and the following disclaimer in the documentation and/or 
* other materials provided with the distribution.
*
* 3. Neither the name of the copyright holder nor the names of its contributors may 
* be used to endorse or promote products derived from this software without specific 
* prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
* IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, 
* INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT 
* NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
* PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
* WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
* ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
* POSSIBILITY OF SUCH DAMAGE.
* (See https://opensource.org/licenses/BSD-3-Clause )
*
***************************************************************************
* ALPS Project: Algorithms and Libraries for Physics Simulations
* ALPS Libraries
* Copyright (C) 1997-2012 by Synge Todo <wistaria@comp-phys.org>
*
* This software is part of the ALPS libraries, published under the ALPS
* Library License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Library License along with
* the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
* available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/
// heisenberg_2D3DPAF.h
// modified from   example/parapack/heisenberg/heisenberg.h

#ifndef PARAPACK_EXAMPLE_HEISENBERG_HEISENBERG_H
#define PARAPACK_EXAMPLE_HEISENBERG_HEISENBERG_H

#include <alps/parapack/worker.h>
#include <alps/random/uniform_on_sphere_n.h>
#include <boost/array.hpp>

class heisenberg_spin : public boost::array<double, 3> {
private:
  typedef boost::array<double, 3> super_type;
public:
  heisenberg_spin() { clear(); }
  explicit heisenberg_spin(std::size_t i) { clear(); }
  heisenberg_spin(double x, double y, double z) {
    super_type::operator[](0) = x;
    super_type::operator[](1) = y;
    super_type::operator[](2) = z;
  }
  void clear() { for (int i = 0; i < 3; ++i) super_type::operator[](i) = 0; }
  double operator*(heisenberg_spin const& v) const {
    double res = 0;
    for (int i = 0; i < 3; ++i) res += super_type::operator[](i) * v[i];
    return res;
  }
  heisenberg_spin& operator+=(heisenberg_spin const& v) {
    for (int i = 0; i < 3; ++i) super_type::operator[](i) += v[i];
    return *this;
  }
};

inline alps::ODump& operator<<(alps::ODump& od, heisenberg_spin const& v) {
  od << v[0] << v[1] << v[2];
  return od;
}

inline alps::IDump& operator>>(alps::IDump& id, heisenberg_spin& v) {
  id >> v[0] >> v[1] >> v[2];
  return id;
}

class heisenberg_worker : public alps::parapack::lattice_mc_worker<> {
private:
  typedef alps::parapack::lattice_mc_worker<> super_type;
  typedef heisenberg_spin spin_type;
  typedef alps::uniform_on_sphere_n<3, double, spin_type> dist_type;

public:
  heisenberg_worker(alps::Parameters const& params) : super_type(params), mcs_(params),
    spins_(num_sites()) {
    beta_ = (params.defined("T") ? 1 / evaluate("T", params) : 0);
    // table of coupling constants
    int num_types = 0;
    bond_iterator itr, itr_end;
    for (boost::tie(itr, itr_end) = bonds(); itr != itr_end; ++itr)
      num_types = std::max(num_types, int(bond_type(*itr)) + 1);
    coupling_.clear();
    coupling_.resize(num_types, (params.defined("J")) ? evaluate("J", params) : 1);
    for (int t = 0; t < num_types; ++t) {
      std::string name = std::string("J") + boost::lexical_cast<std::string>(t);
      if (params.defined(name)) coupling_[t] = evaluate(name, params);
    }
    // std::cerr << "coupling constants = " << alps::write_vector(coupling_) << std::endl;
// coupling_[0] = J  exchange energy for NN pair  =  J_{NN} sigma_z sigma_z
// coupling_[1] = J  J = -1 --> J_{NN} = 1
// coupling_[2] = J
//
//  Lattice size // hexagonal a-axis L = L_a_hex
//  Nspin = 4 * L * L * (3 * H)
    L_a_hex = (params.defined("L") ? evaluate("L", params) : 1);
// local z diections at 4 sites j = 0, 1, 2, 3
    double xx=1.0/sqrt(3.0);
    double zvec_init[4][3]={{xx,xx,xx},{xx,-xx,-xx},{-xx,xx,-xx}, {-xx,-xx,xx}};
    for(int i=0 ; i<4 ; ++i){
        for(int j=0 ; j<3 ; ++j){
            zvec[i][j]=zvec_init[i][j];
        }
    }
// (delta, q) parameters for the Hamiltonian ( Onoda & Tanaka Phys. Rev. B 83, 094411 (2011) )
    coupling_delta = (params.defined("delta") ? evaluate("delta", params) : 0);
    coupling_q = (params.defined("q") ? evaluate("q", params) : 0);
    xx = 2.0*(2.0/3.0)*(4.0*atan(1.0) );
    coupling_cos2phi[0]=cos(-xx); //cos(-2*2*pi/3)
    coupling_cos2phi[1]=cos( xx); //cos( 2*2*pi/3)
    coupling_cos2phi[2]=1.0; //cos( 0 );
    coupling_sin2phi[0]=sin(-xx); //sin(-2*2*pi/3)
    coupling_sin2phi[1]=sin( xx); //sin( 2*2*pi/3)
    coupling_sin2phi[2]=0.0; //sin( 0 );
//
// magnetic field h111 along [111] direction
// 1.0 T (Tb2Ti2O7 ~ 5 mu_B) -->  h111 = 3.35856 K
// 5 * (1 mu_B) * (1 T) / k_{B} = 5 * 9.27400915*10^(-24) / (1.3806488*10^(-23)) = 3.35856922847 K
    hh[0] = (params.defined("h111") ? evaluate("h111", params) : 0);
    hh[1] = hh[0]/sqrt(3.0);
    hh[0] = hh[1];
    hh[2] = hh[1];
//
    // random initial spins
    for (int s = 0; s < num_sites(); ++s) spins_[s] = dist_(engine());
//
    update_energy();
  }
  virtual ~heisenberg_worker() {}

  void init_observables(alps::Parameters const&, alps::ObservableSet& obs) {
    obs << alps::SimpleRealObservable("Temperature")
        << alps::SimpleRealObservable("Inverse Temperature")
        << alps::SimpleRealObservable("Number of Sites")
        << alps::SimpleRealObservable("Number of 2D Sites")
        << alps::RealObservable("Energy")
        << alps::RealObservable("Energy Density")
        << alps::RealObservable("Energy^2")
        << alps::RealObservable("Magnetization 111 abs")
        << alps::RealObservable("Magnetization 111^2")
        << alps::RealObservable("Magnetization 111^4")
        << alps::RealObservable("Mag quadrupole uu3D[0]")
        << alps::RealObservable("Mag quadrupole uu3D[1]")
        << alps::RealObservable("Mag quadrupole uu3D[2]")
        << alps::RealObservable("Mag quadrupole uu3D")
        << alps::RealObservable("Mag quadrupole uu3D E")
        << alps::RealObservable("Mag quadrupole uu3D^2")
        << alps::RealObservable("Mag quadrupole uu3D^4")
        << alps::RealObservable("Mag quadrupole uu2D[0]")
        << alps::RealObservable("Mag quadrupole uu2D[1]")
        << alps::RealObservable("Mag quadrupole uu2D[2]")
        << alps::RealObservable("Mag quadrupole uu2D")
        << alps::RealObservable("Mag quadrupole uu2D E")
        << alps::RealObservable("Mag quadrupole uu2D^2")
        << alps::RealObservable("Mag quadrupole uu2D^4");
  }

  bool is_thermalized() const { return mcs_.is_thermalized(); }
  double progress() const { return mcs_.progress(); }

  void run(alps::ObservableSet& obs) {
    ++mcs_;
    for (int s = 0; s < num_sites(); ++s) {
      spin_type spin_new = dist_(engine());
      double diff = 0;
// diff of Zeemen energy
      diff  = hh[0]*zvec[site_type(s)][0];
      diff += hh[1]*zvec[site_type(s)][1];
      diff += hh[2]*zvec[site_type(s)][2];
      diff = diff * (spins_[s][2] - spin_new[2]);
//
      neighbor_bond_iterator itr, itr_end;
      for (boost::tie(itr, itr_end) = neighbor_bonds(s); itr != itr_end; ++itr) {
        int t = s ^ source(*itr) ^ target(*itr); // calculate index of neighbor spin
//        diff += coupling_[bond_type(*itr)] * (spins_[s][2] * spins_[t][2] - spin_new[2] * spins_[t][2]);
//
// diff of exchange energy
        double xx,yy,zz,xy,yx;
        xx= spins_[s][0] * spins_[t][0] - spin_new[0] * spins_[t][0];
        yy= spins_[s][1] * spins_[t][1] - spin_new[1] * spins_[t][1];
        zz= spins_[s][2] * spins_[t][2] - spin_new[2] * spins_[t][2];
        xy= spins_[s][0] * spins_[t][1] - spin_new[0] * spins_[t][1];
        yx= spins_[s][1] * spins_[t][0] - spin_new[1] * spins_[t][0];
        diff += coupling_[0]*zz;
        diff += coupling_[0]*coupling_delta*(xx+yy);
        diff += coupling_[0]*coupling_q*(coupling_cos2phi[bond_type(*itr)]*(xx-yy)-coupling_sin2phi[bond_type(*itr)]*(xy+yx));
      }
      if (uniform_01() < std::exp(- beta_ * diff)) spins_[s] = spin_new;
    }
//
    double magvec[3]={0.0, 0.0, 0.0};
    for (int s = 0; s < num_sites(); ++s) {
        magvec[0] += zvec[site_type(s)][0] * spins_[s][2];
        magvec[1] += zvec[site_type(s)][1] * spins_[s][2];
        magvec[2] += zvec[site_type(s)][2] * spins_[s][2];
    }
// update magnetization along 111
    double mag111;
        mag111=(magvec[0]+magvec[1]+magvec[2])/(sqrt(3.0) * num_sites() );
    double mag111sq; mag111sq = mag111 * mag111 ;
// update mag_quadrupole3D[4][2]
    double mag_quadrupole3D[4][2]={ {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0} };
    for (int s = 0; s < num_sites(); ++s) {
        mag_quadrupole3D[site_type(s)][0] += spins_[s][0];
        mag_quadrupole3D[site_type(s)][1] += spins_[s][1];
    }
    for(int nu = 0; nu < 4; ++nu) {
      for(int al = 0; al < 2; ++al) {  mag_quadrupole3D[nu][al] *= 4.0/num_sites();  }
    }
//
    double xx,yy,zz;
    xx=sqrt(3.0)/4.0; yy=1.0/4.0;  zz=1.0/2.0;
    double vv[3][4][2]={ {{0.0, -zz}, {0.0,  zz}, {0.0,  zz}, {0.0, -zz}},
                         {{-xx,  yy}, { xx, -yy}, {-xx,  yy}, { xx, -yy}},
                         {{ xx,  yy}, { xx,  yy}, {-xx, -yy}, {-xx, -yy}} };
    double uu3D[3]={0.0, 0.0, 0.0};
      for(int nu = 0; nu < 4; ++nu) {
        for(int al = 0; al < 2; ++al) {
          uu3D[0] += mag_quadrupole3D[nu][al] * vv[0][nu][al];
          uu3D[1] += mag_quadrupole3D[nu][al] * vv[1][nu][al];
          uu3D[2] += mag_quadrupole3D[nu][al] * vv[2][nu][al];
        }
      }
      uu3D[0] = 0.5 * uu3D[0];
      uu3D[1] = 0.5 * uu3D[1];
      uu3D[2] = 0.5 * uu3D[2];
// update uusq
    double uu3Dsq;
      uu3Dsq = uu3D[0] * uu3D[0] + uu3D[1] * uu3D[1] + uu3D[2] * uu3D[2] ;
// update mag_quadrupole2D[4][2]
    double mag_quadrupole2D[4][2]={ {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0} };
    for (int s = 0; s < num_sites(); ++s) {
      if ( ( coordinate(s)[2] * sqrt(6.0) )  < 1.5) {
//           coordinate(s)[2]/(1.0/sqrt(6.0)) 
        mag_quadrupole2D[ site_type(s) ][0] += spins_[s][0];
        mag_quadrupole2D[ site_type(s) ][1] += spins_[s][1];
      }
    }
    for(int nu = 0; nu < 4; ++nu) {
      for(int al = 0; al < 2; ++al) { mag_quadrupole2D[nu][al] *= 1.0/(L_a_hex * L_a_hex); }
    }
//
    double uu2D[3]={0.0, 0.0, 0.0};
      for(int nu = 0; nu < 4; ++nu) {
        for(int al = 0; al < 2; ++al) {
          uu2D[0] += mag_quadrupole2D[nu][al] * vv[0][nu][al];
          uu2D[1] += mag_quadrupole2D[nu][al] * vv[1][nu][al];
          uu2D[2] += mag_quadrupole2D[nu][al] * vv[2][nu][al];
        }
      }
      uu2D[0] = 0.5 * uu2D[0];
      uu2D[1] = 0.5 * uu2D[1];
      uu2D[2] = 0.5 * uu2D[2];
// update uu2Dsq
    double uu2Dsq;
      uu2Dsq = (uu2D[0] + uu2D[1] + uu2D[2])*2.0/3.0 ;
      uu2Dsq = uu2Dsq * uu2Dsq ;
//
    update_energy();
    add_constant(obs["Temperature"], 1/beta_);
    add_constant(obs["Inverse Temperature"], beta_);
    add_constant(obs["Number of Sites"], (double)num_sites());
    add_constant(obs["Number of 2D Sites"], (L_a_hex * L_a_hex) );
    obs["Energy"] << energy_;
    obs["Energy Density"] << energy_ / num_sites();
    obs["Energy^2"] << energy_ * energy_;
    obs["Magnetization 111 abs"] << fabs(mag111) ;
    obs["Magnetization 111^2"] << mag111sq ;
    obs["Magnetization 111^4"] << mag111sq * mag111sq ;
    obs["Mag quadrupole uu3D[0]"] << uu3D[0] ;
    obs["Mag quadrupole uu3D[1]"] << uu3D[1] ;
    obs["Mag quadrupole uu3D[2]"] << uu3D[2] ;
    obs["Mag quadrupole uu3D"] << sqrt(uu3Dsq) ;
    obs["Mag quadrupole uu3D E"] << sqrt(uu3Dsq) * energy_ ;
    obs["Mag quadrupole uu3D^2"] << uu3Dsq ;
    obs["Mag quadrupole uu3D^4"] << uu3Dsq * uu3Dsq;
    obs["Mag quadrupole uu2D[0]"] << uu2D[0] ;
    obs["Mag quadrupole uu2D[1]"] << uu2D[1] ;
    obs["Mag quadrupole uu2D[2]"] << uu2D[2] ;
    obs["Mag quadrupole uu2D"] << sqrt(uu2Dsq) ;
    obs["Mag quadrupole uu2D E"] << sqrt(uu2Dsq) * energy_ ;
    obs["Mag quadrupole uu2D^2"] << uu2Dsq;
    obs["Mag quadrupole uu2D^4"] << uu2Dsq * uu2Dsq;
  }

  // for exchange Monte Carlo
  typedef double weight_parameter_type;
  void set_beta(double beta) { beta_ = beta; }
  weight_parameter_type weight_parameter() const { return energy_; }
  static double log_weight(weight_parameter_type gw, double beta) { return - beta * gw; }

  void save(alps::ODump& dp) const { dp << mcs_ << spins_ << energy_; }
  void load(alps::IDump& dp) { dp >> mcs_ >> spins_ >> energy_; }

protected:
  void update_energy() {
    bond_iterator itr, itr_end;
//
    double xx,yy,zz,xy,yx;
//
    energy_ = 0;
    for (boost::tie(itr, itr_end) = bonds(); itr != itr_end; ++itr){
//      energy_ -= coupling_[bond_type(*itr)] * (spins_[source(*itr)][2] * spins_[target(*itr)][2]);
//
      xx = spins_[source(*itr)][0] * spins_[target(*itr)][0];
      yy = spins_[source(*itr)][1] * spins_[target(*itr)][1];
      zz = spins_[source(*itr)][2] * spins_[target(*itr)][2];
      xy = spins_[source(*itr)][0] * spins_[target(*itr)][1];
      yx = spins_[source(*itr)][1] * spins_[target(*itr)][0];
      energy_ -=coupling_[0]*zz;
      energy_ -=coupling_[0]*coupling_delta*(xx+yy);
      energy_ -=coupling_[0]*coupling_q*(coupling_cos2phi[bond_type(*itr)]*(xx-yy)-coupling_sin2phi[bond_type(*itr)]*(xy+yx));
    }
//
    double magvec[3]={0.0, 0.0, 0.0};
    for (int s = 0; s < num_sites(); ++s) {
        magvec[0] += zvec[site_type(s)][0] * spins_[s][2];
        magvec[1] += zvec[site_type(s)][1] * spins_[s][2];
        magvec[2] += zvec[site_type(s)][2] * spins_[s][2];
    }
    energy_ -= (hh[0]*magvec[0] + hh[1]*magvec[1] + hh[2]*magvec[2]);
//
  }

private:
  // parameteters
  double beta_;
  std::vector<double> coupling_;
//
  double L_a_hex;
  double zvec[4][3];
  double hh[3];
  double coupling_delta;
  double coupling_q;
  double coupling_cos2phi[3];
  double coupling_sin2phi[3];
//
  // configuration (need checkpointing)
  alps::mc_steps mcs_;
  std::vector<spin_type> spins_;
  double energy_;
  // random number distribution  
  dist_type dist_;
};

class heisenberg_evaluator : public alps::parapack::simple_evaluator {
public:
  heisenberg_evaluator(alps::Parameters const&) {}
  virtual ~heisenberg_evaluator() {}
//
  void evaluate(alps::ObservableSet& obs) const {
    if (obs.has("Inverse Temperature") && obs.has("Number of Sites") &&
        obs.has("Energy") && obs.has("Energy^2")) {
      alps::RealObsevaluator beta = obs["Inverse Temperature"];
      alps::RealObsevaluator n = obs["Number of Sites"];
      alps::RealObsevaluator ene = obs["Energy"];
      alps::RealObsevaluator ene2 = obs["Energy^2"];
      if (beta.count() && n.count() && ene.count() && ene2.count()) {
        alps::RealObsevaluator c("Specific Heat");
        c = beta.mean() * beta.mean() * (ene2 - ene * ene) / n.mean();
        obs.addObservable(c);
      }
    }
//
    if (obs.has("Magnetization 111 abs") && obs.has("Magnetization 111^2")) {
      alps::RealObsevaluator beta = obs["Inverse Temperature"];
      alps::RealObsevaluator n = obs["Number of Sites"];
      alps::RealObsevaluator m1 = obs["Magnetization 111 abs"];
      alps::RealObsevaluator m2 = obs["Magnetization 111^2"];
      if (beta.count() && n.count() && m1.count() && m2.count()) {
        alps::RealObsevaluator chi111("Suscptibility 111");
        chi111 = beta.mean() * (m2 - m1 * m1) * n.mean();
        obs.addObservable(chi111);
      }
    }
//
    if (obs.has("Magnetization 111^2") && obs.has("Magnetization 111^4")) {
      alps::RealObsevaluator m2 = obs["Magnetization 111^2"];
      alps::RealObsevaluator m4 = obs["Magnetization 111^4"];
      if (m2.count() && m4.count()) {
        alps::RealObsevaluator binder("X Binder Ratio of Magnetization 111");
        binder = m4 / (m2 * m2) ;
        obs.addObservable(binder);
      }
    }
//
    if (obs.has("Mag quadrupole uu3D") && obs.has("Mag quadrupole uu3D^2")) {
      alps::RealObsevaluator beta = obs["Inverse Temperature"];
      alps::RealObsevaluator n = obs["Number of Sites"];
      alps::RealObsevaluator m1 = obs["Mag quadrupole uu3D"];
      alps::RealObsevaluator m2 = obs["Mag quadrupole uu3D^2"];
      if (beta.count() && n.count() && m1.count() && m2.count()) {
        alps::RealObsevaluator chi_uu3D("Suscptibility uu3D");
        chi_uu3D = beta.mean() * (m2 - m1 * m1) * n.mean();
        obs.addObservable(chi_uu3D);
      }
    }
//
    if (obs.has("Mag quadrupole uu3D") && obs.has("Mag quadrupole uu3D E") &&
        obs.has("Energy")) {
      alps::RealObsevaluator beta = obs["Inverse Temperature"];
      alps::RealObsevaluator m1 = obs["Mag quadrupole uu3D"];
      alps::RealObsevaluator m2 = obs["Mag quadrupole uu3D E"];
      alps::RealObsevaluator ene = obs["Energy"];
      if (beta.count() && m1.count() && m2.count() && ene.count()) {
        alps::RealObsevaluator uu3D_log_prime("Tslope log uu3D");
        uu3D_log_prime = beta.mean() * beta.mean() * (m2/m1 - ene) ;
        obs.addObservable(uu3D_log_prime);
      }
    }
//
    if (obs.has("Mag quadrupole uu3D^2") && obs.has("Mag quadrupole uu3D^4")) {
      alps::RealObsevaluator m2 = obs["Mag quadrupole uu3D^2"];
      alps::RealObsevaluator m4 = obs["Mag quadrupole uu3D^4"];
      if ( m2.count() && m4.count()) {
        alps::RealObsevaluator binder3D("X Binder Ratio of quadrupole_uu3D");
        binder3D = m4 / (m2 * m2);
        obs.addObservable(binder3D);
      }
    }
//
    if (obs.has("Mag quadrupole uu2D") && obs.has("Mag quadrupole uu2D^2")) {
      alps::RealObsevaluator beta = obs["Inverse Temperature"];
      alps::RealObsevaluator n = obs["Number of 2D Sites"];
      alps::RealObsevaluator m1 = obs["Mag quadrupole uu2D"];
      alps::RealObsevaluator m2 = obs["Mag quadrupole uu2D^2"];
      if (beta.count() && n.count() && m1.count() && m2.count()) {
        alps::RealObsevaluator chi_uu2D("Suscptibility uu2D");
        chi_uu2D = beta.mean() * (m2 - m1 * m1) * n.mean();
        obs.addObservable(chi_uu2D);
      }
    }
//
    if (obs.has("Mag quadrupole uu2D") && obs.has("Mag quadrupole uu2D E") &&
        obs.has("Energy")) {
      alps::RealObsevaluator beta = obs["Inverse Temperature"];
      alps::RealObsevaluator m1 = obs["Mag quadrupole uu2D"];
      alps::RealObsevaluator m2 = obs["Mag quadrupole uu2D E"];
      alps::RealObsevaluator ene = obs["Energy"];
      if (beta.count() && m1.count() && m2.count() && ene.count()) {
        alps::RealObsevaluator uu2D_log_prime("Tslope log uu2D");
        uu2D_log_prime = beta.mean() * beta.mean() * (m2/m1 - ene) ;
        obs.addObservable(uu2D_log_prime);
      }
    }
//
    if (obs.has("Mag quadrupole uu2D^2") && obs.has("Mag quadrupole uu2D^4")) {
      alps::RealObsevaluator m2 = obs["Mag quadrupole uu2D^2"];
      alps::RealObsevaluator m4 = obs["Mag quadrupole uu2D^4"];
      if ( m2.count() && m4.count()) {
        alps::RealObsevaluator binder2D("X Binder Ratio of quadrupole_uu2D");
        binder2D = m4 / (m2 * m2);
        obs.addObservable(binder2D);
      }
    }
//
  }
};

#endif // PARAPACK_EXAMPLE_HEISENBERG_HEISENBERG_H
