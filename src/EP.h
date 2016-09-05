#ifndef Physics_Analysis_EP_H
#define Physics_Analysis_EP_H 

#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"    //No NTuple!
//#include "VertexFit/ReadBeamParFromDb.h"
#include "PartPropSvc/PartPropSvc.h"  

#include <string>
#include <TTree.h>
#include <TROOT.h>
#include <TLorentzVector.h>
#include <TFile.h>

#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "RootCnvSvc/RootCnvSvc.h"

#include "TROOT.h"
#include "TBenchmark.h"

#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "VertexFit/WTrackParameter.h"
#include "MdcRecEvent/RecMdcTrack.h"
#include "MdcRecEvent/RecMdcKalTrack.h"

class EP : public Algorithm {

	public:
		EP(const std::string& name, ISvcLocator* pSvcLocator);
		StatusCode initialize();
		StatusCode execute();
		StatusCode finalize();  

	private:
		IPartPropSvc *p_PartPropSvc;
		HepPDT::ParticleDataTable* m_particleTable;

               void corgen(HepMatrix &, HepVector &, int );
               void corset(HepSymMatrix &, HepMatrix &, int );
               void calibration(RecMdcKalTrack * , HepVector &, int );
		//Declare r0, z0 and cos cut for charged tracks
		Double_t m_vr0cut;
		Double_t m_vz0cut;
		Double_t m_ccoscut;
		Double_t m_ptcut;
		Double_t m_pcut;
		Double_t emchit_p;
		Double_t emchit_m;
		//Declare energy, dphi, dthe cuts for fake gamma's
		Double_t m_isoAngleCut;
		
		//Declare flag for cut flow
		Int_t m_cutFlow;

		//Declare flag for analysis
		Int_t m_readSig;

		//Declare flag for MCTruth
		Int_t m_mcTruth;

		//Declare type of psi. For Jpsi psitype=1, for psip psitype=0
		Int_t m_psiType;
                
                Int_t m_useVxfitCut;
                Int_t m_useKmfitCut;
		//Declare name of output file
		std::string m_OutputFileName;
		TFile *saveFile;

		//Define TCutFlow here
		TTree *TCutFlow;
		Int_t cutflow;

		//Define TreeAna here
		TTree *TreeAna;

		//For header
		Int_t runid;
		Int_t evtid;
		Int_t nevt;
                Int_t n_gamma;
                Int_t n_endgamma;
                Int_t n_charged;
               
  		Double_t m_dang[500];
                Double_t m_Rxy[500];
                Double_t m_tdc[500];
                Double_t m_z0[500];
		//For storing 4-mom of particles. Change names below according to your channel
                //After 4c
                TLorentzVector *gamma1;
                TLorentzVector *gamma2;
                TLorentzVector *gamma3;
                TLorentzVector *gamma4;
                TLorentzVector *pip;
                TLorentzVector *pim;
                TLorentzVector *kp;
                TLorentzVector *km;

                //Before 4c
		TLorentzVector *gamma1_unfitted;
                TLorentzVector *gamma2_unfitted;
                TLorentzVector *gamma3_unfitted;
                TLorentzVector *gamma4_unfitted;
                TLorentzVector *pip_unfitted;
                TLorentzVector *pim_unfitted;
                TLorentzVector *kp_unfitted;
                TLorentzVector *km_unfitted;

                TLorentzVector *ep_unfitted;
                TLorentzVector *em_unfitted;
		//For further cut. Add your own cut here
 		TLorentzVector *etap_mc;
 		TLorentzVector *pi0_mc;
 		TLorentzVector *eta_mc;
 		TLorentzVector *gamma2_mc;
		TLorentzVector *ep_mc;
		TLorentzVector *em_mc;
  		TLorentzVector *gamma1_mc;
  		TLorentzVector *gamma_mc;
   		TLorentzVector *Rgamma_mc;
  		TLorentzVector *gamma3_mc;

        	Double_t V_xy;
		Double_t vxchisq;
		Double_t kmchisq_4c;
                Double_t Chisq_low;
		Double_t EOP_p;
		Double_t Edpop;
		Double_t Edpom;
                Double_t EOP_m;
		Double_t Prob_ep;
		Double_t Prob_em;
		Double_t Prob_pip;
                Double_t Prob_pim;
                Double_t P_ep;
                Double_t P_ep2;
                Double_t P_em;
                Double_t P_em2;
		Double_t PD_p;
                Double_t PD_m;
//                Double_t kmchisq_eta[500];
          //gamma conversion
	 	Double_t m_xconv1;
                Double_t m_yconv1;
	 	Double_t m_zconv1;
		Double_t m_rconv1;
                Double_t m_xconv2;
                Double_t m_yconv2;
                Double_t m_zconv2;
                Double_t m_rconv2;
                Double_t m_xiep;
                Double_t m_deltaxy;
               
   		Double_t m_deltaz1;
                Double_t m_deltaz2;

                Double_t m_lep;
                Double_t m_psipair;
//                Double_t m_dgamma;
                Double_t MEE;
                Double_t m_vx_x;
                Double_t m_vx_y;
                Double_t m_vx_r;
	        Double_t m_thetaeg1;
                Double_t m_thetaeg2; 
 		Double_t m_cthep;
                Double_t m_ptrkp;
                Double_t m_ptrkm;
                Double_t m_mgamma;
                Double_t m_egamma;
                Double_t m_theta;
                Double_t m_cosTheta;
                Double_t m_phi;
		Double_t m_rp;
                Double_t m_re;
                Double_t m_deltaeq;
      		Double_t m_case;
             		//Define TMCTruth Here

		TTree *TMCTruth;
//MC truth is stored in TreeAna.
		Int_t pdgid[500];
		Int_t indexmc;
		Int_t motheridx[500];

                TLorentzVector *p_f0_mctruth;
		Double_t mass;

		double m_epos_ratio;
		double m_eneg_ratio;
		
		double moment1_p;
		double moment2_p;
		double moment1_m;
		double moment2_m;
		double miss;

		double m_ppos_ratio;
		double m_depos_ratio;
		double m_eppos_ratio;
		double m_pneg_ratio;
		double m_deneg_ratio;
		double m_epneg_ratio;
		double m_tmass_ratio;
		double m_hitspos_ratio;
		double m_hitsneg_ratio;
		double m_hits;		

		double chiDeDx[4][5];
		double chiTof1[4][5];
		double chiTof2[4][5];
		double chisq_pid[4][5];
		double prob[4][5];

};

#endif 
