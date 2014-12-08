// -*- C++ -*-
#include "Rivet/Analyses/MC_JetAnalysis.hh"
#include "Rivet/Math/MathUtils.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/RivetAIDA.hh"
#include <map>
#include <sstream>

namespace Rivet {


	/// @brief MC validation analysis for higgs [-> tau tau] + jets events
	class MC_PHOTONPHOTON_JOEY : public Analysis {
		private:
			double _phy, _ph1pT, _ph2pT, _phisoE, _phisodR;
			double _jphdR, _jeta, _jR, _jpTtag, _jpT, _jjdy, _jjmass;
			std::vector<double> _j1pTs;
			std::map<std::string,AIDA::IHistogram1D *> histos, jetvetohistos;
		public:
			MC_PHOTONPHOTON_JOEY() : 
				Analysis("MC_PHOTONPHOTON_JOEY"),
				_phy(2.37), _ph1pT(43.75), _ph2pT(31.25), _phisoE(14.), _phisodR(0.4),
				_jphdR(0.4), _jeta(4.4), _jR(0.4), _jpT(30.)
		{}

			double sqr(const double& x) { return x*x; }

			void init() {
				FinalState fs;
				IdentifiedFinalState photons(-_phy,_phy,_ph2pT*GeV);
				photons.acceptId(PHOTON);
				addProjection(fs, "FS");
				addProjection(photons, "Photons");

				inithistos();
			}



			void inithistos() {

				histos["XS"] = bookHistogram1D("XS",1,0.,1.);
				histos["m_gammagamma"] = bookHistogram1D("m_gammagamma",BWspace(21,124.9,125.1,125.,0.00407));

				std::vector<double> H_pT_bins,H_y_bins,H_y_fine_bins,jet1_pT_bins,jet2_pT_bins,jet3_pT_bins, deltaphi_jj_bins,deltaphi_Hjj_bins,deltay_jj_bins,deltay_yy_bins,deltay_H_jj_bins;
				std::vector<double> Hjj_pT_bins,Hjj_pT_fine_bins,HT_jets_bins,HT_jets,H_pT_jj_fine_bins,H_pT_0j_bins;
				std::vector<double> H_pT_1j_bins,H_pT_2j_bins,jet1_y_bins,H_pT_j_bins,H_pT_jj_bins,Hj_pT_bins,jet1_y_fine_bins,deltaphi_jj_fine_bins,pTt_bins;
				std::vector<double> deltay_jj_fine_bins,dijet_mass_bins,dijet_mass_fine_bins,H_dijet_mass_bins,tau_jet_bins;

				histos["NJet_incl"] = bookHistogram1D("NJet_incl",4,-0.5,3.5);
				histos["NJet_excl"] = bookHistogram1D("NJet_excl",4,-0.5,3.5);

				histos["NJet_incl_50"] = bookHistogram1D("NJet_incl_50",4,-0.5,3.5);
				histos["NJet_excl_50"] = bookHistogram1D("NJet_excl_50",4,-0.5,3.5);

				H_pT_bins += 0.,20.,30.,40.,50.,60.,80.,100.,200.;
				H_pT_j_bins += 0.,50.,70.,100.,145.,200.;
				H_pT_jj_bins += 0.,90.,130.,170.,200;
				H_pT_jj_fine_bins += 0.,5.,10.,15.,20.,25.,30.,35.,40.,45.,50.,55.,60.,65.,70.,75.,80.,90.,100.,110.,120.,130.,140.,150,160.,170.,180.,190.,200.;

				Hj_pT_bins += 0.,18.,35.,60.,150.;
				Hjj_pT_bins += 0.,30.,55.,80.,140.,200;

				H_pT_0j_bins += 0.,20.,30.,45.,200.;
				H_pT_1j_bins += 0.,40.,60.,95.,200.;
				H_pT_2j_bins += 0.,90.,140.,200.;

				histos["H_pT"] = bookHistogram1D("H_pT",H_pT_bins);
				histos["H_pT_fine"] = bookHistogram1D("H_pT_fine",H_pT_jj_fine_bins);

				histos["H_j_pT"] = bookHistogram1D("H_j_pT",H_pT_j_bins);
				histos["H_j_pT_fine"]= bookHistogram1D("H_j_pT_fine",H_pT_jj_fine_bins);

				histos["H_jj_pT"] = bookHistogram1D("H_jj_pT",H_pT_jj_bins);
				histos["H_jj_pT_fine"] = bookHistogram1D("H_jj_pT_fine",H_pT_jj_fine_bins);

				histos["Hj_pT"] = bookHistogram1D("Hj_pT",Hj_pT_bins);
				histos["Hj_pT_fine"]= bookHistogram1D("Hj_pT_fine",H_pT_jj_fine_bins);

				histos["Hjj_pT"] = bookHistogram1D("Hjj_pT",Hjj_pT_bins);
				histos["Hjj_pT_fine"] = bookHistogram1D("Hjj_pT_fine",H_pT_jj_fine_bins);

				histos["H_pT_excl"] = bookHistogram1D("H_pT_excl",H_pT_0j_bins);
				histos["H_pT_fine_excl"]= bookHistogram1D("H_pT_fine_excl",H_pT_jj_fine_bins);

				histos["Hj_pT_excl"] = bookHistogram1D("Hj_pT_excl",H_pT_1j_bins);
				histos["Hj_pT_fine_excl"]= bookHistogram1D("Hj_pT_fine_excl",H_pT_jj_fine_bins);

				histos["Hjj_pT_excl"] = bookHistogram1D("Hjj_pT_excl",H_pT_2j_bins);
				histos["Hjj_pT_fine_excl"] = bookHistogram1D("Hjj_pT_fine_excl",H_pT_jj_fine_bins);

				histos["H_j_pT_excl"] = bookHistogram1D("H_j_pT_excl",H_pT_1j_bins);
				histos["H_j_pT_fine_excl"]= bookHistogram1D("H_j_pT_fine_excl",H_pT_jj_fine_bins);

				histos["H_jj_pT_excl"] = bookHistogram1D("H_jj_pT_excl",H_pT_2j_bins);
				histos["H_jj_pT_fine_excl"] = bookHistogram1D("H_jj_pT_fine_excl",H_pT_jj_fine_bins);

				jet1_pT_bins += 0.,30.,50.,70.,100.,140.,500.;
				jet2_pT_bins += 0,30,40,50,140,500.;
				jet3_pT_bins += 0,30,50,150,500.;

				histos["jet1_pT"] = bookHistogram1D("jet1_pT",jet1_pT_bins);
				histos["jet1_pT_fine"] = bookHistogram1D("jet1_pT_fine",H_pT_jj_fine_bins);
				histos["jet1_pT_excl"] = bookHistogram1D("jet1_pT_excl",jet1_pT_bins);
				histos["jet1_pT_excl_fine"] = bookHistogram1D("jet1_pT_excl_fine",H_pT_jj_fine_bins);

				histos["jet2_pT"] = bookHistogram1D("jet2_pT",jet2_pT_bins);
				histos["jet2_pT_fine"] = bookHistogram1D("jet2_pT_fine",H_pT_jj_fine_bins);

				histos["jet3_pT"] = bookHistogram1D("jet3_pT",jet3_pT_bins);
				histos["jet3_pT_fine"] = bookHistogram1D("jet3_pT_fine",H_pT_jj_fine_bins);


				H_y_bins += 0.,0.3,0.6,0.9,1.2,1.6,2.4;
				H_y_fine_bins += 0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.,2.1,2.2,2.3,2.4;

				histos["H_y"] = bookHistogram1D("H_y",H_y_bins);
				histos["H_y_fine"] = bookHistogram1D("H_y_fine",H_y_fine_bins);

				jet1_y_bins += 0.,1.,2.,3.,4.4,5;
				jet1_y_fine_bins += 0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.4;

				histos["jet1_y"] = bookHistogram1D("jet1_y",jet1_y_bins);
				histos["jet1_y_fine"] = bookHistogram1D("jet1_y_fine",jet1_y_fine_bins);

				histos["jet2_y"] = bookHistogram1D("jet2_y",jet1_y_bins);
				histos["jet2_y_fine"] = bookHistogram1D("jet2_y_fine",jet1_y_fine_bins);

				histos["jet3_y"] = bookHistogram1D("jet3_y",jet1_y_bins);
				histos["jet3_y_fine"] = bookHistogram1D("jet3_y_fine",jet1_y_fine_bins);

				histos["cos_theta_star"] = bookHistogram1D("cos_theta_star",10,0.,1.);
				histos["cos_theta_star_80"] = bookHistogram1D("cos_theta_star_80",4,0.,1.);
				histos["cos_theta_star_200"] = bookHistogram1D("cos_theta_star_200",4,0.,1.);
				histos["cos_theta_star_gt200"] = bookHistogram1D("cos_theta_star_gt200",4,0.,1.);

				histos["loose"] = bookHistogram1D("loose",3,0.,3.);
				histos["tight"] = bookHistogram1D("tight",3,0.,3.);

				deltaphi_jj_bins += 0.,PI/3.,2.*PI/3.,5.*PI/6.,PI;
				deltaphi_jj_fine_bins += 0.,0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2;

				histos["deltaphi_jj"] = bookHistogram1D("deltaphi_jj",deltaphi_jj_bins);
				histos["deltaphi_jj_fine"] = bookHistogram1D("deltaphi_jj_fine",deltaphi_jj_fine_bins);
				histos["deltaphi_jj_excl"] = bookHistogram1D("deltaphi_jj_excl",deltaphi_jj_bins);
				histos["deltaphi_jj_VBF"] = bookHistogram1D("deltaphi_jj_VBF",deltaphi_jj_fine_bins);

				deltay_jj_bins += 0.,2.,4.,5.5,8.8,10;
				deltay_jj_fine_bins += 0.,2.,4.,5.5,8.8;
				//0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5.,6.,6.5,7.,7.5,8.,8.5,9.;

				histos["deltay_jj"] = bookHistogram1D("deltay_jj",deltay_jj_bins);
				histos["deltay_jj_fine"] = bookHistogram1D("deltay_j_fine",deltay_jj_fine_bins);

				deltay_yy_bins +=0,0.3,0.6,0.9,1.2,1.5,2.0,2.55;
				histos["deltay_yy"] = bookHistogram1D("deltay_yy",deltay_yy_bins);

				dijet_mass_bins += 0,200.,400.,650.,1000.;
				dijet_mass_fine_bins += 0,50.,100.,150.,200.,250.,300.,350.,400.,450.,500.,550.,600.,650.,700.,750.,800.,850.,900.,950.,1000.;

				histos["dijet_mass"] = bookHistogram1D("dijet_mass",dijet_mass_bins);
				histos["dijet_mass_fine"] = bookHistogram1D("dijet_mass_fine",dijet_mass_fine_bins);

				H_dijet_mass_bins += 0,450.,700.,1100.,1200.;

				histos["H_dijet_mass"]= bookHistogram1D("H_dijet_mass",H_dijet_mass_bins);

				deltaphi_Hjj_bins += 0.,2.8,3,3.1,PI;

				histos["deltaphi_Hjj"] = bookHistogram1D("deltaphi_Hjj",deltaphi_Hjj_bins);
				histos["deltaphi_Hjj_excl"] = bookHistogram1D("deltaphi_Hjj_excl",deltaphi_Hjj_bins);

				deltay_H_jj_bins +=0,0.4,0.7,1.3,4.5;

				histos["deltay_H_jj"] = bookHistogram1D("deltay_H_jj",deltay_H_jj_bins);
				histos["deltay_H_jj_fine"] = bookHistogram1D("deltay_H_jj_fine",deltay_jj_fine_bins);

				HT_jets_bins +=0,30,50,70,150,250,500;
				//      HT_jets_bins +=0,100,200,300,400,500,600,700,800,900,1000;
				histos["HT_jets_hist"] = bookHistogram1D("HT_jets_hist",HT_jets_bins);

				tau_jet_bins +=0,8,16,30,85;
				histos["tau_jet1"] = bookHistogram1D("tau_jet1",tau_jet_bins);
				histos["tau_jet2"] = bookHistogram1D("tau_jet2",tau_jet_bins);
				histos["tau_jet3"] = bookHistogram1D("tau_jet3",tau_jet_bins);
				histos["tau_jet_max"]  = bookHistogram1D("tau_jet_max",tau_jet_bins);
				histos["sum_tau_jet"]  = bookHistogram1D("sum_tau_jet",tau_jet_bins);


				pTt_bins += 0,10,20,30,40,60,80,150,500;
				histos["pTt"] = bookHistogram1D("pTt",pTt_bins);

			}

			/// Do the analysis
			void analyze(const Event & e) {
				const double weight = e.weight();
				ParticleVector photons =
					applyProjection<IdentifiedFinalState>(e, "Photons").particlesByPt();
				// require at least two photons
				if (photons.size() < 2) vetoEvent;
				// photon isolation
				ParticleVector fs = applyProjection<FinalState>(e, "FS").particles();
				ParticleVector jetparts;
				std::vector<FourMomentum> phs;
				// add all particles not in the photons to jetparts
				foreach (const Particle& p, fs) {
					bool addit(true);
					if (p.pdgId()==PHOTON) {
						foreach (const Particle& ph, photons) {
							if (p.momentum()==ph.momentum()) {
								addit=false;
								break;
							}
						}
					}
					if (addit) jetparts.push_back(p);
				}
				// energy in cone
				// add to jet-input everything but the isolated photons
				foreach (const Particle& ph, photons) {
					FourMomentum mom_in_Cone = -ph.momentum();
					foreach (const Particle& p, fs) {
						if (deltaR(ph,p) < _phisodR) mom_in_Cone += p.momentum();
					}
					if (mom_in_Cone.Et()<_phisoE) phs.push_back(ph.momentum());
					else jetparts.push_back(ph);
				}

				// at least two isolated photons
				// TODO: use subleading photons if leading ones are not isolated?
				if (phs.size() < 2) vetoEvent;
				// photon pT cuts
				if (phs[0].pT() < _ph1pT) vetoEvent;
				if (phs[1].pT() < _ph2pT) vetoEvent;
				// isolate photons from one-another
				if (deltaR(phs[0],phs[1]) < _phisodR) vetoEvent;
				// TODO: feed all but the two isolated photons back into jets?

				// define diphoton (Higgs) momentum
				FourMomentum hmom(phs[0]+phs[1]);

				// calculate jets
				FastJets jetpro(FastJets::ANTIKT, _jR);
				jetpro.calc(jetparts);
				Jets jets;
				// apply etamax and dRmin to photons
				foreach (const Jet& jetcand, jetpro.jetsByPt(_jpT*GeV)) {
					FourMomentum jmom = jetcand.momentum();
					if (fabs(jetcand.momentum().rapidity()) < _jeta && 
							deltaR(phs[0],jetcand.momentum()) > _jphdR &&
							deltaR(phs[1],jetcand.momentum()) > _jphdR) {
						jets.push_back(jetcand);
					}
				}

				// fill histograms

				histos["XS"]->fill(0.5,weight);
				histos["m_gammagamma"]->fill(hmom.mass()/GeV,weight);

				histos["H_pT"]->fill(hmom.pT()/GeV,weight);
				histos["H_pT_fine"]->fill(hmom.pT()/GeV,weight);

				histos["H_y"]->fill(fabs(hmom.rapidity()),weight);
				histos["H_y_fine"]->fill(fabs(hmom.rapidity()),weight);

				// |costheta*| from 1307.1432
				double cts(abs(sinh(phs[0].eta()-phs[1].eta()))/
						sqrt(1.+sqr(hmom.pT()/hmom.mass()))
						* 2.*phs[0].pT()*phs[1].pT()/sqr(hmom.mass()));
				histos["cos_theta_star"]->fill(cts,weight);

				if(hmom.pT()<80){
					histos["cos_theta_star_80"]->fill(cts,weight);
				}

				if(hmom.pT()>80 && hmom.pT()<200){
					histos["cos_theta_star_200"]->fill(cts,weight);
				}

				if(hmom.pT()>200){
					histos["cos_theta_star_gt200"]->fill(cts,weight);
				}

				histos["NJet_excl"]->fill(jets.size(),weight);
				for (size_t i(0);i<4;++i) {
					if (jets.size()>=i) histos["NJet_incl"]->fill(i,weight);
				}

				double pTt = fabs(phs[0].px()*phs[1].py()-phs[1].px()*phs[0].py())/((phs[0]-phs[1]).pT()*2);
				histos["pTt"]->fill(pTt,weight);

				double deltay_yy = fabs(phs[0].rapidity() - phs[1].rapidity());
				histos["deltay_yy"]->fill(deltay_yy,weight);

				// njets == 0;
				if (jets.size()==0){
					histos["H_pT_excl"]->fill(hmom.pT()/GeV,weight);
					histos["H_pT_fine_excl"]->fill(hmom.pT()/GeV,weight);
					// 6/2 added fill for jet1_pT for 0-30 GeV bin, i.e. no jets
					histos["jet1_pT"]->fill(10,weight);
					// 6/2 added fill for overflow bin for 0 jets in event
					histos["deltay_jj"]->fill(9,weight);
					// 6/2 added fill for overflow bin for 0 jets in event
					histos["Hjj_pT"]->fill(160,weight);
					//6/2 added fill for overflow bin for 0 jets in event
					histos["jet2_y"]->fill(4.6,weight);
					//6/2 added fill for overflow bin for 0 jets in event
					histos["jet2_pT"]->fill(400,weight);
				}

				// njets > 0;
				if (jets.size()>0) {
					const FourMomentum& j1(jets[0].momentum());

					histos["jet1_pT"]->fill(j1.pT()/GeV,weight);
					histos["jet1_pT_fine"]->fill(j1.pT()/GeV,weight);
					histos["jet1_y"]->fill(j1.rapidity(),weight);
					histos["jet1_y_fine"]->fill(j1.rapidity(),weight);
					histos["Hj_pT"]->fill((hmom+j1).pT()/GeV,weight);
					histos["Hj_pT_fine"]->fill((hmom+j1).pT()/GeV,weight);
					histos["H_j_pT"]->fill(hmom.pT()/GeV,weight);
					histos["H_j_pT_fine"]->fill(hmom.pT()/GeV,weight);

					double tauJet1 = sqrt( pow(j1.pT(),2) + pow(j1.mass(),2))/(2.0*cosh(j1.rapidity() - (phs[0]+phs[1]).rapidity()));
					histos["tau_jet1"]->fill(tauJet1,weight);

					// njets == 1;
					if (jets.size()==1){
						histos["Hj_pT_excl"]->fill((hmom+j1).pT()/GeV,weight);
						histos["Hj_pT_fine_excl"]->fill((hmom+j1).pT()/GeV,weight);
						histos["H_j_pT_excl"]->fill(hmom.pT()/GeV,weight);
						histos["H_j_pT_fine_excl"]->fill(hmom.pT()/GeV,weight);
						histos["jet1_pT_excl"]->fill(j1.pT()/GeV,weight);
						histos["jet1_pT_excl_fine"]->fill(j1.pT()/GeV,weight);
						// 6/2 added fill for j2_pT for 0-30 GeV bins, i.e. no 2nd jet
						histos["jet2_pT"]->fill(10,weight);
						// 6/2 added fill for overflow bin for 1 jet in event
						histos["deltay_jj"]->fill(9,weight);
						// 6/2 added fill for overflow bin for 1 jet in event
						histos["Hjj_pT"]->fill(160,weight);
						//6/2 added fill for overflow bin for 1 jet in event
						histos["jet2_y"]->fill(4.6,weight);
					}
				}

				// njets > 1;
				if (jets.size()>1) {
					const FourMomentum& j1(jets[0].momentum());
					const FourMomentum& j2(jets[1].momentum());

					histos["deltaphi_jj"]->fill(deltaPhi(j1,j2),weight);
					histos["deltaphi_jj_fine"]->fill(deltaPhi(j1,j2),weight);
					histos["deltaphi_Hjj"]->fill(deltaPhi(hmom,j1+j2),weight);
					histos["Hjj_pT"]->fill((hmom+j1+j2).pT()/GeV,weight);
					histos["Hjj_pT_fine"]->fill((hmom+j1+j2).pT()/GeV,weight);
					histos["H_jj_pT"]->fill(hmom.pT()/GeV,weight);
					histos["H_jj_pT_fine"]->fill(hmom.pT()/GeV,weight);
					histos["jet2_pT"]->fill(j2.pT()/GeV,weight);
					histos["jet2_pT_fine"]->fill(j2.pT()/GeV,weight);
					histos["jet2_y"]->fill(fabs(j2.rapidity()),weight);
					histos["jet2_y_fine"]->fill(fabs(j2.rapidity()),weight);
					histos["dijet_mass"]->fill((j1+j2).mass(),weight);
					histos["dijet_mass_fine"]->fill((j1+j2).mass(),weight);
					histos["H_dijet_mass"]->fill((hmom+j1+j2).mass(),weight);
					histos["deltay_jj"]->fill(fabs(j1.rapidity()-j2.rapidity()),weight);
					histos["deltay_jj_fine"]->fill(fabs(j1.rapidity()-j2.rapidity()),weight);
					histos["deltay_H_jj"]->fill(fabs((hmom.rapidity()-(j1+j2).rapidity())),weight);
					histos["deltay_H_jj_fine"]->fill(fabs((hmom.rapidity()-(j1+j2).rapidity())),weight);

					if(fabs(j1.rapidity()-j2.rapidity())>2.8){
						if((j1+j2).mass()>400){
							// 6/13 realized that delta-phi cut should be between Higgs and dijet system
							histos["deltaphi_jj_VBF"]->fill(deltaPhi(hmom,j1+j2),weight);
							// 6/2 added loose and tight histograms to tally cross section as cross-checks
							histos["loose"]->fill(1,weight);
							if(fabs(deltaPhi(hmom,j1+j2)>2.6)){
								histos["tight"]->fill(1,weight);
							}
						}
					}

					double tauJet2 = sqrt( pow(j2.pT(),2) + pow(j2.mass(),2))/(2.0*cosh(j2.rapidity() - (phs[0]+phs[1]).rapidity()));
					histos["tau_jet2"]->fill(tauJet2,weight);

					// njets == 2;
					if (jets.size()==2) {
						histos["deltaphi_jj_excl"]->fill(deltaPhi(j1,j2),weight);
						histos["deltaphi_Hjj_excl"]->fill(deltaPhi(hmom,j1+j2),weight);
						histos["Hjj_pT_excl"]->fill((hmom+j1+j2).pT()/GeV,weight);
						histos["Hjj_pT_fine_excl"]->fill((hmom+j1+j2).pT()/GeV,weight);
						histos["H_jj_pT_excl"]->fill(hmom.pT()/GeV,weight);
						histos["H_jj_pT_fine_excl"]->fill(hmom.pT()/GeV,weight);
						// 6/2 added fill for j3_pT 0-30 GeV, i.e. no jet3
						histos["jet3_pT"]->fill(10,weight);
					}
				}

				// njets > 2;
				if (jets.size()>2) {
					const FourMomentum& j3(jets[2].momentum());

					histos["jet3_pT"]->fill(j3.pT()/GeV,weight);
					histos["jet3_pT_fine"]->fill(j3.pT()/GeV,weight);
					histos["jet3_y"]->fill(j3.rapidity(),weight);
					histos["jet3_y_fine"]->fill(j3.rapidity(),weight);

					double tauJet3 = sqrt( pow(j3.pT(),2) + pow(j3.mass(),2))/(2.0*cosh(j3.rapidity() - (phs[0]+phs[1]).rapidity()));
					histos["tau_jet3"]->fill(tauJet3,weight);
				}

				double HT_jets = 0; 

				for(size_t i(0);i<jets.size();i++){
					const FourMomentum& jetlocal(jets[i].momentum());
					HT_jets += jetlocal.pT();
				}

				histos["HT_jets_hist"]->fill(HT_jets,weight);

				double tau_jet_cut=8;
				double max_tj=0;
				double sum_tj=0;
				double tauJet;
				for(size_t i(0);i<jets.size();i++) {
					tauJet = sqrt( pow(jets[i].momentum().pT(),2) + pow(jets[i].momentum().mass(),2))/(2.0*cosh(jets[i].momentum().rapidity() - (phs[0]+phs[1]).rapidity()));
					if ( tauJet > tau_jet_cut ) {
						if (tauJet > max_tj)	
							max_tj=tauJet;
					}
				}
				histos["tau_jet_max"]->fill(max_tj,weight);

				for (size_t i=0;i<jets.size();++i) {
					tauJet = sqrt( pow(jets[i].momentum().pT(),2) + pow(jets[i].momentum().mass(),2))/(2.0*cosh(jets[i].momentum().rapidity() - (phs[0]+phs[1]).rapidity()));
					if ( tauJet > tau_jet_cut )
						sum_tj+=tauJet;
				}

				histos["sum_tau_jet"]->fill(max_tj,weight);


				if (jets.size()>=0) {
					if(jets.size()>=1) {
						const FourMomentum& j1(jets[0].momentum());
						if(j1.pT()>50.*GeV) {
							histos["NJet_incl_50"]->fill(1,weight);
						}
						else {
							histos["NJet_incl_50"]->fill(0,weight);
							histos["NJet_excl_50"]->fill(0,weight);
						}
					}
					if(jets.size()==1) {
						const FourMomentum& j1(jets[0].momentum());
						if(j1.pT()>50.*GeV) {
							histos["NJet_excl_50"]->fill(1,weight);
						}
					}
					if(jets.size()>=2) {
						const FourMomentum& j2(jets[1].momentum());
						if(j2.pT()>50.*GeV) {
							histos["NJet_incl_50"]->fill(2,weight);
						}
					}
					if(jets.size()==2) {
						const FourMomentum& j2(jets[1].momentum());
						if(j2.pT()>50.*GeV) {
							histos["NJet_excl_50"]->fill(2,weight);
						}
					}
					if(jets.size()>=3) {
						const FourMomentum& j3(jets[2].momentum());
						if(j3.pT()>50.*GeV) {
							histos["NJet_incl_50"]->fill(3,weight);
						}
					}
					if(jets.size()==3) {
						const FourMomentum& j3(jets[2].momentum());
						if(j3.pT()>50.*GeV) {
							histos["NJet_excl_50"]->fill(3,weight);
						}
					}
				}
				if (jets.size()==0) {
					histos["NJet_incl_50"]->fill(0,weight);
					histos["NJet_excl_50"]->fill(0,weight);
				}
			}

			/// Finalize
			void finalize() {
				double scalefactor(crossSection()/sumOfWeights());
				for (std::map<std::string,AIDA::IHistogram1D *>::iterator 
						hit=histos.begin(); hit!=histos.end();hit++) 
					scale(hit->second,scalefactor);
			}
	};



	// The hook for the plugin system
	DECLARE_RIVET_PLUGIN(MC_PHOTONPHOTON_JOEY);
}
