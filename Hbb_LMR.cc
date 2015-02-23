#include <TH1F.h>
#include <string>
#include <sstream>
#include <cmath>
#include <TH2F.h>
#include <TH3F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TChain.h>
#include <TLorentzVector.h>
#include <TLegend.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <vector>
double pi=3.14159265;
double Jet_pt_cut=30.;
double Jet_eta_cut=2.5;
double H_mass=160.0;
double mH_diff_cut=50.;
double mH_mean_cut=15.;


int withinRegion(double mH1, double mH2, double r1=15., double r2=30., double mH1_c=H_mass, double mH2_c=H_mass)
{
	double r=pow(pow(mH1-mH1_c, 2)+pow(mH2-mH2_c, 2), 0.5);
	double angle=atan2(mH2-mH2_c, mH1-mH1_c);
	int ret=-1;
	if (r<r1) ret=0;
	else if (r>r1 && r<r2)
	{
		if (angle>=0 && angle<pi/2.) ret=1;
		else if (angle>=pi/2. && angle<pi) ret=4;
		else if (angle<0 && angle>=-pi/2.) ret=2;
		else if (angle<pi/2.&& angle>=-pi) ret=3;
		else std::cout<<"This is within annulus but not within any CR!"<<std::endl;
	}
	else ret=5;
	return ret;
}

typedef std::map<double, int> JetList;
int Hbb_LMR(std::string dir, std::string sample)
{

	std::string inputfilename=dir+sample+".root";
	TFile * f = new TFile(inputfilename.c_str());
	f->cd();	
	TTree *tree=(TTree*)f->Get("tree");
	std::cout<<"Opened input file "<<inputfilename<<std::endl;
	// Book variables
	double Vtype_;
	int evt, run;
	int nJet;               
	double Jet_id[5],Jet_puId[5],Jet_btagCSV[5],Jet_btagCMVA[5],Jet_rawPt[5],Jet_mcPt[5],Jet_mcFlavour[5],Jet_mcMatchId[5],Jet_pt[5],Jet_eta[5],Jet_phi[5],Jet_mass[5],Jet_hadronFlavour[5],Jet_btagProb[5],Jet_btagBProb[5],Jet_btagnew[5],Jet_btagCSVV0[5],Jet_chHEF[5],Jet_neHEF[5],Jet_chEmEF[5],Jet_neEmEF[5],Jet_chMult[5],Jet_leadTrackPt[5],Jet_mcEta[5],Jet_mcPhi[5],Jet_mcM[5], ht;

	//vector<float> Jet_id(0),Jet_puId(0),Jet_btagCSV(0),Jet_btagCMVA(0),Jet_rawPt(0),Jet_mcPt(0),Jet_mcFlavour(0),Jet_mcMatchId(0),Jet_pt(0),Jet_eta(0),Jet_phi(0),Jet_mass(0),Jet_hadronFlavour(0),Jet_btagProb(0),Jet_btagBProb(0),Jet_btagnew(0),Jet_btagCSVV0(0),Jet_chHEF(0),Jet_neHEF(0),Jet_chEmEF(0),Jet_neEmEF(0),Jet_chMult(0),Jet_leadTrackPt(0),Jet_mcEta(0),Jet_mcPhi(0),Jet_mcM(0);




	// GenParticleInfo genX, genH1, genH2,genH1B, genH1Bbar, genH2B, genH2Bbar;

	tree->SetBranchAddress("Vtype", &(Vtype_));
	tree->SetBranchAddress("evt",&evt); 
	tree->SetBranchAddress("run",&run);
	tree->SetBranchAddress("nJet",&nJet);              
	tree->SetBranchAddress("Jet_id",&Jet_id);            
	tree->SetBranchAddress("Jet_puId",&Jet_puId);          
	tree->SetBranchAddress("Jet_btagCSV",&Jet_btagCSV);       
	tree->SetBranchAddress("Jet_btagCMVA",&Jet_btagCMVA);      
	tree->SetBranchAddress("Jet_rawPt",&Jet_rawPt);         
	tree->SetBranchAddress("Jet_mcPt",&Jet_mcPt);          
	tree->SetBranchAddress("Jet_mcFlavour",&Jet_mcFlavour);     
	tree->SetBranchAddress("Jet_mcMatchId",&Jet_mcMatchId);     
	tree->SetBranchAddress("Jet_pt",&(Jet_pt));            
	tree->SetBranchAddress("Jet_eta",&Jet_eta);           
	tree->SetBranchAddress("Jet_phi",&Jet_phi);           
	tree->SetBranchAddress("Jet_mass",&Jet_mass);          
	tree->SetBranchAddress("Jet_hadronFlavour",&Jet_hadronFlavour); 
	tree->SetBranchAddress("Jet_btagProb",&Jet_btagProb);      
	tree->SetBranchAddress("Jet_btagBProb",&Jet_btagBProb);     
	tree->SetBranchAddress("Jet_btagnew",&Jet_btagnew);       
	tree->SetBranchAddress("Jet_btagCSVV0",&Jet_btagCSVV0);     
	tree->SetBranchAddress("Jet_chHEF",&Jet_chHEF);         
	tree->SetBranchAddress("Jet_neHEF",&Jet_neHEF);         
	tree->SetBranchAddress("Jet_chEmEF",&Jet_chEmEF);        
	tree->SetBranchAddress("Jet_neEmEF",&Jet_neEmEF);        
	tree->SetBranchAddress("Jet_chMult",&Jet_chMult);        
	tree->SetBranchAddress("Jet_leadTrackPt",&Jet_leadTrackPt);   
	tree->SetBranchAddress("Jet_mcEta",&Jet_mcEta);         
	tree->SetBranchAddress("Jet_mcPhi",&Jet_mcPhi);         
	tree->SetBranchAddress("Jet_mcM",&Jet_mcM);           

	//std::cout<<Jet_pt[0]<<"  "<<Jet_pt[1]<<"  "<<Jet_pt[2]<< "   "<<std::endl;

	TH1F *h_nJets=new TH1F("h_nJets", "# Central Jets with p_{T}> 30 GeV; n", 10, 0., 10.);
	TH1F *h_n8Jets=new TH1F("h_n8Jets", "# Central Jets with p_{T}> 90 GeV;  n", 10, 0., 10.); 
	TH1F *h_nCand=new TH1F("h_nCand", "# HH candidates", 20, 0., 20.);
	TH1F *h_nCand_true=new TH1F("h_nCand_true", "# HH candidates if true", 20, 0., 20.);	
	TH1F *h_nPV=new TH1F("h_nPV", "# of Primary Vertices; nPV", 51, 0., 50.); h_nPV->Sumw2();
	TH1F *h_nPV_weighted=new TH1F("h_nPV_weighted", "# of Primary Vertices after Reweighting; nPV", 51, 0., 50.); h_nPV_weighted->Sumw2();

	TH1F *h_JetpT1=new TH1F("h_JetpT1", "JetpT1", 50, 0., 900.);
	TH1F *h_JetpT2=new TH1F("h_JetpT2", "JetpT2", 50, 0., 500.);
	TH1F *h_JetpT3=new TH1F("h_JetpT3", "JetpT3", 50, 0., 350.);
	TH1F *h_JetpT4=new TH1F("h_JetpT4", "JetpT4", 50, 0., 250.);

	TH1F *h_CSV1=new TH1F("h_CSV1", "JetCSV1", 50, 0., 1.);
	TH1F *h_CSV2=new TH1F("h_CSV2", "JetCSV2", 50, 0., 1.);
	TH1F *h_CSV3=new TH1F("h_CSV3", "JetCSV3", 50, 0.,1.);
	TH1F *h_CSV4=new TH1F("h_CSV4", "JetCSV4", 50, 0., 1.);
        TH1F *h_CMVA1=new TH1F("h_CMVA1", "JetCMVA1", 50, 0., 1.);
        TH1F *h_CMVA2=new TH1F("h_CMVA2", "JetCMVA2", 50, 0., 1.);
        TH1F *h_CMVA3=new TH1F("h_CMVA3", "JetCMVA3", 50, 0.,1.);
        TH1F *h_CMVA4=new TH1F("h_CMVA4", "JetCMVA4", 50, 0., 1.);

	TH1F *h_CSV_T=new TH1F("h_CSV_T", "JetCSV Tight", 50, -1., 1.);


	TH1F *h_MET=new TH1F("h_MET", "MET; MET (GeV)", 100, 0, 200.);
	TH1F *h_MET_sig=new TH1F("h_MET_sig", "MET Significance; Sig", 20, 0., 20.);
	TH1F *h_HT=new TH1F("h_HT", "HT Distribution", 100, 0., 3000.);

	TH1F *h_H1_mass=new TH1F("h_H1_mass", "H mass; mass (GeV)", 100, 50., 180.);
	TH1F *h_H1_pT=new TH1F("h_H1_pT", "H1 p_{T}; p_{T} (GeV/c)", 50, 0., 900.);
	TH1F *h_H2_mass=new TH1F("h_H2_mass", "H2 mass; mass (GeV)", 100, 50., 180.);
	TH1F *h_H2_pT=new TH1F("h_H2_pT", "H2 p_{T}; p_{T} (GeV/c)", 50, 0., 900.);
	TH2F *h_H1_H2_mass = new TH2F("h_H1_H2_mass", " all comb if min(#sigma (#delta m)^{2}) ", 100, 30., 200., 100, 30., 200.);
	TH1F *h_HH_mass_diff_cand=new TH1F("h_HH_mass_diff_cand", "|#Deltam| between Higgs masses - all candidates", 50, 0., 200.);
	TH1F *h_HH_massNorm_diff=new TH1F("h_HH_massNorm_diff", "|#Deltam| between Higgs masses", 50, 0., 2.);
	TH1F *h_HH_CSV = new TH1F("h_HH_CSV"," Sum CSV | between the two higgs ", 70, -4., 4.); 
	TH1F *h_genX_mass=new TH1F("h_genX_mass", "Generator X mass; m_{X} GeV", 100, 0., 1000.);
	TH1F *h_genH1_mass=new TH1F("h_genH1_mass", "Generator H1 mass; m_{X} GeV", 100, 0., 1000.);
	TH1F *h_genH2_mass=new TH1F("h_genH2_mass", "Generator H2 mass; m_{X} GeV", 100, 0., 1000.);
	TH1F *h_Xmass = new TH1F("h_Xmass"," h_Xmass" , 100, 0., 1000.); 
	TH1F *h_Xmass_right = new TH1F("h_Xmass_right"," h_Xmass right" , 100, 0., 1000.); 
	TH1F *h_Xmass_wrong = new TH1F("h_Xmass_wrong"," h_Xmass wrong" , 100, 0., 1000.);
	TH1F *h_Xmass_first = new TH1F("h_Xmass_first"," h_Xmass first" , 100, 0., 1000.);
	TH1F *h_Xmass_other = new TH1F("h_Xmass_other"," h_Xmass other" , 100, 0., 1000.);		
	TH2F *h_Xmass2 = new TH2F("h_Xmass2"," h_Xmass vs njets " , 10, 4, 14., 100, 0., 1000.); 	
	TH1F *h_Xmass_comb = new TH1F("h_Xmass_comb"," h_Xmass SR" , 200, 0., 2000.);
	TH1F *h_Cuts=new TH1F("h_Cuts", "Cut flow", 16, 0, 16);
	TH1F *h_angle_bbar = new TH1F("h_angle_bbar","h_angle_bbar", 100, 0., 5.);
	TH1F *h_angle_hh = new TH1F("h_angle_hh","h_angle_hh", 100, -1.2, 1.2);


	std::string outfilename=sample+"_selected"+".root";
	TFile *outfile=new TFile(outfilename.c_str(), "recreate");

	TTree *outtree=tree->CloneTree(0);
	int H1jet1_i, H1jet2_i;
	int H2jet1_i, H2jet2_i;
	int    nCJets=0, n8Jets=0, nBJets=0, nBMJets=0, nAJets=0;

	outtree->Branch("H1jet1_i", &H1jet1_i, "H1jet1_i/I");
	outtree->Branch("H1jet2_i", &H1jet2_i, "H1jet2_i/I");
	outtree->Branch("H2jet1_i", &H2jet1_i, "H2jet1_i/I");
	outtree->Branch("H2jet2_i", &H2jet2_i, "H2jet2_i/I");
	outtree->Branch("Ht", &ht, "ht/F");
	// Loop over events
	int nEvents=tree->GetEntries();
	double nCut0=0, nCut1=0, nCut2=0, nCut3=0, nCut4=0, nCut5=0, nCut6=0, nCut7=0;
	int nCandidate =0;
	for (int i=0; i<nEvents; ++i)
	{
		++nCut0;
		tree->GetEvent(i);
		//std::cout<<Vtype_<<"  " <<Jet_pt[0]<<"  "<<Jet_pt[1]<<"  "<<Jet_pt[2]<< "   "<<std::endl;

		TLorentzVector ajet1_p4, ajet2_p4, ajet3_p4, ajet4_p4, p4;
		JetList jetList_CSV, jetList_pT; 
		nCJets=0;
		n8Jets=0;	
		nBJets=0;	
		nBMJets=0;	
		nAJets=0;
		ht=0; 
		for (int j=0; j<nJet; ++j)
		{
			if (fabs(Jet_eta[j])<2.5) 
			{
				jetList_pT[Jet_pt[j]]=j;
				if (Jet_pt[j]>Jet_pt_cut ) 
					{
					++nAJets;  
					jetList_CSV[Jet_btagCMVA[j]]=j;
					++nCJets;
					ht+=Jet_pt[j];
				}

			}
		}

		// Analysis begins here
		//if (triggerFlags[54])//||triggerFlags[57])
		{
			//std::cout<<Vtype<<std::endl;
			//if(triggerFlags[54]) ++nCut1;
			if (Vtype_==-1)
			{
				++nCut2;
				double weightPU=1.;
				h_nJets->Fill(nCJets, weightPU);

				if(nCJets>3){
					int c =0;  
					for (JetList::reverse_iterator iJet=jetList_CSV.rbegin(); iJet!=jetList_CSV.rend(); ++iJet)
					{
						if(Jet_pt[iJet->second]>90.) ++n8Jets;
						if(iJet->first > 0.71) ++nBMJets;

						++c;

					} //end for jet iterator
					h_CMVA1->Fill(Jet_btagCMVA[0]);
					h_CMVA2->Fill(Jet_btagCMVA[1]);
					h_CMVA3->Fill(Jet_btagCMVA[2]);
					h_CMVA4->Fill(Jet_btagCMVA[3]);
					h_CSV1->Fill(Jet_btagCSV[0]);
                                        h_CSV2->Fill(Jet_btagCSV[1]);
                                        h_CSV3->Fill(Jet_btagCSV[2]);
                                        h_CSV4->Fill(Jet_btagCSV[3]);
					h_JetpT1->Fill(Jet_pt[0]);
					h_JetpT2->Fill(Jet_pt[1]);
					h_JetpT3->Fill(Jet_pt[2]);
					h_JetpT4->Fill(Jet_pt[3]);	
					h_HT->Fill(ht);
					h_n8Jets->Fill(n8Jets, weightPU);
					++nCut3;
					bool foundHH=false; 
					bool match=false;
					int cont=0; double deltaR1; double deltaR2;
					TLorentzVector jet1_p4, jet2_p4, jet3_p4, jet4_p4, jet5;
					TLorentzVector H1_p4, H2_p4;
					double min_diff=50;
					for (int j=0; j<nCJets; ++j)
					{
						jet1_p4.SetPtEtaPhiM(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_mass[j]);  
						for (int k=0; k<nCJets; ++k)
						{
							if (k!=j)
							{
								jet2_p4.SetPtEtaPhiM(Jet_pt[k], Jet_eta[k], Jet_phi[k], Jet_mass[k]);
								cont++;        
								{

									for (int l=0; l<nCJets; ++l)	
									{     
										if(l!=j && l!=k )
										{
											jet3_p4.SetPtEtaPhiM(Jet_pt[l], Jet_eta[l], Jet_phi[l], Jet_mass[l]);
											for (int m=0; m<nCJets; ++m)
											{
												if (m!=l && m!=j && m!=k )
												{
													jet4_p4.SetPtEtaPhiM(Jet_pt[m], Jet_eta[m], Jet_phi[m], Jet_mass[m]);
													TLorentzVector diJet2_p4, diJet1_p4;	

													diJet2_p4=jet3_p4+jet4_p4;
													diJet1_p4=jet1_p4+jet2_p4;
													double diff=fabs(diJet1_p4.M()-diJet2_p4.M());     

													deltaR1=jet3_p4.DeltaR(jet4_p4);
													deltaR2=jet1_p4.DeltaR(jet2_p4);

 if (diff<min_diff && deltaR1<1.5 && deltaR2<1.5) {
{
													//if((diff<min_diff)&&(((diJet2_p4.M()<160 && diJet2_p4.M()> 90.) && (diJet1_p4.M()<160 && diJet1_p4.M()> 90.)))){
													//	if(((jet2_p4.Pt()>90&&jet4_p4.Pt()>90) || (jet3_p4.Pt()>90&&jet1_p4.Pt()>90)||(jet1_p4.Pt()>90&&jet4_p4.Pt()>90) || (jet3_p4.Pt()>90&&jet2_p4.Pt()>90)|| (jet1_p4.Pt()>90&&jet2_p4.Pt()>90) || (jet3_p4.Pt()>90&&jet4_p4.Pt()>90)  )){

															if(jet1_p4.Pt()>90) h_CSV_T->Fill(Jet_btagCMVA[j]);
															if(jet2_p4.Pt()>90) h_CSV_T->Fill(Jet_btagCMVA[k]);
															if(jet3_p4.Pt()>90) h_CSV_T->Fill(Jet_btagCMVA[l]);
															if(jet4_p4.Pt()>90) h_CSV_T->Fill(Jet_btagCMVA[m]);

															//if((diJet2_p4.M()<160 && diJet2_p4.M()> 90.) && (diJet1_p4.M()<160 && diJet1_p4.M()> 90.)) pass125125 =1;
															//if((diJet2_p4.M()<125 && diJet2_p4.M()> 55.) && (diJet1_p4.M()<125 && diJet1_p4.M()>  55.)) pass9090=1;

															foundHH=true;	
															H1jet1_i=j;
															H1jet2_i=k;
															H2jet1_i=l;
															H2jet2_i=m;
															ajet1_p4.SetPtEtaPhiM(Jet_pt[H1jet1_i], Jet_eta[H1jet1_i], Jet_phi[H1jet1_i], Jet_mass[H1jet1_i]);
															ajet2_p4.SetPtEtaPhiM(Jet_pt[H1jet2_i], Jet_eta[H1jet2_i], Jet_phi[H1jet2_i], Jet_mass[H1jet2_i]);
															ajet3_p4.SetPtEtaPhiM(Jet_pt[H2jet1_i], Jet_eta[H2jet1_i], Jet_phi[H2jet1_i], Jet_mass[H2jet1_i]);
															ajet4_p4.SetPtEtaPhiM(Jet_pt[H2jet2_i], Jet_eta[H2jet2_i], Jet_phi[H2jet2_i], Jet_mass[H2jet2_i]);


															min_diff=diff;
														} // Loop over 4th jet
													} // Conditions on 3rd jet
												} // Loop over 3rd jet
											} // Conditions on 2nd jet
										}
									} // Loop over 1st jet
								}
							}
						}
					}  if (foundHH)
					{	
						TLorentzVector diJet2_p4, diJet1_p4;
						h_angle_bbar->Fill(deltaR1);
						h_angle_bbar->Fill(deltaR2);
						diJet2_p4=ajet3_p4+ajet4_p4;
						diJet1_p4=ajet1_p4+ajet2_p4;
						H1_p4=diJet1_p4;
						H2_p4=diJet2_p4;
						float h11 = H1_p4.M();
						float h22 = H2_p4.M();
						h_H1_mass->Fill(h11);
						h_H2_mass->Fill(h22);
						h_Xmass->Fill((H1_p4+H2_p4).M());
						int region=withinRegion(h11, h22, 17.5, 37.5, 125, 125);
						if(region==0) h_Xmass_comb->Fill((H1_p4+H2_p4).M());
						h_H1_H2_mass->Fill(h11, h22);

						outtree->Fill();			
						++nCut4;
					} // found HH
				} // nCJets>3
			} // vType==5
		} // Trigger
	} // Event loop

	outtree->Write();
	outfile->Close();

	std::string histfilename="Histograms_"+sample+".root";
	TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
	h_nCand->Write();
	h_nCand_true->Write();
	h_Xmass_right->Write();
	h_Xmass_wrong->Write();	
	h_Xmass_first->Write();
	h_Xmass_other->Write();
	h_nJets->Write();
	h_n8Jets->Write();
	h_nPV->Write();
	h_nPV_weighted->Write();
	h_JetpT1->Write();
	h_JetpT2->Write();
	h_JetpT3->Write();
	h_JetpT4->Write();
	h_MET->Write();
	h_MET_sig->Write();
	h_HT->Write();
	h_H1_mass->Write();
	h_H1_pT->Write();
	h_H2_mass->Write();
	h_H2_pT->Write();
	h_HH_mass_diff_cand->Write();
	h_HH_massNorm_diff->Write();
	h_genX_mass->Write();
	h_genH1_mass->Write();
	h_genH2_mass->Write();
	h_H1_H2_mass->Write();
	h_HH_CSV->Write();
	h_Xmass_comb->Write();
	h_Xmass->Write();
	h_Xmass2->Write(); 
	h_Cuts->Write();
	h_CSV1->Write();
	h_CSV2->Write();
	h_CSV3->Write();
	h_CSV4->Write();
	h_CMVA1->Write();
	h_CMVA2->Write();
	h_CMVA3->Write();
	h_CMVA4->Write();
	h_CSV_T->Write(); 
	h_angle_bbar->Write();
	h_angle_hh->Write();
	tFile->Write();
	tFile->Close();
	std::cout<<"Wrote output file "<<histfilename<<std::endl;
	std::cout<<"Number of events at the end of step 2 = "<<nCut0<<std::endl;
	std::cout<<"Number of events after trigger = "<<(float)nCut1/nCut0<<std::endl;
	//std::cout<<"Number of events after trigger 50 = "<<(float)nCut50/nCut0<<std::endl;
	//std::cout<<"Number of events after trigger  csv = "<<(float)nCutT/nCut0<<std::endl;
	//std::cout<<"Number of events after trigger  or = "<<(float)nCutOR/nCut0<<std::endl;	
	std::cout<<"Number of events after Vtype==0 = "<<(float)nCut2/nCut0<<std::endl;
	std::cout<<"Number of events after nCJets>3 = "<<(float)nCut3/nCut2<<std::endl;
	std::cout<<"Number of events after finding HH kinematic candidate = "<<(float)nCut4/nCut3<<" ncut3 "<<nCut3<<std::endl;
	//			 std::cout<<"Number of matched events "<<(float)nCutGen/nCut4<<"   "<<contComb<<"  nCut3 : "<<nCut3<<"  fail "<<ncutFail<<" ratio "<<(float)ncutFail/nCutGen<< std::endl; 	
	std::cout<<"Number of candidates "<<(float)nCandidate<<"Numeber of event passing HHfound criterion "<<nCut4<<std::endl;	
	return 0;
}




