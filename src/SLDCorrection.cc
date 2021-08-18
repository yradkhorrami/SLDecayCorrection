#include "SLDCorrection.h"
#include <iostream>
#include <EVENT/LCCollection.h>
#include "EVENT/LCCollection.h"
#include "IMPL/LCCollectionVec.h"
#include "EVENT/MCParticle.h"
#include "EVENT/Vertex.h"
#include "EVENT/ReconstructedParticle.h"
#include <IMPL/ReconstructedParticleImpl.h>
#include "IMPL/ParticleIDImpl.h"
#include "UTIL/PIDHandler.h"
#include "marlin/VerbosityLevels.h"
#include <GeometryUtil.h>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1I.h"
#include "TH2I.h"
#include "TF1.h"
#include "TTree.h"
#include "TPaveStats.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TRatioPlot.h"
#include "TAxis.h"
#include "TLine.h"

#include "visibleFourMomentum.h"
#include "flightDirection.h"

using namespace lcio ;
using namespace marlin ;

SLDCorrection aSLDCorrection;

SLDCorrection::SLDCorrection() :

	Processor("SLDCorrection"),
	m_nRun(0),
	m_nEvt(0),
	m_nRunSum(0),
	m_nEvtSum(0),
	m_Bfield(0.f),
	c(0.),
	mm2m(0.),
	eV2GeV(0.),
	eB(0.),
	foundFlightDirection(true),
	m_nTauSLDecay(0),
	m_nTauNeutrino(0),
	m_nNeutrino(0),
	m_nChargedPFOwoTrack(0),
	n_NuPxResidual(0),
	n_NuPyResidual(0),
	n_NuPzResidual(0),
	n_NuEResidual(0),
	n_NuPxNormalizedResidual(0),
	n_NuPyNormalizedResidual(0),
	n_NuPzNormalizedResidual(0),
	n_NuENormalizedResidual(0)
{
	_description = "SLDCorrection finds semi-leptonic decays within jets and performs a correction to 4-momentum of the jet due to the missing neutrino(s)";

	registerInputCollection(	LCIO::MCPARTICLE,
					"MCParticleCollection" ,
					"Name of the MCParticle collection"  ,
					m_mcParticleCollection,
					std::string("MCParticle")
				);

	registerInputCollection(	LCIO::RECONSTRUCTEDPARTICLE,
					"PfoCollection",
					"Name of input pfo collection",
					m_inputPfoCollection,
					std::string("PandoraPFOs")
				);

	registerInputCollection(	LCIO::RECONSTRUCTEDPARTICLE,
					"JetCollection",
					"Name of input jet collection",
					m_inputJetCollection,
					std::string("Durham_nJets")
				);

	registerInputCollection(	LCIO::MCPARTICLE,
					"jetFlavour" ,
					"Name of input Jet Flavour collection",
					m_jetFlavour ,
					std::string("jetFlavour")
				);

	registerInputCollection(	LCIO::VERTEX,
					"PrimaryVertex",
					"Name of Primary Vertex Collection",
					m_inputPrimaryVertex,
					std::string("PrimaryVertex")
				);

	registerInputCollection(	LCIO::VERTEX,
					"BuildUpVertex",
					"Name of BuildUp Vertex Collection",
					m_inputBuildUpVertex,
					std::string("BuildUpVertex")
				);

	registerInputCollection(	LCIO::LCRELATION,
					"RecoMCTruthLinkCollection",
					"Name of input RecoMCTruthLink Collection",
					m_RecoMCTruthLinkCollection,
					std::string("RecoMCTruthLinkCollection")
				);

	registerInputCollection(	LCIO::LCRELATION,
					"MCTruthRecoLinkCollection",
					"Name of input MCTruthRecoLink Collection",
					m_MCTruthRecoLinkCollection,
					std::string("MCTruthRecoLinkCollection")
				);

	registerInputCollection(	LCIO::LCRELATION,
					"TrackMCTruthLinkCollection",
					"Name of input TrackMCTruthLink Collection",
					m_TrackMCTruthLinkCollection,
					std::string("MarlinTrkTracksMCTruthLink")
				);

	registerInputCollection(	LCIO::LCRELATION,
					"MCTruthTrackLinkCollection",
					"Name of input MCTruthTrackLink Collection",
					m_MCTruthTrackLinkCollection,
					std::string("MCTruthMarlinTrkTracksLink")
				);

	registerInputCollection(	LCIO::LCRELATION,
					"ClusterMCTruthLinkCollection",
					"Name of input m_ClusterMCTruthLink Collection",
					m_ClusterMCTruthLinkCollection,
					std::string("ClusterMCTruthLink")
				);

	registerInputCollection(	LCIO::LCRELATION,
					"MCTruthClusterLinkCollection",
					"Name of input MCTruthClusterLink Collection",
					m_MCTruthClusterLinkCollection,
					std::string("MCTruthClusterLink")
				);

	registerOutputCollection(	LCIO::RECONSTRUCTEDPARTICLE,
					"SLDNeutrinoCollection",
					"Name of semi-leptonic Neutrino collection",
					m_SLDNuCollection,
					std::string("SLDNeutrinoCollection")
				);

	registerProcessorParameter(	"includbJets",
					"Include b-jets for semi-decay correction",
					m_includbJets,
					bool(false)
				);

	registerProcessorParameter(	"includcJets",
					"Include c-jets for semi-decay correction",
					m_includcJets,
					bool(false)
				);

	registerProcessorParameter(	"includgJets",
					"Include g-jets for semi-decay correction",
					m_includgJets,
					bool(false)
				);

	registerProcessorParameter(	"includOthers",
					"Include Other final states for semi-decay correction",
					m_includOthers,
					bool(false)
				);

	registerProcessorParameter(	"includeBSLD",
					"do correction for semi-leptonic decays of B-Hadrons",
					m_includeBSLD,
					bool(true)
				);

	registerProcessorParameter(	"includeCSLD",
					"do correction for semi-leptonic decays of C-Hadrons",
					m_includeCSLD,
					bool(true)
				);

	registerProcessorParameter(	"includeTSLD",
					"do correction for semi-leptonic decays of Tau-Leptons",
					m_includeTSLD,
					bool(true)
				);

	registerProcessorParameter(	"cheatSLDLeptons",
					"Cheat semi-leptonic decays lepton from MCTruth",
					m_cheatSLDLeptons,
					bool(true)
				);

	registerProcessorParameter(	"cheatFlightDirection",
					"Cheat Flight direction of mother hadron",
					m_cheatFlightDirection,
					bool(true)
				);

	registerProcessorParameter(	"useJetAxisAsFlightDirection",
					"use jet axis as flight direction of mother hadron",
					m_useJetAxisAsFlightDirection,
					bool(true)
				);

	registerProcessorParameter(	"considerParentCharge",
					"Consider charge of parent hadron for calculating flight direction",
					m_considerParentCharge,
					bool(true)
				);

	registerProcessorParameter(	"cheatVertices",
					"Cheat vertices of mother hadron for calculating flight direction",
					m_cheatVertices,
					bool(true)
				);

	registerProcessorParameter(	"cheatLepton4momentum",
					"Cheat FourMomentum of lepton in semi-leptonic decays",
					m_cheatLepton4momentum,
					bool(true)
				);

	registerProcessorParameter(	"cheatCharged4momentum",
					"Cheat FourMomentum of charged visibles in semi-leptonic decays",
					m_cheatCharged4momentum,
					bool(true)
				);

	registerProcessorParameter(	"cheatNeutral4momentum",
					"Cheat FourMomentum of neutral visibles in semi-leptonic decays",
					m_cheatNeutral4momentum,
					bool(true)
				);

	registerProcessorParameter(	"nIterFlightDirCorrection",
					"Number of iterations for correcting flight direction of CHARGED parent hadron",
					m_nIterFlightDirCorrection,
					int(1)
				);

	registerProcessorParameter(	"recoFourMomentumOfVisibles",
					"0: get 4p from linked track/cluster to MCP, 1: get 4p from PFO with linked track/cluster, 2: get 4p from linked PFO",
					m_recoFourMomentumOfVisibles,
					int(0)
				);

	registerProcessorParameter(	"fillRootTree",
					"Fill root tree to check processor performance",
					m_fillRootTree,
					bool(true)
				);

	registerProcessorParameter(	"RootFile",
	                                "Name of the output root file",
					m_rootFile,
					std::string("Output.root")
				);

}


void SLDCorrection::init()
{

	streamlog_out(DEBUG) << "	init called  " << std::endl;
	m_Bfield = 3.5;
	c = 2.99792458e8;
	mm2m = 1e-3;
	eV2GeV = 1e-9;
	eB = m_Bfield * c * mm2m * eV2GeV;
	printParameters();
	if ( m_fillRootTree )
	{
		m_pTFile = new TFile(m_rootFile.c_str(), "recreate");
		m_pTTree = new TTree("SLDCorrection", "SLDCorrection");
		m_pTTree->SetDirectory(m_pTFile);
		m_pTTree->Branch("event", &m_nEvt, "event/I");
		m_pTTree->Branch("nTauSLDecay",&m_nTauSLDecay,"nTauSLDecay/I");
		m_pTTree->Branch("nTauNeutrino",&m_nTauNeutrino,"nTauNeutrino/I");
		m_pTTree->Branch("nNeutrino",&m_nNeutrino,"nNeutrino/I");
		m_pTTree->Branch("jetFlavourPDG",&m_jetFlavourPDG);
		m_pTTree->Branch("nSLD_chargedMCPwoTrack",&m_nSLD_chargedMCPwoTrack);
		m_pTTree->Branch("GenStatParentHadron",&m_GenStatParentHadron);
		m_pTTree->Branch("ChargeParentHadron",&m_ChargeParentHadron);
		m_pTTree->Branch("foundRecoLepton",&m_foundRecoLepton);
		m_pTTree->Branch("foundBuildUpVertex",&m_foundBuildUpVertex);
		m_pTTree->Branch("foundRecoLeptonInBuildUpVertex",&m_foundRecoLeptonInBuildUpVertex);
		m_pTTree->Branch("foundRecoLeptonInPrimaryVertex",&m_foundRecoLeptonInPrimaryVertex);
		m_pTTree->Branch("lostChargedMCP_CosTheta",&m_lostChargedMCP_CosTheta);
		m_pTTree->Branch("lostChargedMCP_Energy",&m_lostChargedMCP_Energy);
		m_pTTree->Branch("lostChargedMCP_Pt",&m_lostChargedMCP_Pt);
		m_pTTree->Branch("SLDecayXi", &m_SLDecayXi);
		m_pTTree->Branch("SLDecayYi", &m_SLDecayYi);
		m_pTTree->Branch("SLDecayZi", &m_SLDecayZi);
		m_pTTree->Branch("SLDecayRi", &m_SLDecayRi);
		m_pTTree->Branch("SLDecayXf", &m_SLDecayXf);
		m_pTTree->Branch("SLDecayYf", &m_SLDecayYf);
		m_pTTree->Branch("SLDecayZf", &m_SLDecayZf);
		m_pTTree->Branch("SLDecayRf", &m_SLDecayRf);
		m_pTTree->Branch("trueNuPx", &m_trueNuPx);
		m_pTTree->Branch("trueNuPy", &m_trueNuPy);
		m_pTTree->Branch("trueNuPz", &m_trueNuPz);
		m_pTTree->Branch("trueNuE", &m_trueNuE);
		m_pTTree->Branch("recoNuCloseInitialPx", &m_recoNuCloseInitialPx);
		m_pTTree->Branch("recoNuCloseInitialPy", &m_recoNuCloseInitialPy);
		m_pTTree->Branch("recoNuCloseInitialPz", &m_recoNuCloseInitialPz);
		m_pTTree->Branch("recoNuCloseInitialE", &m_recoNuCloseInitialE);
		m_pTTree->Branch("recoNuClosePx", &m_recoNuClosePx);
		m_pTTree->Branch("recoNuClosePy", &m_recoNuClosePy);
		m_pTTree->Branch("recoNuClosePz", &m_recoNuClosePz);
		m_pTTree->Branch("recoNuCloseE", &m_recoNuCloseE);
		m_pTTree->Branch("recoNuPosPx", &m_recoNuPosPx);
		m_pTTree->Branch("recoNuPosPy", &m_recoNuPosPy);
		m_pTTree->Branch("recoNuPosPz", &m_recoNuPosPz);
		m_pTTree->Branch("recoNuPosE", &m_recoNuPosE);
		m_pTTree->Branch("recoNuNegPx", &m_recoNuNegPx);
		m_pTTree->Branch("recoNuNegPy", &m_recoNuNegPy);
		m_pTTree->Branch("recoNuNegPz", &m_recoNuNegPz);
		m_pTTree->Branch("recoNuNegE", &m_recoNuNegE);
		m_pTTree->Branch("NuPxResidual", &m_NuPxResidual);
		m_pTTree->Branch("NuPyResidual", &m_NuPyResidual);
		m_pTTree->Branch("NuPzResidual", &m_NuPzResidual);
		m_pTTree->Branch("NuEResidual", &m_NuEResidual);
		m_pTTree->Branch("NuPxNormalizedResidual", &m_NuPxNormalizedResidual);
		m_pTTree->Branch("NuPyNormalizedResidual", &m_NuPyNormalizedResidual);
		m_pTTree->Branch("NuPzNormalizedResidual", &m_NuPzNormalizedResidual);
		m_pTTree->Branch("NuENormalizedResidual", &m_NuENormalizedResidual);
		m_pTTree->Branch("flightDirectionError", &m_FlightDirectionError);
		m_pTTree->Branch("distRecoLeptonToDownStreamVertex", &m_distRecoLeptonToDownStreamVertex);
		h_NuPxResidual = new TH1F( "PxResidual" , "; _{}p_{x,#nu}^{REC} - p_{x,#nu}^{MC} [GeV]; Normalized Entries / 0.1" , 2000 , -10.0 , 10.0 ); n_NuPxResidual = 0;
		h_NuPyResidual = new TH1F( "PyResidual" , "; _{}p_{y,#nu}^{REC} - p_{y,#nu}^{MC} [GeV]; Normalized Entries / 0.1" , 2000 , -10.0 , 10.0 ); n_NuPyResidual = 0;
		h_NuPzResidual = new TH1F( "PzResidual" , "; _{}p_{z,#nu}^{REC} - p_{z,#nu}^{MC}  [GeV]; Normalized Entries / 0.1" , 2000 , -10.0 , 10.0 ); n_NuPzResidual = 0;
		h_NuEResidual = new TH1F( "EResidual" , "; _{}E_{#nu}^{REC} - E_{#nu}^{MC} [GeV]; Normalized Entries / 0.1" , 2000 , -10.0 , 10.0 ); n_NuEResidual = 0;
		h_NuPxNormalizedResidual = new TH1F( "PxNormalizedResidual" , "; ( _{}p_{x,#nu}^{REC} - p_{x,#nu}^{MC} ) / #sigma_{p_{x,#nu}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_NuPxNormalizedResidual = 0;
		h_NuPyNormalizedResidual = new TH1F( "PyNormalizedResidual" , "; ( _{}p_{y,#nu}^{REC} - p_{y,#nu}^{MC} ) / #sigma_{p_{y,#nu}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_NuPyNormalizedResidual = 0;
		h_NuPzNormalizedResidual = new TH1F( "PzNormalizedResidual" , "; ( _{}p_{z,#nu}^{REC} - p_{z,#nu}^{MC} ) / #sigma_{p_{z,#nu}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_NuPzNormalizedResidual = 0;
		h_NuENormalizedResidual = new TH1F( "ENormalizedResidual" , "; ( _{}E_{#nu}^{REC} - E_{#nu}^{MC} ) / #sigma_{E_{#nu}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_NuENormalizedResidual = 0;
		h_recoNuPx_mcNuPx = new TH2F( "p_{x}^{#nu}" , "; _{}p_{x,#nu}^{MC} [GeV] ;  _{}p_{x,#nu}^{REC} [GeV]" , 120 , -60.0 , 60.0 , 120 , -60.0 , 60.0 );
		h_recoNuPy_mcNuPy = new TH2F( "p_{y}^{#nu}" , "; _{}p_{y,#nu}^{MC} [GeV] ;  _{}p_{y,#nu}^{REC} [GeV]" , 120 , -60.0 , 60.0 , 120 , -60.0 , 60.0 );
		h_recoNuPz_mcNuPz = new TH2F( "p_{z}^{#nu}" , "; _{}p_{z,#nu}^{MC} [GeV] ;  _{}p_{z,#nu}^{REC} [GeV]" , 120 , -60.0 , 60.0 , 120 , -60.0 , 60.0 );
		h_recoNuE_mcNuE = new TH2F( "E^{#nu}" , "; _{}E_{#nu}^{MC} [GeV] ;  _{}E_{#nu}^{REC} [GeV]" , 100 , 0.0 , 100.0 , 100 , 0.0 , 100.0 );
		h_parentPx_daughtersPx = new TH2F( "p_{x} conservation" , "; _{}Px_{parentHadron}^{MC} [GeV] ;  _{}Px_{daughters}^{MC} [GeV]" , 200 , -100.0 , 100.0 , 200 , -100.0 , 100.0 );
		h_parentPy_daughtersPy = new TH2F( "p_{y} conservation" , "; _{}Py_{parentHadron}^{MC} [GeV] ;  _{}Py_{daughters}^{MC} [GeV]" , 200 , -100.0 , 100.0 , 200 , -100.0 , 100.0 );
		h_parentPz_daughtersPz = new TH2F( "p_{z} conservation" , "; _{}Pz_{parentHadron}^{MC} [GeV] ;  _{}Pz_{daughters}^{MC} [GeV]" , 200 , -100.0 , 100.0 , 200 , -100.0 , 100.0 );
		h_parentE_daughtersE = new TH2F( "E conservation" , "; _{}E_{parentHadron}^{MC} [GeV] ;  _{}E_{daughters}^{MC} [GeV]" , 100 , 0.0 , 100.0 , 100 , 0.0 , 100.0 );
		h_recoPFOLinkedToElectron_Type = new TH1I( "PFOTypeofTrueElectron" , "; PFO Type" , 8 , 0 , 8 );
		h_recoPFOLinkedToElectron_Type->GetXaxis()->SetBinLabel(1,"e^{#pm}");
		h_recoPFOLinkedToElectron_Type->GetXaxis()->SetBinLabel(2,"#mu^{#pm}");
		h_recoPFOLinkedToElectron_Type->GetXaxis()->SetBinLabel(3,"#pi^{#pm}");
		h_recoPFOLinkedToElectron_Type->GetXaxis()->SetBinLabel(4,"#gamma");
		h_recoPFOLinkedToElectron_Type->GetXaxis()->SetBinLabel(5,"K^{0}_{S}");
		h_recoPFOLinkedToElectron_Type->GetXaxis()->SetBinLabel(6,"#Lambda");
		h_recoPFOLinkedToElectron_Type->GetXaxis()->SetBinLabel(7,"Other");
		h_recoPFOLinkedToElectron_Type->GetXaxis()->SetBinLabel(8,"Not Found");
		h_recoPFOLinkedToMuon_Type = new TH1I( "PFOTypeofTrueMuon" , "; PFO Type" , 8 , 0 , 8 );
		h_recoPFOLinkedToMuon_Type->GetXaxis()->SetBinLabel(1,"e^{#pm}");
		h_recoPFOLinkedToMuon_Type->GetXaxis()->SetBinLabel(2,"#mu^{#pm}");
		h_recoPFOLinkedToMuon_Type->GetXaxis()->SetBinLabel(3,"#pi^{#pm}");
		h_recoPFOLinkedToMuon_Type->GetXaxis()->SetBinLabel(4,"#gamma");
		h_recoPFOLinkedToMuon_Type->GetXaxis()->SetBinLabel(5,"K^{0}_{S}");
		h_recoPFOLinkedToMuon_Type->GetXaxis()->SetBinLabel(6,"#Lambda");
		h_recoPFOLinkedToMuon_Type->GetXaxis()->SetBinLabel(7,"Other");
		h_recoPFOLinkedToMuon_Type->GetXaxis()->SetBinLabel(8,"Not Found");
		h_SLDecayOrder = new TH1I( "SLDecayOrder" , "; SLDecay Type" , 3 , 0 , 3 );
		h_SLDecayOrder->GetXaxis()->SetBinLabel(1,"upStream");
		h_SLDecayOrder->GetXaxis()->SetBinLabel(2,"primary");
		h_SLDecayOrder->GetXaxis()->SetBinLabel(3,"downStream");
		h_foundVertex = new TH2I( "Viertices" , "; primary vertex ; secondary vertex" , 2 , 0 , 2 , 2 , 0 , 2 );
		h_foundVertex->GetXaxis()->SetBinLabel(1,"vertex not found");
		h_foundVertex->GetXaxis()->SetBinLabel(2,"vertex found");
		h_foundVertex->GetYaxis()->SetBinLabel(1,"vertex not found");
		h_foundVertex->GetYaxis()->SetBinLabel(2,"vertex found");
		h_secondaryVertex = new TH1I( "secondary vertices" , ";" , 6 , 0 , 6 );
		h_secondaryVertex->GetXaxis()->SetBinLabel(1,"lep in BUp vtx");
		h_secondaryVertex->GetXaxis()->SetBinLabel(2,"lep in Prim vtx");
		h_secondaryVertex->GetXaxis()->SetBinLabel(3,"SLD with downStream vtx");
		h_secondaryVertex->GetXaxis()->SetBinLabel(4,"Sec. vtx not found");
		h_secondaryVertex->GetXaxis()->SetBinLabel(5,"reco lep not found");
		h_secondaryVertex->GetXaxis()->SetBinLabel(6,"other(?)");
		h_parentHadronCharge = new TH1I( "parentHadronCharge" , "; Parent Hadron Charge" , 5 , 0 , 5 );
		h_parentHadronCharge->GetXaxis()->SetBinLabel(1,"-2");
		h_parentHadronCharge->GetXaxis()->SetBinLabel(2,"-1");
		h_parentHadronCharge->GetXaxis()->SetBinLabel(3,"0");
		h_parentHadronCharge->GetXaxis()->SetBinLabel(4,"1");
		h_parentHadronCharge->GetXaxis()->SetBinLabel(5,"2");
		h_MCPTracks = new TH1I( "chargedMCPTracks" , ";" , 2 , 0 , 2 );
		h_MCPTracks->GetXaxis()->SetBinLabel(1,"track found");
		h_MCPTracks->GetXaxis()->SetBinLabel(2,"track lost");
		h_MCPTracks_Eweighted = new TH1I( "chargedMCPTracks" , "Energy weighted;" , 2 , 0 , 2 );
		h_MCPTracks_Eweighted->GetXaxis()->SetBinLabel(1,"track found");
		h_MCPTracks_Eweighted->GetXaxis()->SetBinLabel(2,"track lost");
		h_MCPTracks_Ptweighted = new TH1I( "chargedMCPTracks" , "p_{T} weighted;" , 2 , 0 , 2 );
		h_MCPTracks_Ptweighted->GetXaxis()->SetBinLabel(1,"track found");
		h_MCPTracks_Ptweighted->GetXaxis()->SetBinLabel(2,"track lost");
		h_FlightDirectionError = new TH1F ( "Flight Direction Error" , "; cos #alpha" , 200 , -1.0 , 1.0 );
		h_distRecoLeptonToDownStreamVertex = new TH1F( "distance of RecoLepton to DownStreamVertex" , ";r^{3D} [mm]", 100 , 0.0 , 10.0 );
	}

}

void SLDCorrection::Clear()
{
	m_jetFlavourPDG.clear();	//4: c-jet, 5: b-jet, 21: g-jet, 0:Other
	m_nSLD_chargedMCPwoTrack.clear();
	m_GenStatParentHadron.clear();
	m_ChargeParentHadron.clear();
	m_foundRecoLepton.clear();
	m_foundBuildUpVertex.clear();
	m_foundRecoLeptonInBuildUpVertex.clear();
	m_foundRecoLeptonInPrimaryVertex.clear();
	m_lostChargedMCP_CosTheta.clear();
	m_lostChargedMCP_Energy.clear();
	m_lostChargedMCP_Pt.clear();
	m_nTauSLDecay = 0;
	m_nTauNeutrino = 0;
	m_nNeutrino = 0;
	m_nChargedPFOwoTrack = 0;
	m_SLDecayXi.clear();
	m_SLDecayYi.clear();
	m_SLDecayZi.clear();
	m_SLDecayRi.clear();
	m_SLDecayXf.clear();
	m_SLDecayYf.clear();
	m_SLDecayZf.clear();
	m_SLDecayRf.clear();
	m_trueNuPx.clear();
	m_trueNuPy.clear();
	m_trueNuPz.clear();
	m_trueNuE.clear();
	m_recoNuCloseInitialPx.clear();
	m_recoNuCloseInitialPy.clear();
	m_recoNuCloseInitialPz.clear();
	m_recoNuCloseInitialE.clear();
	m_recoNuClosePx.clear();
	m_recoNuClosePy.clear();
	m_recoNuClosePz.clear();
	m_recoNuCloseE.clear();
	m_recoNuPosPx.clear();
	m_recoNuPosPy.clear();
	m_recoNuPosPz.clear();
	m_recoNuPosE.clear();
	m_recoNuNegPx.clear();
	m_recoNuNegPy.clear();
	m_recoNuNegPz.clear();
	m_recoNuNegE.clear();
	m_NuPxResidual.clear();
	m_NuPyResidual.clear();
	m_NuPzResidual.clear();
	m_NuEResidual.clear();
	m_NuPxNormalizedResidual.clear();
	m_NuPyNormalizedResidual.clear();
	m_NuPzNormalizedResidual.clear();
	m_NuENormalizedResidual.clear();
	m_FlightDirectionError.clear();
	m_distRecoLeptonToDownStreamVertex.clear();

}

void SLDCorrection::processRunHeader()
{
	m_nRun = 0;
	m_nEvt = 0;
	++m_nRunSum;
}

void SLDCorrection::processEvent( EVENT::LCEvent *pLCEvent )
{
	m_nRun = pLCEvent->getRunNumber();
	m_nEvt = pLCEvent->getEventNumber();

	LCCollection *MCParticleCollection{};
	LCCollection *jetFlacourCollection{};
	int nTauNeutrino = 0;
	int m_bJet = 0;
	int m_cJet = 0;
	int m_gJet = 0;
	int m_other = 0;
	int jetFlavourPDG = 0;
	++m_nEvtSum;
	this->Clear();
	streamlog_out(MESSAGE) << "" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////	Processing event 	" << m_nEvt << "	////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;

	try
	{
		jetFlacourCollection = pLCEvent->getCollection( m_jetFlavour );
		m_bJet = jetFlacourCollection->getParameters().getIntVal("isDecayedTob"); if ( m_bJet == 1 ) jetFlavourPDG = 5;
		m_cJet = jetFlacourCollection->getParameters().getIntVal("isDecayedToc"); if ( m_cJet == 1 ) jetFlavourPDG = 4;
		m_gJet = jetFlacourCollection->getParameters().getIntVal("isDecayedTog"); if ( m_gJet == 1 ) jetFlavourPDG = 21;
		m_other = jetFlacourCollection->getParameters().getIntVal("isDecayedToother"); if ( m_other == 1 ) jetFlavourPDG = 0;
		if ( m_bJet == 1 && !m_includbJets ) return;
		if ( m_cJet == 1 && !m_includcJets ) return;
		if ( m_gJet == 1 && !m_includgJets ) return;
		if ( m_other == 1 && !m_includOthers ) return;

		MCParticleCollection = pLCEvent->getCollection( m_mcParticleCollection );
		int nMCP = MCParticleCollection->getNumberOfElements();
		for ( int i_mcp = 0 ; i_mcp < nMCP ; ++i_mcp )
		{
			MCParticle *testLepton = dynamic_cast<EVENT::MCParticle*>( MCParticleCollection->getElementAt( i_mcp ) );
			if ( abs( testLepton->getPDG() ) == 16 && ( testLepton->getGeneratorStatus() ) == 1 ) ++nTauNeutrino;
			bool primarySLDecay = false;
			bool downStreamSLDecay = false;
			bool upStreamSLDecay = false;
			bool isBHadronSLDecay = false;
			bool isCHadronSLDecay = false;
			bool isTauLeptonSLDecay = false;
			if ( ( abs( testLepton->getPDG() ) == 11 || abs( testLepton->getPDG() ) == 13 || abs( testLepton->getPDG() ) == 15 ) && ( testLepton->getGeneratorStatus() ) == 1 )
			{
				for ( long unsigned int i_parent = 0 ; i_parent < ( testLepton->getParents() ).size() ; ++i_parent )
				{
					MCParticle *parent = testLepton->getParents()[ i_parent ];
					primarySLDecay = hasPrimarySLDecay( parent );
					if ( primarySLDecay ) downStreamSLDecay = hasDownStreamSLDecay( parent );
					if ( primarySLDecay ) upStreamSLDecay = hasUpStreamSLDecay( parent );
				}
				if ( primarySLDecay )
				{
					std::vector< TLorentzVector > recoNeutrinoFourMomentum;
					std::vector< std::vector< float > > recoNeutrinoCovMat;
					isBHadronSLDecay = checkBHadronSLDecay( testLepton );
					isCHadronSLDecay = checkCHadronSLDecay( testLepton );
					isTauLeptonSLDecay = checkTauLeptonSLDecay( testLepton );
					if ( isBHadronSLDecay && !m_includeBSLD ) continue;
					if ( isCHadronSLDecay && !m_includeCSLD ) continue;
					if ( isTauLeptonSLDecay && !m_includeTSLD ) continue;
					streamlog_out(DEBUG3) << "" << std::endl;
					streamlog_out(DEBUG3) << "	<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
					streamlog_out(DEBUG3) << "	<<<<<<<<<<<<<<<<<<<<<<<<<<< Found a primary semi-leptonic decay >>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
					if ( downStreamSLDecay )
					{
						streamlog_out(DEBUG3) << "	There is/are downstream semi-leptonic(s) decay in primary semi-leptonic decay products" << std::endl;
						h_SLDecayOrder->Fill( 2.5 );
					}
					if ( upStreamSLDecay )
					{
						streamlog_out(DEBUG3) << "	There is/are upstream semi-leptonic(s) decay in primary semi-leptonic decay products" << std::endl;
						h_SLDecayOrder->Fill( 0.5 );
					}
					if ( !downStreamSLDecay && !upStreamSLDecay )
					{
						streamlog_out(DEBUG3) << "	<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
						streamlog_out(DEBUG3) << "	<<<<<<<<<<<<<<<< There are no upstream and downstream semi-leptonic decay >>>>>>>>>>>>>>>>>" << std::endl;
						h_SLDecayOrder->Fill( 1.5 );
						m_jetFlavourPDG.push_back( jetFlavourPDG );
						doSLDCorrection( pLCEvent , testLepton );
					}
				}
			}
		}
		m_nTauNeutrino = nTauNeutrino;
		m_pTTree->Fill();
	}
	catch(DataNotAvailableException &e)
        {
        	streamlog_out(MESSAGE) << "	Input collection not found in event " << m_nEvt << std::endl;
        }

}

bool SLDCorrection::hasPrimarySLDecay( MCParticle *parentHadron )
{
	bool hasSLDecay = false;
	if ( parentHadron->getGeneratorStatus() == 2 && ( floor( abs( parentHadron->getPDG() ) / 100 ) == 5 || ( floor( abs( parentHadron->getPDG() ) / 1000 ) == 5 ) || floor( abs( parentHadron->getPDG() ) / 100 ) == 4 || ( floor( abs( parentHadron->getPDG() ) / 1000 ) == 4 ) || ( abs( parentHadron->getPDG() ) == 15 ) ) )
	{
		for ( long unsigned int i_daughter = 0 ; i_daughter < ( parentHadron->getDaughters() ).size() ; ++i_daughter )
		{
			MCParticle *daughter = parentHadron->getDaughters()[ i_daughter ];
			if ( daughter->getGeneratorStatus() == 1 )
			{
				if ( abs( daughter->getPDG() ) == 11 || abs( daughter->getPDG() ) == 13 || abs( daughter->getPDG() ) == 15 )
				{
					for ( long unsigned int i_NuCandidate = 0 ; i_NuCandidate < ( parentHadron->getDaughters() ).size() ; ++i_NuCandidate )
					{
						MCParticle *secondDaughter = parentHadron->getDaughters()[ i_NuCandidate ];
						if ( ( abs( secondDaughter->getPDG() ) == abs( daughter->getPDG() ) + 1 ) && secondDaughter->getGeneratorStatus() == 1 )
						{
							hasSLDecay = true;
						}
					}
				}
			}
			else if( abs( daughter->getPDG() ) == 15 )
			{
				for ( long unsigned int i_NuCandidate = 0 ; i_NuCandidate < ( parentHadron->getDaughters() ).size() ; ++i_NuCandidate )
				{
					MCParticle *secondDaughter = parentHadron->getDaughters()[ i_NuCandidate ];
					if ( ( abs( secondDaughter->getPDG() ) == abs( daughter->getPDG() ) + 1 ) && secondDaughter->getGeneratorStatus() == 1 )
					{
						hasSLDecay = true;
						m_nTauSLDecay += 1;
					}
				}
			}
		}
	}
	return hasSLDecay;
}

bool SLDCorrection::hasDownStreamSLDecay( MCParticle *parentHadron )
{
	bool hasSLDecay = false;
	bool primarySLDecay = false;
	bool downStreamSLDecay = false;
	for ( long unsigned int i_primaryDaughter = 0 ; i_primaryDaughter < ( parentHadron->getDaughters() ).size() ; ++i_primaryDaughter )
	{
		MCParticle *primaryDaughter = parentHadron->getDaughters()[ i_primaryDaughter ];
		primarySLDecay = primarySLDecay || hasPrimarySLDecay( primaryDaughter );
		downStreamSLDecay = downStreamSLDecay || hasDownStreamSLDecay( primaryDaughter );
	}
	hasSLDecay = primarySLDecay || downStreamSLDecay ;
	return hasSLDecay;
}

bool SLDCorrection::hasUpStreamSLDecay( MCParticle *parentHadron )
{
	bool hasSLDecay = false;
	bool primarySLDecay = false;
	bool upStreamSLDecay = false;
	for ( long unsigned int i_upperParent = 0 ; i_upperParent < ( parentHadron->getParents() ).size() ; ++i_upperParent )
	{
		MCParticle *upperParent = parentHadron->getParents()[ i_upperParent ];
		primarySLDecay = primarySLDecay || hasPrimarySLDecay( upperParent );
		upStreamSLDecay = upStreamSLDecay || hasUpStreamSLDecay( upperParent );
	}
	hasSLDecay = primarySLDecay || upStreamSLDecay ;
	return hasSLDecay;
}

bool SLDCorrection::checkBHadronSLDecay( MCParticle *SLDLepton )
{
	bool isBHadronSLDecay = false;
	MCParticle *parentHadron = SLDLepton->getParents()[ 0 ];
	if ( floor( abs( parentHadron->getPDG() ) / 100 ) == 5 || floor( abs( parentHadron->getPDG() ) / 1000 ) == 5 ) isBHadronSLDecay = true;
	return isBHadronSLDecay;
}

bool SLDCorrection::checkCHadronSLDecay( MCParticle *SLDLepton )
{
	bool isCHadronSLDecay = false;
	MCParticle *parentHadron = SLDLepton->getParents()[ 0 ];
	if ( floor( abs( parentHadron->getPDG() ) / 100 ) == 4 || floor( abs( parentHadron->getPDG() ) / 1000 ) == 4 ) isCHadronSLDecay = true;
	return isCHadronSLDecay;
}

bool SLDCorrection::checkTauLeptonSLDecay( MCParticle *SLDLepton )
{
	bool TauLeptonSLDecay = false;
	MCParticle *parent = SLDLepton->getParents()[ 0 ];
	if ( abs( parent->getPDG() ) == 15 ) TauLeptonSLDecay = true;
	return TauLeptonSLDecay;
}

void SLDCorrection::doSLDCorrection( EVENT::LCEvent *pLCEvent , MCParticle *SLDLepton )
{
	TLorentzVector recoNeutrinoFourMomentumPos( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector recoNeutrinoFourMomentumNeg( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector recoNeutrinoFourMomentumClose( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector trueLeptonFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector trueChargedFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector trueNeutralFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector leptonFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector chargedFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector neutralFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	TVector3 flightDirection( 0.0 , 0.0 , 0.0 );
	float restCharge = 0.0;
	double parentHadronMass = 0.0;
	streamlog_out(DEBUG2) << "			     (  PDG	, Mass		, Px		, Py		, Pz		, E		, Charge	)" << std::endl;
	getLeptonFourMomentum( pLCEvent , SLDLepton  , m_cheatLepton4momentum , leptonFourMomentum , trueLeptonFourMomentum , m_RecoMCTruthLinkCollection , m_MCTruthRecoLinkCollection );
	streamlog_out(DEBUG2) << "			Used:(	" << SLDLepton->getPDG() << "	, " << leptonFourMomentum.M() << "	, " << leptonFourMomentum.Px() << "	, " << leptonFourMomentum.Py() << "	, " << leptonFourMomentum.Pz() << "	, " << leptonFourMomentum.E() << "	, " << (int)SLDLepton->getCharge() << "		)" << std::endl;
	streamlog_out(DEBUG2) << "" << std::endl;
	getChargedFourMomentum( pLCEvent , SLDLepton  , m_cheatCharged4momentum , chargedFourMomentum , trueChargedFourMomentum , restCharge , m_RecoMCTruthLinkCollection , m_MCTruthRecoLinkCollection );
	streamlog_out(DEBUG2) << "			Used:(	" << "qqq" << "	, " << chargedFourMomentum.M() << "	, " << chargedFourMomentum.Px() << "	, " << chargedFourMomentum.Py() << "	, " << chargedFourMomentum.Pz() << "	, " << chargedFourMomentum.E() << "	, " << (int)restCharge << "		)" << std::endl;
	streamlog_out(DEBUG2) << "" << std::endl;
	getNeutralFourMomentum( pLCEvent , SLDLepton  , m_cheatNeutral4momentum , neutralFourMomentum , trueNeutralFourMomentum , m_RecoMCTruthLinkCollection , m_MCTruthRecoLinkCollection );
	streamlog_out(DEBUG2) << "			Used:(	" << "nnn" << "	, " << chargedFourMomentum.M() << "	, " << chargedFourMomentum.Px() << "	, " << chargedFourMomentum.Py() << "	, " << chargedFourMomentum.Pz() << "	, " << chargedFourMomentum.E() << "	, " << "0" << "		)" << std::endl;
	streamlog_out(DEBUG2) << "" << std::endl;
//	getParentHadronFlightDirection( pLCEvent , SLDLepton , m_inputPrimaryVertex , m_inputBuildUpVertex , m_inputJetCollection );

//	recoNeutrinoFourMomentumPos = getNeutrinoFourMomentum( flightDirection , leptonFourMomentum , chargedFourMomentum , neutralFourMomentum , parentHadronMass , +1 );
//	recoNeutrinoFourMomentumNeg = getNeutrinoFourMomentum( flightDirection , leptonFourMomentum , chargedFourMomentum , neutralFourMomentum , parentHadronMass , -1 );

//	TLorentzVector fourMomentumLepton = getLeptonFourMomentum( pLCEvent , SLDLepton );
//	TLorentzVector visibleFourMomentumCharged = getVisibleFourMomentum( pLCEvent , SLDLepton , SLDLepton->getParents()[ 0 ] , true , false );
//	TLorentzVector visibleFourMomentumNeutral = getVisibleFourMomentum( pLCEvent , SLDLepton , SLDLepton->getParents()[ 0 ] , false , true );
}

TLorentzVector SLDCorrection::getNeutrinoFourMomentum( TVector3 flightDirection , TLorentzVector FourMomentumLepton , TLorentzVector VisibleFourMomentumCharged , TLorentzVector VisibleFourMomentumNeutral , double ParentHadronMass , int solutionSign )
{
	int sign = ( solutionSign != 0 ? solutionSign / abs( solutionSign ) : 1 );
	const char *solSign = ( sign >= 0 ? "+" : "-" );
	streamlog_out(DEBUG1) << "" << std::endl;
	streamlog_out(DEBUG1) << "		--------------------------------------------" << std::endl;
	streamlog_out(DEBUG1) << "		Calculate Neutrino 4-Momentum for " << solSign << " solution" << std::endl;
	streamlog_out(DEBUG1) << "		--------------------------------------------" << std::endl;

	flightDirection.SetMag( 1.0 );
	streamlog_out(DEBUG1) << "		flightDirection:			( " << flightDirection.Px() << "	, " << flightDirection.Py() << "	, " << flightDirection.Pz() << "	)" << std::endl;
	streamlog_out(DEBUG1) << "		Parent Hadron Mass =	 " << ParentHadronMass << std::endl;

	TLorentzVector visible_tlv	= VisibleFourMomentumCharged + VisibleFourMomentumNeutral + FourMomentumLepton;
	streamlog_out(DEBUG1) << "		Visible 4-Momentum:			( " << visible_tlv.Px() << "	, " << visible_tlv.Py() << "	, " << visible_tlv.Pz() << "	, " << visible_tlv.E() << " )" << std::endl;
	streamlog_out(DEBUG1) << "				Lepton			( " << FourMomentumLepton.Px() << "	, " << FourMomentumLepton.Py() << "	, " << FourMomentumLepton.Pz() << "	, " << FourMomentumLepton.E() << " )" << std::endl;
	streamlog_out(DEBUG1) << "				Charged		( " << VisibleFourMomentumCharged.Px() << "	, " << VisibleFourMomentumCharged.Py() << "	, " << VisibleFourMomentumCharged.Pz() << "	, " << VisibleFourMomentumCharged.E() << " )" << std::endl;
	streamlog_out(DEBUG1) << "				Neutral		( " << VisibleFourMomentumNeutral.Px() << "	, " << VisibleFourMomentumNeutral.Py() << "	, " << VisibleFourMomentumNeutral.Pz() << "	, " << VisibleFourMomentumNeutral.E() << " )" << std::endl;

	double visible_mass		= visible_tlv.M();
	streamlog_out(DEBUG1) << "		Visible Inv Mass:	" << visible_mass << std::endl;

	double visible_E		= visible_tlv.E();
	streamlog_out(DEBUG1) << "		Visible Energy:									" << visible_E << std::endl;

	TVector3 visible_p		= TVector3( visible_tlv.Px() , visible_tlv.Py() , visible_tlv.Pz() );
	streamlog_out(DEBUG1) << "		Visible Momentum:			( " << visible_p.Px() << "	, " << visible_p.Py() << "	, " << visible_p.Pz() << "	)" << std::endl;
	streamlog_out(DEBUG1) << "		(Visible Momentum).(FlightDirection):							" << visible_p.Dot( flightDirection ) << std::endl;

	TVector3 visible_p_par		= visible_p.Dot( flightDirection ) * flightDirection;
	streamlog_out(DEBUG1) << "		Visible Momentum (par):		( " << visible_p_par.Px() << "	, " << visible_p_par.Py() << "	, " << visible_p_par.Pz() << "	)" << std::endl;

	TVector3 visible_p_nor		= visible_p - visible_p_par;
	streamlog_out(DEBUG1) << "		Visible Momentum (nor):		( " << visible_p_nor.Px() << "	, " << visible_p_nor.Py() << "	, " << visible_p_nor.Pz() << "	)" << std::endl;

	double visible_E_prime		= ( pow( ParentHadronMass , 2 ) + pow( visible_mass , 2 ) ) / ( 2 * ParentHadronMass );
	streamlog_out(DEBUG1) << "		Visible Energy (prime):								" << visible_E_prime << std::endl;

	TVector3 visible_p_par_prime	= sign * sqrt( pow( ( pow( ParentHadronMass , 2 ) - pow( visible_mass , 2 ) ) / ( 2 * ParentHadronMass ) , 2 ) - visible_p_nor.Mag2() ) * flightDirection;
	if ( pow( ( pow( ParentHadronMass , 2 ) - pow( visible_mass , 2 ) ) / ( 2 * ParentHadronMass ) , 2 ) < visible_p_nor.Mag2() )
	{
		visible_p_par_prime	= std::numeric_limits<double>::min() * flightDirection;
	}
	streamlog_out(DEBUG1) << "		Visible Momentum (par-prime):		( " << visible_p_par_prime.Px() << "	, " << visible_p_par_prime.Py() << "	, " << visible_p_par_prime.Pz() << "	)" << std::endl;

	double parent_hadron_E		= ( ( visible_tlv.E() * ( pow( ParentHadronMass , 2 ) + pow( visible_mass , 2 ) ) / ( 2 * ParentHadronMass ) ) - visible_p.Dot( flightDirection ) * sign * sqrt( pow( ( pow( ParentHadronMass , 2 ) - pow( visible_mass , 2 ) ) / ( 2 * ParentHadronMass ) , 2 ) - visible_p_nor.Mag2() ) ) * ParentHadronMass / ( pow( visible_mass , 2 ) + visible_p_nor.Mag2() );
	TVector3 parent_hadron_p	= sqrt( pow( parent_hadron_E , 2 ) - pow( ParentHadronMass , 2 ) ) * flightDirection;
	streamlog_out(DEBUG1) << "		Parent Hadron Momentum:		( " << parent_hadron_p.Px() << "	, " << parent_hadron_p.Py() << "	, " << parent_hadron_p.Pz() << "	, " << parent_hadron_E << " )" << std::endl;

	double Neutrino_E		= parent_hadron_E - visible_E;
	streamlog_out(DEBUG1) << "		Neutrino Energy:									" << Neutrino_E << std::endl;

	TVector3 Neutrino_p_nor		= -1 * visible_p_nor;
	streamlog_out(DEBUG1) << "		Neutrino Momentum (nor):		( " << Neutrino_p_nor.Px() << "	, " << Neutrino_p_nor.Py() << "	, " << Neutrino_p_nor.Pz() << "	)" << std::endl;

	TVector3 Neutrino_p_par		= sqrt( pow( Neutrino_E , 2 ) - Neutrino_p_nor.Mag2() ) * flightDirection;
	streamlog_out(DEBUG1) << "		Neutrino Momentum (par):		( " << Neutrino_p_nor.Px() << "	, " << Neutrino_p_nor.Py() << "	, " << Neutrino_p_nor.Pz() << "	)" << std::endl;

	TVector3 Neutrino_p		= Neutrino_p_nor + Neutrino_p_par;
	streamlog_out(DEBUG1) << "		Neutrino Momentum:			( " << Neutrino_p_par.Px() << "	, " << Neutrino_p_par.Py() << "	, " << Neutrino_p_par.Pz() << "	)" << std::endl;

	TLorentzVector Neutrino_tlv( Neutrino_p , Neutrino_E );
	streamlog_out(DEBUG2) << "" << std::endl;
	streamlog_out(DEBUG2) << "	Reconstructed Neutrino Mass:	" << Neutrino_tlv.Mag() << std::endl;
	streamlog_out(DEBUG2) << "	Reconstructed Neutrino 4-Momentum:		( " << Neutrino_tlv.Px() << "	, " << Neutrino_tlv.Py() << "	, " << Neutrino_tlv.Pz() << "	, " << Neutrino_tlv.E() << " )" << std::endl;
	return Neutrino_tlv;
}

void SLDCorrection::check( EVENT::LCEvent *pLCEvent )
{
	LCCollection *SLDNuCollection{};
	try
	{
		SLDNuCollection = pLCEvent->getCollection(m_SLDNuCollection);
		int n_SLDNeutrinos = SLDNuCollection->getNumberOfElements();
		streamlog_out(DEBUG4) << "	CHECK : processed events: " << m_nEvtSum << " , Number of Neutrinos after semi-leptonic correction (considering n-fold ambiguity)" << n_SLDNeutrinos << std::endl;
	}
	catch(DataNotAvailableException &e)
        {
          streamlog_out(MESSAGE) << "	Input/Output collection not found in event " << m_nEvt << std::endl;
        }

}

void SLDCorrection::end()
{
	if ( m_fillRootTree )
	{
		m_pTFile->cd();
		m_pTTree->Write();
		h_SLDecayOrder->Write();

		m_pTFile->Close();
		delete m_pTFile;
	}
}
