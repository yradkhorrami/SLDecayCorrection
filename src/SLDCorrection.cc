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
#include <IMPL/VertexImpl.h>
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
#include "DDMarlinCED.h"

#include "visibleFourMomentum.h"
#include "flightDirection.h"
#include "AssignParticlestoSLD.h"
//#include "linkedPFO.h"
#include "FindParticle.h"

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
	n_SLDStatus(0),
	n_NuPxResidual(0),
	n_NuPyResidual(0),
	n_NuPzResidual(0),
	n_NuEResidual(0),
	n_NuPxNormalizedResidual(0),
	n_NuPyNormalizedResidual(0),
	n_NuPzNormalizedResidual(0),
	n_NuENormalizedResidual(0),
	n_secondaryVertex(0)
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

	registerProcessorParameter(	"vertexingScenario",
					"Scenario for finding flight direction of mother hadron: 1 = default , 2 = assign jet axis , 3 = assign flight direction of leading particle in the jet",
					m_vertexingScenario,
					int(1)
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

	registerProcessorParameter(	"cheatPIDcharged",
					"Cheat Particle ID for Charged Decay Products",
					m_cheatPIDcharged,
					bool(true)
				);

	registerProcessorParameter(	"cheatPIDneutral",
					"Cheat Particle ID for Neutral Decay Products",
					m_cheatPIDneutral,
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

	registerProcessorParameter(	"displayEvent",
					"Display recoLepton, downstraem vertex and RP and reconstructed secondary vertex in Event Display",
					m_displayEvent,
					bool(true)
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
	m_Bfield = MarlinUtil::getBzAtOrigin();
//	m_Bfield = 3.5;
	c = 2.99792458e8;
	mm2m = 1e-3;
	eV2GeV = 1e-9;
	eB = m_Bfield * c * mm2m * eV2GeV;
	printParameters();
	DDMarlinCED::init(this);

	if ( m_fillRootTree )
	{
		m_pTFile = new TFile(m_rootFile.c_str(), "recreate");
		m_pTTree1 = new TTree("SLDCorrection", "SLDCorrection");
		m_pTTree1->SetDirectory(m_pTFile);
		m_pTTree1->Branch("event", &m_nEvt, "event/I");
		m_pTTree1->Branch("nTauSLDecay",&m_nTauSLDecay,"nTauSLDecay/I");
		m_pTTree1->Branch("nTauNeutrino",&m_nTauNeutrino,"nTauNeutrino/I");
		m_pTTree1->Branch("nNeutrino",&m_nNeutrino,"nNeutrino/I");
		m_pTTree1->Branch("jetFlavourPDG",&m_jetFlavourPDG);
		m_pTTree1->Branch("nSLD_chargedMCPwoTrack",&m_nSLD_chargedMCPwoTrack);
		m_pTTree1->Branch("GenStatParentHadron",&m_GenStatParentHadron);
		m_pTTree1->Branch("ChargeParentHadron",&m_ChargeParentHadron);
		m_pTTree1->Branch("foundRecoLepton",&m_foundRecoLepton);
		m_pTTree1->Branch("foundBuildUpVertex",&m_foundBuildUpVertex);
		m_pTTree1->Branch("foundRecoLeptonInBuildUpVertex",&m_foundRecoLeptonInBuildUpVertex);
		m_pTTree1->Branch("foundRecoLeptonInPrimaryVertex",&m_foundRecoLeptonInPrimaryVertex);
		m_pTTree1->Branch("lostChargedMCP_CosTheta",&m_lostChargedMCP_CosTheta);
		m_pTTree1->Branch("lostChargedMCP_Energy",&m_lostChargedMCP_Energy);
		m_pTTree1->Branch("lostChargedMCP_Pt",&m_lostChargedMCP_Pt);
		m_pTTree1->Branch("SLDecayXi", &m_SLDecayXi);
		m_pTTree1->Branch("SLDecayYi", &m_SLDecayYi);
		m_pTTree1->Branch("SLDecayZi", &m_SLDecayZi);
		m_pTTree1->Branch("SLDecayRi", &m_SLDecayRi);
		m_pTTree1->Branch("SLDecayXf", &m_SLDecayXf);
		m_pTTree1->Branch("SLDecayYf", &m_SLDecayYf);
		m_pTTree1->Branch("SLDecayZf", &m_SLDecayZf);
		m_pTTree1->Branch("SLDecayRf", &m_SLDecayRf);
		m_pTTree1->Branch("trueNuPx", &m_trueNuPx);
		m_pTTree1->Branch("trueNuPy", &m_trueNuPy);
		m_pTTree1->Branch("trueNuPz", &m_trueNuPz);
		m_pTTree1->Branch("trueNuE", &m_trueNuE);
		m_pTTree1->Branch("recoNuCloseInitialPx", &m_recoNuCloseInitialPx);
		m_pTTree1->Branch("recoNuCloseInitialPy", &m_recoNuCloseInitialPy);
		m_pTTree1->Branch("recoNuCloseInitialPz", &m_recoNuCloseInitialPz);
		m_pTTree1->Branch("recoNuCloseInitialE", &m_recoNuCloseInitialE);
		m_pTTree1->Branch("recoNuClosePx", &m_recoNuClosePx);
		m_pTTree1->Branch("recoNuClosePy", &m_recoNuClosePy);
		m_pTTree1->Branch("recoNuClosePz", &m_recoNuClosePz);
		m_pTTree1->Branch("recoNuCloseE", &m_recoNuCloseE);
		m_pTTree1->Branch("recoNuPosPx", &m_recoNuPosPx);
		m_pTTree1->Branch("recoNuPosPy", &m_recoNuPosPy);
		m_pTTree1->Branch("recoNuPosPz", &m_recoNuPosPz);
		m_pTTree1->Branch("recoNuPosE", &m_recoNuPosE);
		m_pTTree1->Branch("recoNuNegPx", &m_recoNuNegPx);
		m_pTTree1->Branch("recoNuNegPy", &m_recoNuNegPy);
		m_pTTree1->Branch("recoNuNegPz", &m_recoNuNegPz);
		m_pTTree1->Branch("recoNuNegE", &m_recoNuNegE);
		m_pTTree1->Branch("NuPxResidual", &m_NuPxResidual);
		m_pTTree1->Branch("NuPyResidual", &m_NuPyResidual);
		m_pTTree1->Branch("NuPzResidual", &m_NuPzResidual);
		m_pTTree1->Branch("NuEResidual", &m_NuEResidual);
		m_pTTree1->Branch("NuPxNormalizedResidual", &m_NuPxNormalizedResidual);
		m_pTTree1->Branch("NuPyNormalizedResidual", &m_NuPyNormalizedResidual);
		m_pTTree1->Branch("NuPzNormalizedResidual", &m_NuPzNormalizedResidual);
		m_pTTree1->Branch("NuENormalizedResidual", &m_NuENormalizedResidual);
		m_pTTree1->Branch("solutionSign", &m_solutionSign);
		m_pTTree1->Branch("E_vis", &m_E_vis);
		m_pTTree1->Branch("E_vis_prime", &m_E_vis_prime);
		m_pTTree1->Branch("P_vis_par", &m_P_vis_par);
		m_pTTree1->Branch("P_vis_par_prime", &m_P_vis_par_prime);
		m_pTTree1->Branch("P_vis_nor", &m_P_vis_nor);
		m_pTTree1->Branch("P_vis_nor_prime", &m_P_vis_nor_prime);
		m_pTTree1->Branch("flightDirectionStatus", &m_flightDirectionStatus);
		m_pTTree1->Branch("flightDirectionErrorCosAlpha", &m_FlightDirectionErrorCosAlpha);
		m_pTTree1->Branch("flightDirectionErrorSinAlpha", &m_FlightDirectionErrorSinAlpha);
		m_pTTree1->Branch("flightDirectionErrorAlpha", &m_FlightDirectionErrorAlpha);
		m_pTTree1->Branch("distRecoLeptonToDownStreamVertex", &m_distRecoLeptonToDownStreamVertex);
		m_pTTree1->Branch("dsVertexResidualX", &m_dsVertexResidualX);
		m_pTTree1->Branch("dsVertexResidualY", &m_dsVertexResidualY);
		m_pTTree1->Branch("dsVertexResidualZ", &m_dsVertexResidualZ);
		m_pTTree1->Branch("secVertexResidualX", &m_SecVertexResidualX);
		m_pTTree1->Branch("secVertexResidualY", &m_SecVertexResidualY);
		m_pTTree1->Branch("secVertexResidualZ", &m_SecVertexResidualZ);
		m_pTTree1->Branch("parentHadronMass", &m_parentHadronMass );
		m_pTTree1->Branch("parentHadronFlightDistance", &m_parentHadronFlightDistance );
		m_pTTree1->Branch("daughterHadronMass", &m_daughterHadronMass );
		m_pTTree1->Branch("SLDStatus", &m_SLDStatus );
		m_pTTree1->Branch("nTrueNeutralDecayProducts",&m_nTrueNeutralDecayProducts);
		m_pTTree1->Branch("nTrueAloneChargedDecayProducts",&m_nTrueAloneChargedDecayProducts);
		m_pTTree1->Branch("nTrueVertices",&m_nTrueVertices);
		m_pTTree1->Branch("nJetsVerticesDistributedIn",&m_nJetsVerticesDistributedIn);
		m_pTTree1->Branch("nRecoVerticesInJet",&m_nRecoVerticesInJet);
		m_pTTree1->Branch("nAloneChargedPFOs",&m_nAloneChargedPFOs);
		m_pTTree1->Branch("nAloneChargedPFOsFromSLD",&m_nAloneChargedPFOsFromSLD);
		m_pTTree1->Branch("distLeptonAlonePFOsNotFromSLD",&m_distLeptonAlonePFOsNotFromSLD);
		m_pTTree1->Branch("distLeptonAlonePFOsFromSLD",&m_distLeptonAlonePFOsFromSLD);
		m_pTTree1->Branch("widestConeAlphaNeutrals", &m_widestConeAlphaNeutrals );
		m_pTTree1->Branch("widestConeAlphaVertices", &m_widestConeAlphaVertices );
		m_pTTree1->Branch("widestConeAlphaAloneCharged", &m_widestConeAlphaCharged );
		m_pTTree1->Branch("widestConeCosAlphaNeutrals", &m_widestConeCosAlphaNeutrals );
		m_pTTree1->Branch("widestConeCosAlphaVertices", &m_widestConeCosAlphaVertices );
		m_pTTree1->Branch("widestConeCosAlphaAloneCharged", &m_widestConeCosAlphaCharged );
		m_pTTree1->Branch("widestConeAlphaAlonePFOsFromSLDwrtLepton", &m_widestConeAlphaAlonePFOsFromSLDwrtLepton );
		m_pTTree1->Branch("widestConeCosAlphaAlonePFOsFromSLDwrtLepton", &m_widestConeCosAlphaAlonePFOsFromSLDwrtLepton );
		m_pTTree1->Branch("widestConeAlphaAlonePFOsFromSLDwrtFD", &m_widestConeAlphaAlonePFOsFromSLDwrtFD );
		m_pTTree1->Branch("widestConeCosAlphaAlonePFOsFromSLDwrtFD", &m_widestConeCosAlphaAlonePFOsFromSLDwrtFD );
		m_pTTree1->Branch("widestConeAlphaAlonePFOsFromSLDwrtJet", &m_widestConeAlphaAlonePFOsFromSLDwrtJet );
		m_pTTree1->Branch("widestConeCosAlphaAlonePFOsFromSLDwrtJet", &m_widestConeCosAlphaAlonePFOsFromSLDwrtJet );
		m_pTTree1->Branch("widestConeAlphaAlonePFOsNotFromSLDwrtLepton", &m_widestConeAlphaAlonePFOsNotFromSLDwrtLepton );
		m_pTTree1->Branch("widestConeCosAlphaAlonePFOsNotFromSLDwrtLepton", &m_widestConeCosAlphaAlonePFOsNotFromSLDwrtLepton );
		m_pTTree1->Branch("widestConeAlphaAlonePFOsNotFromSLDwrtFD", &m_widestConeAlphaAlonePFOsNotFromSLDwrtFD );
		m_pTTree1->Branch("widestConeCosAlphaAlonePFOsNotFromSLDwrtFD", &m_widestConeCosAlphaAlonePFOsNotFromSLDwrtFD );
		m_pTTree1->Branch("widestConeAlphaAlonePFOsNotFromSLDwrtJet", &m_widestConeAlphaAlonePFOsNotFromSLDwrtJet );
		m_pTTree1->Branch("widestConeCosAlphaAlonePFOsNotFromSLDwrtJet", &m_widestConeCosAlphaAlonePFOsNotFromSLDwrtJet );
		m_pTTree1->Branch("weightPFOtoMCP_Lepton", &m_weightPFOtoMCP_Lepton );
		m_pTTree1->Branch("weightMCPtoPFO_Lepton", &m_weightMCPtoPFO_Lepton );
		m_pTTree1->Branch("weightPFOtoMCP_Neutral", &m_weightPFOtoMCP_Neutral );
		m_pTTree1->Branch("weightMCPtoPFO_Neutral", &m_weightMCPtoPFO_Neutral );
		m_pTTree1->Branch("weightPFOtoMCP_Charged", &m_weightPFOtoMCP_Charged );
		m_pTTree1->Branch("weightMCPtoPFO_Charged", &m_weightMCPtoPFO_Charged );

		m_pTTree2 = new TTree("FourMomentums", "FourMomentums");
		m_pTTree2->Branch("SLDStatus", &m_SLDStatus );
		m_pTTree2->Branch("trueNeutralPx", &m_trueNeutralPx );
		m_pTTree2->Branch("trueNeutralPy", &m_trueNeutralPy );
		m_pTTree2->Branch("trueNeutralPz", &m_trueNeutralPz );
		m_pTTree2->Branch("trueNeutralE", &m_trueNeutralE );
		m_pTTree2->Branch("trueChargedPx", &m_trueChargedPx );
		m_pTTree2->Branch("trueChargedPy", &m_trueChargedPy );
		m_pTTree2->Branch("trueChargedPz", &m_trueChargedPz );
		m_pTTree2->Branch("trueChargedE", &m_trueChargedE );
		m_pTTree2->Branch("trueLeptonPx", &m_trueLeptonPx );
		m_pTTree2->Branch("trueLeptonPy", &m_trueLeptonPy );
		m_pTTree2->Branch("trueLeptonPz", &m_trueLeptonPz );
		m_pTTree2->Branch("trueLeptonE", &m_trueLeptonE );
		m_pTTree2->Branch("trueVisiblePx", &m_trueVisiblePx );
		m_pTTree2->Branch("trueVisiblePy", &m_trueVisiblePy );
		m_pTTree2->Branch("trueVisiblePz", &m_trueVisiblePz );
		m_pTTree2->Branch("trueVisibleE", &m_trueVisibleE );
		m_pTTree2->Branch("trueNeutrinoPx", &m_trueNeutrinoPx );
		m_pTTree2->Branch("trueNeutrinoPy", &m_trueNeutrinoPy );
		m_pTTree2->Branch("trueNeutrinoPz", &m_trueNeutrinoPz );
		m_pTTree2->Branch("trueNeutrinoE", &m_trueNeutrinoE );
		m_pTTree2->Branch("trueHadronPx", &m_trueHadronPx );
		m_pTTree2->Branch("trueHadronPy", &m_trueHadronPy );
		m_pTTree2->Branch("trueHadronPz", &m_trueHadronPz );
		m_pTTree2->Branch("trueHadronE", &m_trueHadronE );
		m_pTTree2->Branch("recoNeutralPx", &m_recoNeutralPx );
		m_pTTree2->Branch("recoNeutralPy", &m_recoNeutralPy );
		m_pTTree2->Branch("recoNeutralPz", &m_recoNeutralPz );
		m_pTTree2->Branch("recoNeutralE", &m_recoNeutralE );
		m_pTTree2->Branch("recoChargedPx", &m_recoChargedPx );
		m_pTTree2->Branch("recoChargedPy", &m_recoChargedPy );
		m_pTTree2->Branch("recoChargedPz", &m_recoChargedPz );
		m_pTTree2->Branch("recoChargedE", &m_recoChargedE );
		m_pTTree2->Branch("recoLeptonPx", &m_recoLeptonPx );
		m_pTTree2->Branch("recoLeptonPy", &m_recoLeptonPy );
		m_pTTree2->Branch("recoLeptonPz", &m_recoLeptonPz );
		m_pTTree2->Branch("recoLeptonE", &m_recoLeptonE );
		m_pTTree2->Branch("recoVisiblePx", &m_recoVisiblePx );
		m_pTTree2->Branch("recoVisiblePy", &m_recoVisiblePy );
		m_pTTree2->Branch("recoVisiblePz", &m_recoVisiblePz );
		m_pTTree2->Branch("recoVisibleE", &m_recoVisibleE );
		m_pTTree2->Branch("recoNeutrinoPx", &m_recoNeutrinoPx );
		m_pTTree2->Branch("recoNeutrinoPy", &m_recoNeutrinoPy );
		m_pTTree2->Branch("recoNeutrinoPz", &m_recoNeutrinoPz );
		m_pTTree2->Branch("recoNeutrinoE", &m_recoNeutrinoE );
		m_pTTree2->Branch("recoHadronPx", &m_recoHadronPx );
		m_pTTree2->Branch("recoHadronPy", &m_recoHadronPy );
		m_pTTree2->Branch("recoHadronPz", &m_recoHadronPz );
		m_pTTree2->Branch("recoHadronE", &m_recoHadronE );



		h_SLDStatus = new TH1F( "SLDStatus" , ";" , 7 , 0.5 , 7.5 ); n_SLDStatus = 0;
		h_SLDStatus->GetXaxis()->SetBinLabel(1,"lep not found");
		h_SLDStatus->GetXaxis()->SetBinLabel(2,"lep not in jet");
		h_SLDStatus->GetXaxis()->SetBinLabel(3,"lep in Prim. Vtx");
		h_SLDStatus->GetXaxis()->SetBinLabel(4,"lep in Sec. Vtx");
		h_SLDStatus->GetXaxis()->SetBinLabel(5,"lep + 3^{rd} Vtx");
		h_SLDStatus->GetXaxis()->SetBinLabel(6,"lep + alone track");
		h_SLDStatus->GetXaxis()->SetBinLabel(7,"other");

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
		h_secondaryVertex = new TH1F( "secondary_vertices" , ";" , 11 , 0 , 11 ); n_secondaryVertex = 0;
		h_secondaryVertex->GetXaxis()->SetBinLabel(1,"reco lep not found");
		h_secondaryVertex->GetXaxis()->SetBinLabel(2,"reco lep not in jet");
		h_secondaryVertex->GetXaxis()->SetBinLabel(3,"lep in Sec. Vtx");
		h_secondaryVertex->GetXaxis()->SetBinLabel(4,"Ch. 3^{rd} Vtx + lep");
		h_secondaryVertex->GetXaxis()->SetBinLabel(5,"N. 3^{rd} Vtx + lep");
		h_secondaryVertex->GetXaxis()->SetBinLabel(6,"jet axis");
		h_secondaryVertex->GetXaxis()->SetBinLabel(7,"lead. par. in jet (Ch.)");
		h_secondaryVertex->GetXaxis()->SetBinLabel(8,"lead. par. in jet (N.)");
		h_secondaryVertex->GetXaxis()->SetBinLabel(9,"lead. par. in jet (#gamma)");
		h_secondaryVertex->GetXaxis()->SetBinLabel(10,"SLD-lep");
		h_secondaryVertex->GetXaxis()->SetBinLabel(11,"lep in Prim. Vtx");
	//	h_secondaryVertex->GetXaxis()->SetBinLabel(6,"other(?)");
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
	m_nTrueNeutralDecayProducts.clear();
	m_nTrueAloneChargedDecayProducts.clear();
	m_nTrueVertices.clear();
	m_nJetsVerticesDistributedIn.clear();
	m_nRecoVerticesInJet.clear();
	m_nAloneChargedPFOs.clear();
	m_nAloneChargedPFOsFromSLD.clear();
	m_distLeptonAlonePFOsNotFromSLD.clear();
	m_distLeptonAlonePFOsFromSLD.clear();
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
	m_solutionSign.clear();
	m_E_vis.clear();
	m_E_vis_prime.clear();
	m_P_vis_par.clear();
	m_P_vis_par_prime.clear();
	m_P_vis_nor.clear();
	m_P_vis_nor_prime.clear();
	m_flightDirectionStatus.clear();
	m_FlightDirectionErrorCosAlpha.clear();
	m_FlightDirectionErrorSinAlpha.clear();
	m_FlightDirectionErrorAlpha.clear();
	m_dsVertexResidualX.clear();
	m_dsVertexResidualY.clear();
	m_dsVertexResidualZ.clear();
	m_SecVertexResidualX.clear();
	m_SecVertexResidualY.clear();
	m_SecVertexResidualZ.clear();
	m_parentHadronMass.clear();
	m_parentHadronFlightDistance.clear();
	m_daughterHadronMass.clear();
	m_SLDStatus.clear();
	m_widestConeAlphaNeutrals.clear();
	m_widestConeAlphaVertices.clear();
	m_widestConeAlphaCharged.clear();
	m_widestConeCosAlphaNeutrals.clear();
	m_widestConeCosAlphaVertices.clear();
	m_widestConeCosAlphaCharged.clear();
	m_widestConeAlphaAlonePFOsFromSLDwrtLepton.clear();
	m_widestConeCosAlphaAlonePFOsFromSLDwrtLepton.clear();
	m_widestConeAlphaAlonePFOsFromSLDwrtFD.clear();
	m_widestConeCosAlphaAlonePFOsFromSLDwrtFD.clear();
	m_widestConeAlphaAlonePFOsFromSLDwrtJet.clear();
	m_widestConeCosAlphaAlonePFOsFromSLDwrtJet.clear();
	m_widestConeAlphaAlonePFOsNotFromSLDwrtLepton.clear();
	m_widestConeCosAlphaAlonePFOsNotFromSLDwrtLepton.clear();
	m_widestConeAlphaAlonePFOsNotFromSLDwrtFD.clear();
	m_widestConeCosAlphaAlonePFOsNotFromSLDwrtFD.clear();
	m_widestConeAlphaAlonePFOsNotFromSLDwrtJet.clear();
	m_widestConeCosAlphaAlonePFOsNotFromSLDwrtJet.clear();
	m_weightPFOtoMCP_Lepton.clear();
	m_weightMCPtoPFO_Lepton.clear();
	m_weightPFOtoMCP_Neutral.clear();
	m_weightMCPtoPFO_Neutral.clear();
	m_weightPFOtoMCP_Charged.clear();
	m_weightMCPtoPFO_Charged.clear();
	m_distRecoLeptonToDownStreamVertex.clear();

	m_trueNeutralPx.clear();
	m_trueNeutralPy.clear();
	m_trueNeutralPz.clear();
	m_trueNeutralE.clear();
	m_trueChargedPx.clear();
	m_trueChargedPy.clear();
	m_trueChargedPz.clear();
	m_trueChargedE.clear();
	m_trueLeptonPx.clear();
	m_trueLeptonPy.clear();
	m_trueLeptonPz.clear();
	m_trueLeptonE.clear();
	m_trueVisiblePx.clear();
	m_trueVisiblePy.clear();
	m_trueVisiblePz.clear();
	m_trueVisibleE.clear();
	m_trueNeutrinoPx.clear();
	m_trueNeutrinoPy.clear();
	m_trueNeutrinoPz.clear();
	m_trueNeutrinoE.clear();
	m_trueHadronPx.clear();
	m_trueHadronPy.clear();
	m_trueHadronPz.clear();
	m_trueHadronE.clear();
	m_recoNeutralPx.clear();
	m_recoNeutralPy.clear();
	m_recoNeutralPz.clear();
	m_recoNeutralE.clear();
	m_recoChargedPx.clear();
	m_recoChargedPy.clear();
	m_recoChargedPz.clear();
	m_recoChargedE.clear();
	m_recoLeptonPx.clear();
	m_recoLeptonPy.clear();
	m_recoLeptonPz.clear();
	m_recoLeptonE.clear();
	m_recoVisiblePx.clear();
	m_recoVisiblePy.clear();
	m_recoVisiblePz.clear();
	m_recoVisibleE.clear();
	m_recoNeutrinoPx.clear();
	m_recoNeutrinoPy.clear();
	m_recoNeutrinoPz.clear();
	m_recoNeutrinoE.clear();
	m_recoHadronPx.clear();
	m_recoHadronPy.clear();
	m_recoHadronPz.clear();
	m_recoHadronE.clear();

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
//					if ( abs( testLepton->getPDG() ) != 11 ) continue;
					streamlog_out(DEBUG3) << "" << std::endl;
					streamlog_out(DEBUG3) << "	<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
					streamlog_out(DEBUG3) << "	<<<<<<<<<<<<<<<<<<<<<<<<<<< Found a primary semi-leptonic decay >>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
					if ( downStreamSLDecay )
					{
						streamlog_out(DEBUG3) << "	There is/are downstream semi-leptonic(s) decay in primary semi-leptonic decay products" << std::endl;
						if ( m_fillRootTree ) h_SLDecayOrder->Fill( 2.5 );
					}
					if ( upStreamSLDecay )
					{
						streamlog_out(DEBUG3) << "	There is/are upstream semi-leptonic(s) decay in primary semi-leptonic decay products" << std::endl;
						if ( m_fillRootTree ) h_SLDecayOrder->Fill( 0.5 );
					}
					if ( !downStreamSLDecay && !upStreamSLDecay )
					{
						streamlog_out(DEBUG3) << "	<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
						streamlog_out(DEBUG3) << "	<<<<<<<<<<<<<<<< There are no upstream and downstream semi-leptonic decay >>>>>>>>>>>>>>>>>" << std::endl;
						streamlog_out(DEBUG3) << "	<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
						if ( m_fillRootTree ) h_SLDecayOrder->Fill( 1.5 );
						m_jetFlavourPDG.push_back( jetFlavourPDG );
						doSLDCorrection( pLCEvent , testLepton );
						m_parentHadronMass.push_back( ( testLepton->getParents()[ 0 ] )->getMass() );
						m_parentHadronFlightDistance.push_back( std::sqrt( std::pow( testLepton->getParents()[ 0 ]->getEndpoint()[ 0 ] - testLepton->getParents()[ 0 ]->getVertex()[ 0 ] , 2 ) + std::pow( testLepton->getParents()[ 0 ]->getEndpoint()[ 1 ] - testLepton->getParents()[ 0 ]->getVertex()[ 1 ] , 2 ) + std::pow( testLepton->getParents()[ 0 ]->getEndpoint()[ 2 ] - testLepton->getParents()[ 0 ]->getVertex()[ 2 ] , 2 ) ) );
						for ( unsigned int i_d = 0 ; i_d < ( testLepton->getParents()[ 0 ] )->getDaughters().size() ; ++i_d )
						{
							MCParticle *mcDaughter = ( testLepton->getParents()[ 0 ] )->getDaughters()[ i_d ];
							int daughterPDG = std::abs( mcDaughter->getPDG() );
							if ( daughterPDG < 11 || daughterPDG > 16 ) m_daughterHadronMass.push_back( mcDaughter->getMass() );
						}
					}
				}
			}
		}
		m_nTauNeutrino = nTauNeutrino;
		if ( m_fillRootTree )
		{
			m_pTTree1->Fill();
			m_pTTree2->Fill();
		}

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

int SLDCorrection::getVertexInJetsDistribution( Vertex* testVertex , pfoVector jetVector )
{
	pfoVector jets{};
	jets.clear();
	ReconstructedParticle* vertexParticle = testVertex->getAssociatedParticle();
	for ( unsigned int i_par = 0 ; i_par < vertexParticle->getParticles().size() ; ++i_par )
	{
		ReconstructedParticle* testParticle = vertexParticle->getParticles()[ i_par ];
		for ( unsigned int i_jet = 0 ; i_jet < jetVector.size() ; ++i_jet )
		{
			ReconstructedParticle* jet = jetVector[ i_jet ];
			for ( unsigned int i_jetPar = 0 ; i_jetPar < jet->getParticles().size() ; ++i_jetPar )
			{
				if ( jet->getParticles()[ i_jetPar ] == testParticle )
				{
					bool jetWasInList = false;
					for ( unsigned int j = 0 ; j < jets.size() ; ++j )
					{
						if ( jets[ j ] == jet ) jetWasInList = true;
					}
					if ( !jetWasInList ) jets.push_back( jet );
				}
			}
		}
	}
	return jets.size();
}

void SLDCorrection::doSLDCorrection( EVENT::LCEvent *pLCEvent , MCParticle *SLDLepton )
{
	++n_SLDStatus;
	ReconstructedParticle *assignedJet = NULL;
	LCCollection *primaryVertexCollection = pLCEvent->getCollection( m_inputPrimaryVertex );
	Vertex* primaryVertex = dynamic_cast<Vertex*>( primaryVertexCollection->getElementAt( 0 ) );
	Vertex* startVertex = dynamic_cast<Vertex*>( primaryVertexCollection->getElementAt( 0 ) );
	LCCollection *jetCollection = pLCEvent->getCollection( m_inputJetCollection );
	pfoVector jetVector{};
	for ( int i_jet = 0 ; i_jet < jetCollection->getNumberOfElements(); ++i_jet)
	{
		ReconstructedParticle* jet = dynamic_cast<ReconstructedParticle*>( jetCollection->getElementAt( i_jet ) );
		jetVector.push_back( jet );
	}
	LCCollection *buildUpVertexCollection = pLCEvent->getCollection( m_inputBuildUpVertex );
	vtxVector buildUpVertexVector{};
	for ( int i_vtx = 0 ; i_vtx < buildUpVertexCollection->getNumberOfElements(); ++i_vtx)
	{
		Vertex* vertex = dynamic_cast<Vertex*>( buildUpVertexCollection->getElementAt( i_vtx ) );
		buildUpVertexVector.push_back( vertex );
	}
	bool recoLeptonIsInJet = false;

	if ( m_displayEvent )
	{
		DDMarlinCED::newEvent( this ); // refresh
		DDMarlinCED::drawDD4hepDetector( this->_theDetector , 0 , std::vector<std::string>{} ); // draw geometry
		DDCEDPickingHandler& pHandler = DDCEDPickingHandler::getInstance();
		pHandler.update(pLCEvent);
		drawMCParticles( SLDLepton->getParents()[ 0 ] );
	}
	LCRelationNavigator RecoMCParticleNav( pLCEvent->getCollection( m_RecoMCTruthLinkCollection ) );
	LCRelationNavigator MCParticleRecoNav( pLCEvent->getCollection( m_MCTruthRecoLinkCollection ) );


	TLorentzVector trueLeptonFourMomentum( SLDLepton->getMomentum() , SLDLepton->getEnergy() );
	TLorentzVector trueChargedFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector trueNeutralFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector trueVisibleFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector trueNeutrinoFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	MCParticle* trueNeutrino = getTrueNeutrino( SLDLepton , trueNeutrinoFourMomentum );
	TLorentzVector trueHadronFourMomentum( ( SLDLepton->getParents()[ 0 ] )->getMomentum() , ( SLDLepton->getParents()[ 0 ] )->getEnergy() );

	TLorentzVector recoLeptonFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector recoChargedFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector recoNeutralFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector recoVisibleFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector recoNeutrinoFourMomentumPos( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector recoNeutrinoFourMomentumNeg( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector recoNeutrinoFourMomentumClose( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector recoHadronFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	std::vector< float > NeutrinoCovMat( 10 , 0.0 );

	TLorentzVector leptonFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector chargedFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector neutralFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector visibleFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );


	TVector3 trueFlightDirection( 0.0 , 0.0 , 0.0 );
	TVector3 recoFlightDirection( 0.0 , 0.0 , 0.0 );
	TVector3 flightDirection( 0.0 , 0.0 , 0.0 );
	double parentHadronMass = 0.0;
	float helicesDistance = 0.0;

	std::vector<double> SecondaryVertexPar;

	MCParticle *parentHadron = SLDLepton->getParents()[ 0 ];
	parentHadronMass = ( SLDLepton->getParents()[ 0 ] )->getMass();
	mcpVector trueNeutralDecayProducts{};
	mcpVector trueChargedDecayProductsAll{};
	mcpVector trueChargedDecayProducts{};
	mcpVector aloneChargedDecayProducts{};
	mcpVector MCParticlesWithVertex{};
	pfoVector PFOswithAloneTracks{};
	pfoVector PFOswithAloneTracksFromSLD{};
	pfoVector PFOswithAloneTracksNotFromSLD{};
	vtxVector SLDVertices{};
	pfoVector SLDVerticesRP{};

	std::vector<float> truePrimaryVertex{};
	truePrimaryVertex.push_back( parentHadron->getVertex()[ 0 ] );
	truePrimaryVertex.push_back( parentHadron->getVertex()[ 1 ] );
	truePrimaryVertex.push_back( parentHadron->getVertex()[ 2 ] );

	streamlog_out(DEBUG2) << "" << std::endl;
	streamlog_out(DEBUG2) << "----------------------------------------------------------------------" << std::endl;
	streamlog_out(DEBUG2) << "--------------------- Stable Neutral MCParticles ---------------------" << std::endl;
	streamlog_out(DEBUG2) << "----------------------------------------------------------------------" << std::endl;
	int nTrueNeutralMCPs = getNeutralMCPs( parentHadron , trueNeutralDecayProducts );
	streamlog_out(DEBUG2) << "----------------------------------------------------------------------" << std::endl;
	streamlog_out(DEBUG2) << "----------------- " << nTrueNeutralMCPs << " Stable Neutral MCParticles Found------------------" << std::endl;
	streamlog_out(DEBUG2) << "----------------------------------------------------------------------" << std::endl;
	streamlog_out(DEBUG2) << "" << std::endl;
	streamlog_out(DEBUG2) << "----------------------------------------------------------------------" << std::endl;
	streamlog_out(DEBUG2) << "--------------------- Stable Charged MCParticles ---------------------" << std::endl;
	streamlog_out(DEBUG2) << "----------------------------------------------------------------------" << std::endl;
	int nTrueChargedMCPs = getChargedMCPs( SLDLepton , parentHadron , trueChargedDecayProductsAll );
	streamlog_out(DEBUG2) << "----------------------------------------------------------------------" << std::endl;
	streamlog_out(DEBUG2) << "----------------- " << nTrueChargedMCPs << " Stable Charged MCParticles Found------------------" << std::endl;
	streamlog_out(DEBUG2) << "----------------------------------------------------------------------" << std::endl;
	streamlog_out(DEBUG2) << "" << std::endl;
	streamlog_out(DEBUG2) << "----------------------------------------------------------------------" << std::endl;
	streamlog_out(DEBUG2) << "-------------------- Vertex of Stable MCParticles --------------------" << std::endl;
	streamlog_out(DEBUG2) << "----------------------------------------------------------------------" << std::endl;
	int nTrueVertices = getTrueVertices( SLDLepton->getParents()[ 0 ] , MCParticlesWithVertex );
	streamlog_out(DEBUG2) << "---------------- " << nTrueVertices << " Vertex of Stable MCParticles Found-----------------" << std::endl;
	m_nTrueNeutralDecayProducts.push_back( nTrueNeutralMCPs );
	m_nTrueVertices.push_back( nTrueVertices );
	for ( unsigned int i_chMCP = 0 ; i_chMCP < trueChargedDecayProductsAll.size() ; ++i_chMCP )
	{
		bool mcpIsInTrueVertex = false;
		bool mcpIsSLDLepton = false;
		MCParticle* ChargedMCPCandidate = trueChargedDecayProductsAll[ i_chMCP ];
		if ( ChargedMCPCandidate == SLDLepton ) mcpIsSLDLepton = true;
		for ( unsigned int i_vtx = 0 ; i_vtx < MCParticlesWithVertex.size() ; ++i_vtx )
		{
			MCParticle* mcpWithVertex = MCParticlesWithVertex[ i_vtx ];
			for ( unsigned int i_d = 0 ; i_d < mcpWithVertex->getDaughters().size() ; ++i_d )
			{
				MCParticle* daughter = mcpWithVertex->getDaughters()[ i_d ];
				if ( daughter == ChargedMCPCandidate ) mcpIsInTrueVertex = true;
			}
		}
		if ( mcpIsInTrueVertex || mcpIsSLDLepton )
		{
			trueChargedDecayProducts.push_back( ChargedMCPCandidate );
		}
		else
		{
			aloneChargedDecayProducts.push_back( ChargedMCPCandidate );
		}
	}
	m_nTrueAloneChargedDecayProducts.push_back( aloneChargedDecayProducts.size() );

	streamlog_out(DEBUG2) << "	In Total, " << nTrueChargedMCPs << " charged MCParticles and " << nTrueNeutralMCPs << " neutral MCParticles (except neutrino) and " << nTrueVertices << " Vertices are assigned to Semi-Leptonic Decay" << std::endl;

	for ( unsigned int i_mcp = 0 ; i_mcp < trueNeutralDecayProducts.size() ; ++i_mcp )
	{
		trueNeutralFourMomentum += TLorentzVector( trueNeutralDecayProducts[ i_mcp ]->getMomentum() , trueNeutralDecayProducts[ i_mcp ]->getEnergy() );
	}
	for ( unsigned int i_mcp = 0 ; i_mcp < trueChargedDecayProductsAll.size() ; ++i_mcp )
	{
		trueChargedFourMomentum += TLorentzVector( trueChargedDecayProductsAll[ i_mcp ]->getMomentum() , trueChargedDecayProductsAll[ i_mcp ]->getEnergy() );
	}

	std::vector<double> trueStartVertex{};
	std::vector<double> trueSLDVertex{};
	getTrueFlightDirection( SLDLepton , trueFlightDirection , trueStartVertex , trueSLDVertex );
	if ( MCParticlesWithVertex.size() != 0 )
	{
		m_widestConeAlphaVertices.push_back( acos( getWidestCosAlphaOfVertices( MCParticlesWithVertex , trueFlightDirection , truePrimaryVertex ) ) * 180.0 / 3.14159265 );
		m_widestConeCosAlphaVertices.push_back( getWidestCosAlphaOfVertices( MCParticlesWithVertex , trueFlightDirection , truePrimaryVertex ) );
	}
	if ( trueNeutralDecayProducts.size() != 0 )
	{
		m_widestConeAlphaNeutrals.push_back( acos( getWidestCosAlphaOfDecayProducts( trueNeutralDecayProducts , trueFlightDirection ) ) * 180.0 / 3.14159265 );
		m_widestConeCosAlphaNeutrals.push_back( getWidestCosAlphaOfDecayProducts( trueNeutralDecayProducts , trueFlightDirection ) );
	}
	if ( aloneChargedDecayProducts.size() != 0 )
	{
		m_widestConeAlphaCharged.push_back( acos( getWidestCosAlphaOfDecayProducts( aloneChargedDecayProducts , trueFlightDirection ) ) * 180.0 / 3.14159265 );
		m_widestConeCosAlphaCharged.push_back( getWidestCosAlphaOfDecayProducts( aloneChargedDecayProducts , trueFlightDirection ) );
	}

	float weightPFOtoMCP = 0.0;
	float weightMCPtoPFO = 0.0;

	for ( unsigned int i_vtx = 0 ; i_vtx < buildUpVertexVector.size() ; ++i_vtx )
	{
		m_nJetsVerticesDistributedIn.push_back( getVertexInJetsDistribution( buildUpVertexVector[ i_vtx ] , jetVector ) );
	}


	ReconstructedParticle* linkedRecoLepton = getLinkedPFO( SLDLepton , RecoMCParticleNav , MCParticleRecoNav , true , false , weightPFOtoMCP , weightMCPtoPFO );
	if ( linkedRecoLepton == NULL )
	{
		m_SLDStatus.push_back( 1 );
		streamlog_out(WARNING) << "	||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << std::endl;
		streamlog_out(WARNING) << "	||||||||||||||||| Reconstructed Lepton is not found ||||||||||||||||||" << std::endl;
		streamlog_out(WARNING) << "	||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << std::endl;
		if ( m_fillRootTree ) h_SLDStatus->Fill( 1 );
		return;
	}
	assignedJet = getJetAssignedToParticle( linkedRecoLepton , jetVector , recoLeptonIsInJet );
	if ( !recoLeptonIsInJet )
	{
		m_SLDStatus.push_back( 2 );
		streamlog_out(WARNING) << "	||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << std::endl;
		streamlog_out(WARNING) << "	||||||||||| Reconstructed Lepton doesn't belong to any jet |||||||||||" << std::endl;
		streamlog_out(WARNING) << "	||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << std::endl;
		if ( m_fillRootTree ) h_SLDStatus->Fill( 2 );
		return;
	}
	if ( isParticleInVertex( linkedRecoLepton , primaryVertex ) )
	{
		m_SLDStatus.push_back( 3 );
		streamlog_out(WARNING) << "	||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << std::endl;
		streamlog_out(WARNING) << "	||||||||||||| Reconstructed Lepton is in primary vertex ||||||||||||||" << std::endl;
		streamlog_out(WARNING) << "	||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << std::endl;
		if ( m_fillRootTree ) h_SLDStatus->Fill( 3 );
		return;
	}

	m_weightPFOtoMCP_Lepton.push_back( weightPFOtoMCP );
	m_weightMCPtoPFO_Lepton.push_back( weightMCPtoPFO );

	PFOswithAloneTracks = getParticlesWithAloneTracks( linkedRecoLepton , assignedJet , primaryVertex , buildUpVertexVector );
	for ( unsigned int i_pfo = 0 ; i_pfo < PFOswithAloneTracks.size() ; ++i_pfo )
	{
		weightPFOtoMCP = 0.0;
		weightMCPtoPFO = 0.0;
		bool PFOisFromSLD = false;
		MCParticle* linkedMCP = getLinkedMCP( PFOswithAloneTracks[ i_pfo ] , RecoMCParticleNav , MCParticleRecoNav , true , false , weightPFOtoMCP , weightMCPtoPFO );
		for ( unsigned int i_mcp = 0 ; i_mcp < trueChargedDecayProductsAll.size() ; ++i_mcp )
		{
			if ( linkedMCP == trueChargedDecayProductsAll[ i_mcp ] ) PFOisFromSLD = true;
		}
		if ( PFOisFromSLD )
		{
			PFOswithAloneTracksFromSLD.push_back( PFOswithAloneTracks[ i_pfo ] );
		}
		else
		{
			PFOswithAloneTracksNotFromSLD.push_back( PFOswithAloneTracks[ i_pfo ] );
		}
	}
	TVector3 leptonDirection = TVector3( SLDLepton->getMomentum() ); leptonDirection.SetMag( 1.0 );
	TVector3 jetAxis = TVector3( assignedJet->getMomentum() ); jetAxis.SetMag( 1.0 );

	if ( PFOswithAloneTracksFromSLD.size() != 0 )
	{
		for ( unsigned int i_pfo = 0 ; i_pfo < PFOswithAloneTracksFromSLD.size() ; ++i_pfo )
		{
			m_widestConeAlphaAlonePFOsFromSLDwrtLepton.push_back( acos( getWidestCosAlphaOfChargedPFOs( PFOswithAloneTracksFromSLD , leptonDirection ) ) );
			m_widestConeCosAlphaAlonePFOsFromSLDwrtLepton.push_back( getWidestCosAlphaOfChargedPFOs( PFOswithAloneTracksFromSLD , leptonDirection ) );
			m_widestConeAlphaAlonePFOsFromSLDwrtFD.push_back( acos( getWidestCosAlphaOfChargedPFOs( PFOswithAloneTracksFromSLD , trueFlightDirection ) ) );
			m_widestConeCosAlphaAlonePFOsFromSLDwrtFD.push_back( getWidestCosAlphaOfChargedPFOs( PFOswithAloneTracksFromSLD , trueFlightDirection ) );
			m_widestConeAlphaAlonePFOsFromSLDwrtJet.push_back( acos( getWidestCosAlphaOfChargedPFOs( PFOswithAloneTracksFromSLD , jetAxis ) ) );
			m_widestConeCosAlphaAlonePFOsFromSLDwrtJet.push_back( getWidestCosAlphaOfChargedPFOs( PFOswithAloneTracksFromSLD , jetAxis ) );
		}
	}
	if ( PFOswithAloneTracksNotFromSLD.size() != 0 )
	{
		for ( unsigned int i_pfo = 0 ; i_pfo < PFOswithAloneTracksNotFromSLD.size() ; ++i_pfo )
		{
			m_widestConeAlphaAlonePFOsNotFromSLDwrtLepton.push_back( acos( getWidestCosAlphaOfChargedPFOs( PFOswithAloneTracksNotFromSLD , leptonDirection ) ) );
			m_widestConeCosAlphaAlonePFOsNotFromSLDwrtLepton.push_back( getWidestCosAlphaOfChargedPFOs( PFOswithAloneTracksNotFromSLD , leptonDirection ) );
			m_widestConeAlphaAlonePFOsNotFromSLDwrtFD.push_back( acos( getWidestCosAlphaOfChargedPFOs( PFOswithAloneTracksNotFromSLD , trueFlightDirection ) ) );
			m_widestConeCosAlphaAlonePFOsNotFromSLDwrtFD.push_back( getWidestCosAlphaOfChargedPFOs( PFOswithAloneTracksNotFromSLD , trueFlightDirection ) );
			m_widestConeAlphaAlonePFOsNotFromSLDwrtJet.push_back( acos( getWidestCosAlphaOfChargedPFOs( PFOswithAloneTracksNotFromSLD , jetAxis ) ) );
			m_widestConeCosAlphaAlonePFOsNotFromSLDwrtJet.push_back( getWidestCosAlphaOfChargedPFOs( PFOswithAloneTracksNotFromSLD , jetAxis ) );
		}
	}
/*
	Track* leptonTrack = linkedRecoLepton->getTracks()[ 0 ];
	std::vector<double> PCAatTrack1;
	std::vector<double> PCAatTrack2;
	std::vector<double> PCAatTrack;
	std::vector<double> PCAatLine;
	for ( unsigned int i_pfo = 0 ; i_pfo < PFOswithAloneTracksNotFromSLD.size() ; ++i_pfo )
	{
		ReconstructedParticle* testPFO = PFOswithAloneTracksNotFromSLD[ i_pfo ];
		if ( testPFO->getTracks().size() == 1 )
		{
			PCAatTrack1.clear();
			PCAatTrack2.clear();
			m_distLeptonAlonePFOsNotFromSLD.push_back( intersectTrackTrack( leptonTrack , testPFO->getTracks()[ 0 ] , PCAatTrack1 , PCAatTrack2 ) );
		}
		else
		{
			PCAatTrack.clear();
			PCAatLine.clear();
			TVector3 momentumOfLine( testPFO->getMomentum() );
			PCAatTrack1.clear();
			PCAatTrack2.clear();
			Track* track1 = testPFO->getTracks()[ 0 ];
			Track* track2 = testPFO->getTracks()[ 1 ];
			intersectTrackTrack( track1 , track2 , PCAatTrack1 , PCAatTrack2 );
			std::vector<double> pointOnLine{};
			pointOnLine = ( std::fabs( track1->getOmega() ) < std::fabs( track2->getOmega() ) ? PCAatTrack1 : PCAatTrack2 );
			m_distLeptonAlonePFOsNotFromSLD.push_back( intersectTrackLine( leptonTrack , primaryVertex , momentumOfLine , pointOnLine , PCAatTrack , PCAatLine ) );
		}
	}
	for ( unsigned int i_pfo = 0 ; i_pfo < PFOswithAloneTracksFromSLD.size() ; ++i_pfo )
	{
		ReconstructedParticle* testPFO = PFOswithAloneTracksFromSLD[ i_pfo ];
		if ( testPFO->getTracks().size() == 1 )
		{
			PCAatTrack1.clear();
			PCAatTrack2.clear();
			m_distLeptonAlonePFOsFromSLD.push_back( intersectTrackTrack( leptonTrack , testPFO->getTracks()[ 0 ] , PCAatTrack1 , PCAatTrack2 ) );
		}
		else
		{
			PCAatTrack.clear();
			PCAatLine.clear();
			TVector3 momentumOfLine( testPFO->getMomentum() );
			PCAatTrack1.clear();
			PCAatTrack2.clear();
			Track* track1 = testPFO->getTracks()[ 0 ];
			Track* track2 = testPFO->getTracks()[ 1 ];
			intersectTrackTrack( track1 , track2 , PCAatTrack1 , PCAatTrack2 );
			std::vector<double> pointOnLine{};
			pointOnLine = ( std::fabs( track1->getOmega() ) < std::fabs( track2->getOmega() ) ? PCAatTrack1 : PCAatTrack2 );
			m_distLeptonAlonePFOsFromSLD.push_back( intersectTrackLine( leptonTrack , primaryVertex , momentumOfLine , pointOnLine , PCAatTrack , PCAatLine ) );
		}
	}
*/
	vtxVector verticesInJet = getVerticesInJet( assignedJet , buildUpVertexVector );
	m_nRecoVerticesInJet.push_back( verticesInJet.size() );
	m_nAloneChargedPFOs.push_back( PFOswithAloneTracks.size() );
	m_nAloneChargedPFOsFromSLD.push_back( PFOswithAloneTracksFromSLD.size() );
	int SLDStatus = getRecoFlightDirection( linkedRecoLepton , recoFlightDirection , primaryVertex , startVertex , SLDVertices , SLDVerticesRP , assignedJet , verticesInJet , PFOswithAloneTracks , helicesDistance );
	m_distRecoLeptonToDownStreamVertex.push_back( helicesDistance );

	m_SLDStatus.push_back( SLDStatus );
	if ( m_fillRootTree ) h_SLDStatus->Fill( SLDStatus );

	if ( m_cheatFlightDirection )
	{
		flightDirection = trueFlightDirection;
	}
	else
	{
		flightDirection = recoFlightDirection;
	}


	recoLeptonFourMomentum = TLorentzVector( linkedRecoLepton->getMomentum() , linkedRecoLepton->getEnergy() );
	if ( m_cheatLepton4momentum )
	{
		leptonFourMomentum = trueLeptonFourMomentum;
	}
	else
	{
		leptonFourMomentum = recoLeptonFourMomentum;
	}

	pfoVector recoChargedDecayProducts{};
	if ( m_cheatPIDcharged )
	{
		for ( unsigned int i_par = 0 ; i_par < trueChargedDecayProductsAll.size() ; ++i_par )
		{
			ReconstructedParticle* linkedPFO = getLinkedPFO( trueChargedDecayProductsAll[ i_par ] , RecoMCParticleNav , MCParticleRecoNav , true , false , weightPFOtoMCP , weightMCPtoPFO );
			if ( linkedPFO != NULL )
			{
				recoChargedDecayProducts.push_back( linkedPFO );
				if ( m_displayEvent ) drawReconstructedParticle( linkedPFO , primaryVertex , 0x0075df , 0x000000 );
			}
		}
	}
	else if ( SLDStatus == 4 )
	{
		ReconstructedParticle* sldVertexRP = SLDVerticesRP[ 0 ];
		for ( unsigned int i_par = 0 ; i_par < sldVertexRP->getParticles().size() ; ++i_par )
		{
			ReconstructedParticle* chargedDecayProduct = sldVertexRP->getParticles()[ i_par ];
			if ( chargedDecayProduct != linkedRecoLepton )
			{
				recoChargedDecayProducts.push_back( chargedDecayProduct );
				if ( m_displayEvent ) drawReconstructedParticle( chargedDecayProduct , primaryVertex , 0x0075df , 0x000000 );
			}
		}
	}
	else if ( SLDStatus == 5 )
	{
		ReconstructedParticle* sldVertexRP = SLDVerticesRP[ 0 ];
		for ( unsigned int i_par = 0 ; i_par < sldVertexRP->getParticles().size() ; ++i_par )
		{
			ReconstructedParticle* chargedDecayProduct = sldVertexRP->getParticles()[ i_par ];
			recoChargedDecayProducts.push_back( chargedDecayProduct );
			if ( m_displayEvent ) drawReconstructedParticle( chargedDecayProduct , primaryVertex , 0x0075df , 0x2e8e04 );
		}
	}

	for ( unsigned int i_par = 0 ; i_par < recoChargedDecayProducts.size() ; ++i_par )
	{
		recoChargedFourMomentum += TLorentzVector( recoChargedDecayProducts[ i_par ]->getMomentum() , recoChargedDecayProducts[ i_par ]->getEnergy() );
	}

	if ( m_cheatCharged4momentum )
	{
		chargedFourMomentum = trueChargedFourMomentum;
	}
	else
	{
		chargedFourMomentum = recoChargedFourMomentum;
	}

	pfoVector recoNeutralDecayProducts{};
	if ( m_cheatPIDneutral )
	{
		for ( unsigned int i_par = 0 ; i_par < trueNeutralDecayProducts.size() ; ++i_par )
		{
			ReconstructedParticle* linkedPFO = getLinkedPFO( trueNeutralDecayProducts[ i_par ] , RecoMCParticleNav , MCParticleRecoNav , false , true , weightPFOtoMCP , weightMCPtoPFO );
			if ( linkedPFO != NULL )
			{
				recoNeutralDecayProducts.push_back( linkedPFO );
				if ( m_displayEvent ) drawReconstructedParticle( linkedPFO , primaryVertex , 0x0075df , 0x2e8e04 );
			}
		}
	}
	else
	{
		streamlog_out(WARNING) << "Nothing TODO for associating neutral particles to semi-leptonic decay" << std::endl;
	}

	for ( unsigned int i_par = 0 ; i_par < recoNeutralDecayProducts.size() ; ++i_par )
	{
		recoNeutralFourMomentum += TLorentzVector( recoNeutralDecayProducts[ i_par ]->getMomentum() , recoNeutralDecayProducts[ i_par ]->getEnergy() );
	}

	if ( m_cheatNeutral4momentum )
	{
		neutralFourMomentum = trueNeutralFourMomentum;
	}
	else
	{
		neutralFourMomentum = recoNeutralFourMomentum;
	}
	trueVisibleFourMomentum = trueChargedFourMomentum + trueNeutralFourMomentum + trueLeptonFourMomentum;
	recoVisibleFourMomentum = recoChargedFourMomentum + recoNeutralFourMomentum + recoLeptonFourMomentum;
	visibleFourMomentum = chargedFourMomentum + neutralFourMomentum + leptonFourMomentum;

	showTrueParameters( SLDLepton );
	recoNeutrinoFourMomentumPos = getNeutrinoFourMomentum( flightDirection , visibleFourMomentum , parentHadronMass , +1.0 );
	recoNeutrinoFourMomentumNeg = getNeutrinoFourMomentum( flightDirection , visibleFourMomentum , parentHadronMass , -1.0 );

	recoNeutrinoFourMomentumClose = ( fabs( recoNeutrinoFourMomentumPos.E() - trueNeutrinoFourMomentum.E() ) < fabs( recoNeutrinoFourMomentumNeg.E() - trueNeutrinoFourMomentum.E() ) ? recoNeutrinoFourMomentumPos : recoNeutrinoFourMomentumNeg );
	streamlog_out(DEBUG4) << "" << std::endl;
	streamlog_out(DEBUG4) << "	Closest Neutrino 4-Momentum:			( " << recoNeutrinoFourMomentumClose.Px() << "	, " << recoNeutrinoFourMomentumClose.Py() << "	, " << recoNeutrinoFourMomentumClose.Pz() << "	, " << recoNeutrinoFourMomentumClose.E() << " )" << std::endl;
	streamlog_out(DEBUG4) << "	True Neutrino 4-Momentum:			( " << trueNeutrinoFourMomentum.Px() << "	, " << trueNeutrinoFourMomentum.Py() << "	, " << trueNeutrinoFourMomentum.Pz() << "	, " << trueNeutrinoFourMomentum.E() << " )" << std::endl;
	if ( fabs( recoNeutrinoFourMomentumClose.E() - trueNeutrinoFourMomentum.E() ) > 10.0 )
	{
		streamlog_out(DEBUG4) << "	!!! Big Difference between true and reco neutrino Energy : " << recoNeutrinoFourMomentumClose.E() - trueNeutrinoFourMomentum.E() << "  GeV" << std::endl;
	}
	if ( m_fillRootTree ) plotHistograms( trueNeutrinoFourMomentum , recoNeutrinoFourMomentumClose , NeutrinoCovMat );
	m_trueNuPx.push_back( trueNeutrinoFourMomentum.Px() );
	m_trueNuPy.push_back( trueNeutrinoFourMomentum.Py() );
	m_trueNuPz.push_back( trueNeutrinoFourMomentum.Pz() );
	m_trueNuE.push_back( trueNeutrinoFourMomentum.E() );
	m_recoNuClosePx.push_back( recoNeutrinoFourMomentumClose.Px() );
	m_recoNuClosePy.push_back( recoNeutrinoFourMomentumClose.Py() );
	m_recoNuClosePz.push_back( recoNeutrinoFourMomentumClose.Pz() );
	m_recoNuCloseE.push_back( recoNeutrinoFourMomentumClose.E() );
	m_recoNuPosPx.push_back( recoNeutrinoFourMomentumPos.Px() );
	m_recoNuPosPy.push_back( recoNeutrinoFourMomentumPos.Py() );
	m_recoNuPosPz.push_back( recoNeutrinoFourMomentumPos.Pz() );
	m_recoNuPosE.push_back( recoNeutrinoFourMomentumPos.E() );
	m_recoNuNegPx.push_back( recoNeutrinoFourMomentumNeg.Px() );
	m_recoNuNegPy.push_back( recoNeutrinoFourMomentumNeg.Py() );
	m_recoNuNegPz.push_back( recoNeutrinoFourMomentumNeg.Pz() );
	m_recoNuNegE.push_back( recoNeutrinoFourMomentumNeg.E() );


/*
	streamlog_out(DEBUG2) << "			     (  PDG	, Mass		, Px		, Py		, Pz		, E		, Charge	)" << std::endl;
	streamlog_out(DEBUG2) << "		Neutrino" << std::endl;
	streamlog_out(DEBUG2) << "			True:(	" << "***" << "	, " << trueNeutrinoFourMomentum.M() << "	, " << trueNeutrinoFourMomentum.Px() << "	, " << trueNeutrinoFourMomentum.Py() << "	, " << trueNeutrinoFourMomentum.Pz() << "	, " << trueNeutrinoFourMomentum.E() << "	, " << "0" << "	)" << std::endl;
	streamlog_out(DEBUG2) << "" << std::endl;
	streamlog_out(DEBUG2) << "		Hadron" << std::endl;
	streamlog_out(DEBUG2) << "			True:(	" << parentHadron->getPDG() << "	, " << parentHadron->getMass() << "	, " << parentHadron->getMomentum()[ 0 ] << "	, " << parentHadron->getMomentum()[ 1 ] << "	, " << parentHadron->getMomentum()[ 2 ] << "	, " << parentHadron->getEnergy() << "	, " << parentHadron->getCharge() << "	)" << std::endl;
	streamlog_out(DEBUG2) << "" << std::endl;
	getLeptonFourMomentum( pLCEvent , SLDLepton  , m_cheatLepton4momentum , leptonFourMomentum , trueLeptonFourMomentum , m_RecoMCTruthLinkCollection , m_MCTruthRecoLinkCollection );
	streamlog_out(DEBUG2) << "			Used:(	" << SLDLepton->getPDG() << "	, " << leptonFourMomentum.M() << "	, " << leptonFourMomentum.Px() << "	, " << leptonFourMomentum.Py() << "	, " << leptonFourMomentum.Pz() << "	, " << leptonFourMomentum.E() << "	, " << (int)SLDLepton->getCharge() << "		)" << std::endl;
	streamlog_out(DEBUG2) << "" << std::endl;
	getChargedFourMomentum( pLCEvent , SLDLepton  , m_cheatCharged4momentum , chargedFourMomentum , trueChargedFourMomentum , restCharge , m_RecoMCTruthLinkCollection , m_MCTruthRecoLinkCollection );
	streamlog_out(DEBUG2) << "			Used:(	" << "qqq" << "	, " << chargedFourMomentum.M() << "	, " << chargedFourMomentum.Px() << "	, " << chargedFourMomentum.Py() << "	, " << chargedFourMomentum.Pz() << "	, " << chargedFourMomentum.E() << "	, " << (int)restCharge << "		)" << std::endl;
	streamlog_out(DEBUG2) << "" << std::endl;
	getNeutralFourMomentum( pLCEvent , SLDLepton  , m_cheatNeutral4momentum , neutralFourMomentum , trueNeutralFourMomentum , m_RecoMCTruthLinkCollection , m_MCTruthRecoLinkCollection );
	streamlog_out(DEBUG2) << "			Used:(	" << "nnn" << "	, " << neutralFourMomentum.M() << "	, " << neutralFourMomentum.Px() << "	, " << neutralFourMomentum.Py() << "	, " << neutralFourMomentum.Pz() << "	, " << neutralFourMomentum.E() << "	, " << "0" << "		)" << std::endl;
	streamlog_out(DEBUG2) << "" << std::endl;
	streamlog_out(DEBUG2) << "			     (  X		, Y		, Z	)" << std::endl;


	ReconstructedParticle *assignedJet = NULL;

	int flightDirectionStatus = 0;
	if ( linkedRecoLepton != NULL ) flightDirectionStatus = getRecoFlightDirection( linkedRecoLepton , recoFlightDirection , primaryVertex , startVertex , SLDVertex , SLDVertexRP , assignedJet , jetVector , buildUpVertexVector );

	if ( assignedJet != NULL ) streamlog_out(DEBUG1) << *assignedJet << std::endl;
	if ( SLDVertex != NULL ) streamlog_out(DEBUG1) << *SLDVertex << std::endl;

	std::vector<EVENT::MCParticle*> trueNeutralDecayProducts{};
	getNeutralMCPs( SLDLepton->getParents()[ 0 ] , trueNeutralDecayProducts );


	std::vector<EVENT::MCParticle*> trueChargedDecayProducts{};
	getChargedMCPs( SLDLepton->getParents()[ 0 ] , trueChargedDecayProducts );

	std::vector<EVENT::ReconstructedParticle*> aloneChargedPFOs{};
	std::vector<EVENT::ReconstructedParticle*> recoChargedDecayProducts{};
	if ( assignedJet != NULL ) aloneChargedPFOs = getParticlesWithAloneTracks( linkedRecoLepton , assignedJet , primaryVertex , buildUpVertexVector );
	if ( aloneChargedPFOs.size() != 0 )
	{
		for ( unsigned int i_trk = 0 ; i_trk < aloneChargedPFOs.size() ; ++i_trk )
		{
			ReconstructedParticle* testPFO = (ReconstructedParticle*) aloneChargedPFOs[ i_trk ];
			MCParticle* linkedMCP = getLinkedMCP( pLCEvent , testPFO , m_RecoMCTruthLinkCollection , m_MCTruthRecoLinkCollection , true , false );
			for ( unsigned int i_mcp = 0 ; i_mcp < trueChargedDecayProducts.size() ; ++i_mcp )
			{
				MCParticle* testMCP = (MCParticle*) trueChargedDecayProducts[ i_mcp ];
				if ( testMCP == linkedMCP ) recoChargedDecayProducts.push_back( testPFO );
			}
		}
	}




//	int flightDirectionStatus = 0;
//	flightDirectionStatus = getParentHadronFlightDirection( pLCEvent , SLDLepton , trueFlightDirection , recoFlightDirection , m_inputPrimaryVertex , m_inputBuildUpVertex , m_inputJetCollection , m_vertexingScenario , m_RecoMCTruthLinkCollection , m_MCTruthRecoLinkCollection , helicesDistance , SecondaryVertexPar , m_displayEvent );
//	if ( helicesDistance > 400.0 ) getParentHadronFlightDirection( pLCEvent , SLDLepton , trueFlightDirection , recoFlightDirection , m_inputPrimaryVertex , m_inputBuildUpVertex , m_inputJetCollection , m_vertexingScenario , m_RecoMCTruthLinkCollection , m_MCTruthRecoLinkCollection , helicesDistance , SecondaryVertexPar , true , this );
	if ( m_fillRootTree )
	{
		h_secondaryVertex->Fill( flightDirectionStatus - 0.5 );
		++n_secondaryVertex;
	}
*/

/*
	m_flightDirectionStatus.push_back( flightDirectionStatus );
	m_distRecoLeptonToDownStreamVertex.push_back( helicesDistance );
	m_FlightDirectionErrorCosAlpha.push_back( trueFlightDirection.Dot( recoFlightDirection ) );
	m_FlightDirectionErrorSinAlpha.push_back( sqrt( 1 - pow( trueFlightDirection.Dot( recoFlightDirection ) , 2 ) ) );
	m_FlightDirectionErrorAlpha.push_back( acos( trueFlightDirection.Dot( recoFlightDirection ) ) * 180.0 / 3.14159265 );
	m_dsVertexResidualX.push_back( SecondaryVertexPar[ 3 ] - SecondaryVertexPar[ 0 ] );
	m_dsVertexResidualY.push_back( SecondaryVertexPar[ 4 ] - SecondaryVertexPar[ 1 ] );
	m_dsVertexResidualZ.push_back( SecondaryVertexPar[ 5 ] - SecondaryVertexPar[ 2 ] );
	m_SecVertexResidualX.push_back( SecondaryVertexPar[ 9 ] - SecondaryVertexPar[ 6 ] );
	m_SecVertexResidualY.push_back( SecondaryVertexPar[ 10 ] - SecondaryVertexPar[ 7 ] );
	m_SecVertexResidualZ.push_back( SecondaryVertexPar[ 11 ] - SecondaryVertexPar[ 8 ] );
*/

/*
	if ( m_cheatFlightDirection )
	{
		flightDirection = trueFlightDirection;
	}
	else
	{
		flightDirection = recoFlightDirection;
	}
	streamlog_out(DEBUG4) << "" << std::endl;
	streamlog_out(DEBUG4) << "		Flight Direction Error:		CosAlpha = " << trueFlightDirection.Dot( recoFlightDirection ) << "	, Alpha = " << acos( trueFlightDirection.Dot( recoFlightDirection ) ) * 180.0 / 3.14159265 << " deg" << std::endl;
//	if ( flightDirectionStatus == 1 ) return;
	recoNeutrinoFourMomentumPos = getNeutrinoFourMomentum( flightDirection , leptonFourMomentum , chargedFourMomentum , neutralFourMomentum , parentHadronMass , +1 );
	recoNeutrinoFourMomentumNeg = getNeutrinoFourMomentum( flightDirection , leptonFourMomentum , chargedFourMomentum , neutralFourMomentum , parentHadronMass , -1 );
*/

	if ( m_displayEvent )
	{
		DDMarlinCED::draw( this , 1); // draw everything
	}
	recoHadronFourMomentum = recoVisibleFourMomentum + recoNeutrinoFourMomentumClose;
	fillTrueRecoFourMomentum( trueNeutralFourMomentum , trueChargedFourMomentum , trueLeptonFourMomentum , trueVisibleFourMomentum , trueNeutrinoFourMomentum , trueHadronFourMomentum , recoNeutralFourMomentum , recoChargedFourMomentum , recoLeptonFourMomentum , recoVisibleFourMomentum , recoNeutrinoFourMomentumClose , recoHadronFourMomentum );



}

void SLDCorrection::showTrueParameters( MCParticle *SLDLepton )
{
	TLorentzVector true4mom( 0.0 , 0.0 , 0.0 , 0.0 );
	MCParticle* parentHadron = SLDLepton->getParents()[ 0 ];
	streamlog_out(DEBUG4) << "	PARENT HADRON:" << std::endl;
	streamlog_out(DEBUG4) << *parentHadron << std::endl;
	TVector3 trueFliDir = TVector3( parentHadron->getMomentumAtEndpoint() );
	trueFliDir.SetMag( 1.0 );
	for ( unsigned int i_mcp = 0 ; i_mcp < parentHadron->getDaughters().size() ; ++i_mcp )
	{
		MCParticle* daughter = parentHadron->getDaughters()[ i_mcp ];
		streamlog_out(DEBUG4) << *daughter << std::endl;
		if ( std::fabs( daughter->getPDG() ) != 12 && std::fabs( daughter->getPDG() ) != 14 && std::fabs( daughter->getPDG() ) != 16 )
		{
			true4mom += TLorentzVector( daughter->getMomentum() , daughter->getEnergy() );
			streamlog_out(DEBUG4) << " ONE MCPARTICLE IS ADDED" << std::endl;
		}
	}
	streamlog_out(DEBUG4) << "	TRUE VISIBLE FOUR-MOMENTUM (Px,Py,Pz,E):	(	" << true4mom.Px() << "	,	" << true4mom.Py() << "	,	" << true4mom.Pz() << "	,	" << true4mom.E() << "	)" << std::endl;
	streamlog_out(DEBUG4) << "	TRUE FLIHT DIRECTION (x,y,z): 		(	" << trueFliDir.X() << "	,	" << trueFliDir.Y() << "	,	" << trueFliDir.Z() << "	)" << std::endl;
	streamlog_out(DEBUG4) << "	TRUE PARENT HADRON MASS =  			" << parentHadron->getMass() << std::endl;

}

TLorentzVector SLDCorrection::getNeutrinoFourMomentum( TVector3 flightDirection , TLorentzVector visibleFourMomentum , double ParentHadronMass , float solutionSign )
{
	int sign = ( solutionSign != 0 ? solutionSign / abs( solutionSign ) : 1 );
	m_solutionSign.push_back( sign );
	const char *solSign = ( sign >= 0 ? "+" : "-" );
	streamlog_out(DEBUG4) << "" << std::endl;
	streamlog_out(DEBUG4) << "		--------------------------------------------" << std::endl;
	streamlog_out(DEBUG4) << "		Calculate Neutrino 4-Momentum for " << solSign << " solution" << std::endl;
	streamlog_out(DEBUG4) << "		--------------------------------------------" << std::endl;

	streamlog_out(DEBUG4) << "		Test 1, |flightDirection| = " << flightDirection.Mag() << std::endl;
	flightDirection.SetMag( 1.0 );
	streamlog_out(DEBUG4) << "		flightDirection:			( " << flightDirection.Px() << "	, " << flightDirection.Py() << "	, " << flightDirection.Pz() << "	)" << std::endl;
	streamlog_out(DEBUG4) << "		Parent Hadron Mass =	 " << ParentHadronMass << std::endl;

	streamlog_out(DEBUG4) << "		Visible 4-Momentum:			( " << visibleFourMomentum.Px() << "	, " << visibleFourMomentum.Py() << "	, " << visibleFourMomentum.Pz() << "	, " << visibleFourMomentum.E() << " )" << std::endl;

	double visible_mass		= visibleFourMomentum.M();
	streamlog_out(DEBUG4) << "		Visible Inv Mass:	" << visible_mass << std::endl;

	double visible_E		= visibleFourMomentum.E();
	m_E_vis.push_back( visible_E );
	streamlog_out(DEBUG4) << "		Visible Energy:									" << visible_E << std::endl;

	TVector3 visible_p		= TVector3( visibleFourMomentum.Px() , visibleFourMomentum.Py() , visibleFourMomentum.Pz() );
	streamlog_out(DEBUG4) << "		Visible Momentum:			( " << visible_p.Px() << "	, " << visible_p.Py() << "	, " << visible_p.Pz() << "	)" << std::endl;
	streamlog_out(DEBUG4) << "		(Visible Momentum).(FlightDirection):							" << visible_p.Dot( flightDirection ) << std::endl;

	TVector3 visible_p_par		= visible_p.Dot( flightDirection ) * flightDirection;
	m_P_vis_par.push_back( visible_p_par.Mag() );
	streamlog_out(DEBUG4) << "		Visible Momentum (par):			( " << visible_p_par.Px() << "	, " << visible_p_par.Py() << "	, " << visible_p_par.Pz() << "	)" << std::endl;

	TVector3 visible_p_nor		= visible_p - visible_p_par;
	m_P_vis_nor.push_back( visible_p_nor.Mag() );
	m_P_vis_nor_prime.push_back( visible_p_nor.Mag() );
	streamlog_out(DEBUG4) << "		Visible Momentum (nor):			( " << visible_p_nor.Px() << "	, " << visible_p_nor.Py() << "	, " << visible_p_nor.Pz() << "	)" << std::endl;

	double visible_E_prime		= ( pow( ParentHadronMass , 2 ) + pow( visible_mass , 2 ) ) / ( 2 * ParentHadronMass );
	m_E_vis_prime.push_back( visible_E_prime );
	streamlog_out(DEBUG4) << "		Visible Energy (prime):								" << visible_E_prime << std::endl;

	TVector3 visible_p_par_prime	= solutionSign * sqrt( pow( ( pow( ParentHadronMass , 2 ) - pow( visible_mass , 2 ) ) / ( 2 * ParentHadronMass ) , 2 ) - visible_p_nor.Mag2() ) * flightDirection;
	m_P_vis_par_prime.push_back( visible_p_par_prime.Mag() );
	streamlog_out(DEBUG4) << "		Visible Momentum (par-prime):		( " << visible_p_par_prime.Px() << "	, " << visible_p_par_prime.Py() << "	, " << visible_p_par_prime.Pz() << "	)" << std::endl;
	if ( pow( ( pow( ParentHadronMass , 2 ) - pow( visible_mass , 2 ) ) / ( 2 * ParentHadronMass ) , 2 ) < visible_p_nor.Mag2() )
	{
		visible_p_par_prime	= solutionSign * std::numeric_limits<double>::min() * flightDirection;
	}
	streamlog_out(DEBUG4) << "		Visible Momentum (par-prime):		( " << visible_p_par_prime.Px() << "	, " << visible_p_par_prime.Py() << "	, " << visible_p_par_prime.Pz() << "	)" << std::endl;

	double parent_hadron_E		= ( ( visibleFourMomentum.E() * ( pow( ParentHadronMass , 2 ) + pow( visible_mass , 2 ) ) / ( 2 * ParentHadronMass ) ) - visible_p_par.Dot( visible_p_par_prime ) ) * ParentHadronMass / ( pow( visible_mass , 2 ) + visible_p_nor.Mag2() );
	streamlog_out(DEBUG4) << "		Parent Hadron Energy =									" << parent_hadron_E << std::endl;
	TVector3 parent_hadron_p	= sqrt( pow( parent_hadron_E , 2 ) - pow( ParentHadronMass , 2 ) ) * flightDirection;
	streamlog_out(DEBUG4) << "		Parent Hadron Momentum:			( " << parent_hadron_p.Px() << "	, " << parent_hadron_p.Py() << "	, " << parent_hadron_p.Pz() << "	, " << parent_hadron_E << " )" << std::endl;

	double Neutrino_E		= parent_hadron_E - visible_E;
	streamlog_out(DEBUG4) << "		Neutrino Energy:									" << Neutrino_E << std::endl;

	TVector3 Neutrino_p_nor		= -1 * visible_p_nor;
	streamlog_out(DEBUG4) << "		Neutrino Momentum (nor):		( " << Neutrino_p_nor.Px() << "	, " << Neutrino_p_nor.Py() << "	, " << Neutrino_p_nor.Pz() << "	)" << std::endl;

	TVector3 Neutrino_p_par		= sqrt( pow( Neutrino_E , 2 ) - Neutrino_p_nor.Mag2() ) * flightDirection;
	streamlog_out(DEBUG4) << "		Neutrino Momentum (par):		( " << Neutrino_p_par.Px() << "	, " << Neutrino_p_par.Py() << "	, " << Neutrino_p_par.Pz() << "	)" << std::endl;

	TVector3 Neutrino_p		= Neutrino_p_nor + Neutrino_p_par;
	streamlog_out(DEBUG4) << "		Neutrino Momentum:			( " << Neutrino_p_par.Px() << "	, " << Neutrino_p_par.Py() << "	, " << Neutrino_p_par.Pz() << "	)" << std::endl;

	TLorentzVector Neutrino_tlv( Neutrino_p , Neutrino_E );
	streamlog_out(DEBUG2) << "" << std::endl;
	streamlog_out(DEBUG2) << "	Reconstructed Neutrino Mass:	" << Neutrino_tlv.Mag() << std::endl;
	streamlog_out(DEBUG2) << "	Reconstructed Neutrino 4-Momentum:		( " << Neutrino_tlv.Px() << "	, " << Neutrino_tlv.Py() << "	, " << Neutrino_tlv.Pz() << "	, " << Neutrino_tlv.E() << " )" << std::endl;
	return Neutrino_tlv;
}

MCParticle* SLDCorrection::getTrueNeutrino( MCParticle *SLDLepton , TLorentzVector& InisibleFourMomentum )
{
	MCParticle* trueNeutrino{};
//	TLorentzVector InisibleFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	try
	{
		EVENT::MCParticle *MotherHadron = SLDLepton->getParents()[ 0 ];
		int nNeutrinos = 0;
		for ( long unsigned int i_daughter = 0 ; i_daughter < ( MotherHadron->getDaughters() ).size() ; ++i_daughter )
		{
			EVENT::MCParticle *daughter = MotherHadron->getDaughters()[ i_daughter ];
			if ( daughter->getGeneratorStatus() == 1 && ( abs( daughter->getPDG() ) == abs( SLDLepton->getPDG() ) + 1 ) )
			{
				streamlog_out(DEBUG0) << "----------------------------------------------------------------------" << std::endl;
				streamlog_out(DEBUG0) << "------------------------------ Neutrino ------------------------------" << std::endl;
				streamlog_out(DEBUG0) << "----------------------------------------------------------------------" << std::endl;
				streamlog_out(DEBUG0) << *daughter << std::endl;
				trueNeutrino = daughter;
				InisibleFourMomentum = TLorentzVector( daughter->getMomentum()[ 0 ] , daughter->getMomentum()[ 1 ] , daughter->getMomentum()[ 2 ] , daughter->getEnergy() );
			}
		}
		++nNeutrinos;
	}
	catch(DataNotAvailableException &e)
        {
        	streamlog_out(MESSAGE) << "	True Neutrino for semi-leptonic decay not found in MCParticles" << std::endl;
        }
	return trueNeutrino;
}

void SLDCorrection::fillTrueRecoFourMomentum(	TLorentzVector trueNeutralFourMomentum , TLorentzVector trueChargedFourMomentum ,
						TLorentzVector trueLeptonFourMomentum , TLorentzVector trueVisibleFourMomentum ,
						TLorentzVector trueNeutrinoFourMomentum , TLorentzVector trueHadronFourMomentum ,
						TLorentzVector recoNeutralFourMomentum , TLorentzVector recoChargedFourMomentum ,
						TLorentzVector recoLeptonFourMomentum , TLorentzVector recoVisibleFourMomentum ,
						TLorentzVector recoNeutrinoFourMomentum , TLorentzVector recoHadronFourMomentum )
{
	m_trueNeutralPx.push_back( trueNeutralFourMomentum.Px() );
	m_trueNeutralPy.push_back( trueNeutralFourMomentum.Py() );
	m_trueNeutralPz.push_back( trueNeutralFourMomentum.Pz() );
	m_trueNeutralE.push_back( trueNeutralFourMomentum.E() );
	m_trueChargedPx.push_back( trueChargedFourMomentum.Px() );
	m_trueChargedPy.push_back( trueChargedFourMomentum.Py() );
	m_trueChargedPz.push_back( trueChargedFourMomentum.Pz() );
	m_trueChargedE.push_back( trueChargedFourMomentum.E() );
	m_trueLeptonPx.push_back( trueLeptonFourMomentum.Px() );
	m_trueLeptonPy.push_back( trueLeptonFourMomentum.Py() );
	m_trueLeptonPz.push_back( trueLeptonFourMomentum.Pz() );
	m_trueLeptonE.push_back( trueLeptonFourMomentum.E() );
	m_trueVisiblePx.push_back( trueVisibleFourMomentum.Px() );
	m_trueVisiblePy.push_back( trueVisibleFourMomentum.Py() );
	m_trueVisiblePz.push_back( trueVisibleFourMomentum.Pz() );
	m_trueVisibleE.push_back( trueVisibleFourMomentum.E() );
	m_trueNeutrinoPx.push_back( trueNeutrinoFourMomentum.Px() );
	m_trueNeutrinoPy.push_back( trueNeutrinoFourMomentum.Py() );
	m_trueNeutrinoPz.push_back( trueNeutrinoFourMomentum.Pz() );
	m_trueNeutrinoE.push_back( trueNeutrinoFourMomentum.E() );
	m_trueHadronPx.push_back( trueHadronFourMomentum.Px() );
	m_trueHadronPy.push_back( trueHadronFourMomentum.Py() );
	m_trueHadronPz.push_back( trueHadronFourMomentum.Pz() );
	m_trueHadronE.push_back( trueHadronFourMomentum.E() );
	m_recoNeutralPx.push_back( recoNeutralFourMomentum.Px() );
	m_recoNeutralPy.push_back( recoNeutralFourMomentum.Py() );
	m_recoNeutralPz.push_back( recoNeutralFourMomentum.Pz() );
	m_recoNeutralE.push_back( recoNeutralFourMomentum.E() );
	m_recoChargedPx.push_back( recoChargedFourMomentum.Px() );
	m_recoChargedPy.push_back( recoChargedFourMomentum.Py() );
	m_recoChargedPz.push_back( recoChargedFourMomentum.Pz() );
	m_recoChargedE.push_back( recoChargedFourMomentum.E() );
	m_recoLeptonPx.push_back( recoLeptonFourMomentum.Px() );
	m_recoLeptonPy.push_back( recoLeptonFourMomentum.Py() );
	m_recoLeptonPz.push_back( recoLeptonFourMomentum.Pz() );
	m_recoLeptonE.push_back( recoLeptonFourMomentum.E() );
	m_recoVisiblePx.push_back( recoVisibleFourMomentum.Px() );
	m_recoVisiblePy.push_back( recoVisibleFourMomentum.Py() );
	m_recoVisiblePz.push_back( recoVisibleFourMomentum.Pz() );
	m_recoVisibleE.push_back( recoVisibleFourMomentum.E() );
	m_recoNeutrinoPx.push_back( recoNeutrinoFourMomentum.Px() );
	m_recoNeutrinoPy.push_back( recoNeutrinoFourMomentum.Py() );
	m_recoNeutrinoPz.push_back( recoNeutrinoFourMomentum.Pz() );
	m_recoNeutrinoE.push_back( recoNeutrinoFourMomentum.E() );
	m_recoHadronPx.push_back( recoHadronFourMomentum.Px() );
	m_recoHadronPy.push_back( recoHadronFourMomentum.Py() );
	m_recoHadronPz.push_back( recoHadronFourMomentum.Pz() );
	m_recoHadronE.push_back( recoHadronFourMomentum.E() );
}

void SLDCorrection::plotHistograms( TLorentzVector trueFourMomentumNeutrino , TLorentzVector FourMomentumNuClose , std::vector<float> NeutrinoCovMat )
{
	double NuPxResidual = FourMomentumNuClose.Px() - trueFourMomentumNeutrino.Px(); m_NuPxResidual.push_back( NuPxResidual );
	double NuPyResidual = FourMomentumNuClose.Py() - trueFourMomentumNeutrino.Py(); m_NuPyResidual.push_back( NuPyResidual );
	double NuPzResidual = FourMomentumNuClose.Pz() - trueFourMomentumNeutrino.Pz(); m_NuPzResidual.push_back( NuPzResidual );
	double NuEResidual = FourMomentumNuClose.E() - trueFourMomentumNeutrino.E(); m_NuEResidual.push_back( NuEResidual );
//	double NuPxNormalizedResidual = NuPxResidual / sqrt( NeutrinoCovMat[ 0 ] ); m_NuPxNormalizedResidual.push_back( NuPxNormalizedResidual );
//	double NuPyNormalizedResidual = NuPyResidual / sqrt( NeutrinoCovMat[ 2 ] ); m_NuPyNormalizedResidual.push_back( NuPyNormalizedResidual );
//	double NuPzNormalizedResidual = NuPzResidual / sqrt( NeutrinoCovMat[ 5 ] ); m_NuPzNormalizedResidual.push_back( NuPzNormalizedResidual );
//	double NuENormalizedResidual = NuEResidual / sqrt( NeutrinoCovMat[ 9 ] ); m_NuENormalizedResidual.push_back( NuENormalizedResidual );
	h_NuPxResidual->Fill( NuPxResidual ); ++n_NuPxResidual;
	h_NuPyResidual->Fill( NuPyResidual ); ++n_NuPyResidual;
	h_NuPzResidual->Fill( NuPzResidual ); ++n_NuPzResidual;
	h_NuEResidual->Fill( NuEResidual ); ++n_NuEResidual;
//	h_NuPxNormalizedResidual->Fill( NuPxNormalizedResidual ); ++n_NuPxNormalizedResidual;
//	h_NuPyNormalizedResidual->Fill( NuPyNormalizedResidual ); ++n_NuPyNormalizedResidual;
//	h_NuPzNormalizedResidual->Fill( NuPzNormalizedResidual ); ++n_NuPzNormalizedResidual;
//	h_NuENormalizedResidual->Fill( NuENormalizedResidual ); ++n_NuENormalizedResidual;
	h_recoNuPx_mcNuPx->Fill( trueFourMomentumNeutrino.Px() , FourMomentumNuClose.Px() );
	h_recoNuPy_mcNuPy->Fill( trueFourMomentumNeutrino.Py() , FourMomentumNuClose.Py() );
	h_recoNuPz_mcNuPz->Fill( trueFourMomentumNeutrino.Pz() , FourMomentumNuClose.Pz() );
	h_recoNuE_mcNuE->Fill( trueFourMomentumNeutrino.E() , FourMomentumNuClose.E() );
}

void SLDCorrection::InitializeHistogram( TH1F *histogram , int scale , int color , int lineWidth , int markerSize , int markerStyle )
{
	histogram->Scale( 1.0 / scale );
	histogram->SetLineColor( color );
	histogram->SetLineWidth( lineWidth );
	histogram->SetMarkerSize( markerSize );
	histogram->SetMarkerStyle( markerStyle );
	histogram->SetMarkerColor( color );
	float fit_range = 4.0;
	float fit_min = -2.0;
	float fit_max = 2.0;
	doProperGaussianFit( histogram , fit_min , fit_max , fit_range );
	histogram->GetFunction("gaus")->SetLineColor( color );
	float y_max = 1.2 * histogram->GetMaximum();
	histogram->GetYaxis()->SetRangeUser(0.0, y_max);
	histogram->Write();
}

void SLDCorrection::doProperGaussianFit( TH1F *histogram , float fitMin , float fitMax , float fitRange )
{
	float Chi2 = 0.0;
	float NDF = 0.0;
	for ( int i_fit = 0 ; i_fit < 3 ; ++i_fit )
	{
		histogram->Fit( "gaus" , "" , "" , fitMin , fitMax );
		TF1 *fitFunction = (TF1 *)histogram->GetFunction("gaus");
		double fitMean = fitFunction->GetParameter( 1 );
		double fitSigma = fitFunction->GetParameter( 2 );
		fitMin = fitMean - fitRange * fitSigma;
		fitMax = fitMean + fitRange * fitSigma;
		Chi2 = fitFunction->GetChisquare();
		NDF = fitFunction->GetNDF();
	}
	streamlog_out(DEBUG2) << "	FIT : CHI2(" << Chi2 << ") / NDF(" << NDF << ") = " << Chi2 / NDF << " 	, fitrange = " << fitRange << std::endl;
	streamlog_out(DEBUG2) << "" << std::endl;
	if ( Chi2 != 0.0 && NDF != 0.0 && Chi2 / NDF > 2.0 && fitRange >= 0.5 )
	{
		doProperGaussianFit( histogram , fitMin , fitMax , fitRange - 0.1 );
	}
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
		m_pTTree1->Write();
		m_pTTree2->Write();
		h_SLDStatus->Scale( 100.0 / n_SLDStatus );
		h_SLDStatus->GetYaxis()->SetTitle("#SLDecay [%]");
		h_SLDStatus->Write();
		InitializeHistogram( h_NuPxResidual , n_NuPxResidual , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_NuPyResidual , n_NuPyResidual , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_NuPzResidual , n_NuPzResidual , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_NuEResidual , n_NuEResidual , 4 , 1 , 1.0 , 1 );
		h_SLDecayOrder->Write();
		h_NuPxResidual->Write();
		h_NuPyResidual->Write();
		h_NuPzResidual->Write();
		h_NuEResidual->Write();
		h_recoNuPx_mcNuPx->Write();
		h_recoNuPy_mcNuPy->Write();
		h_recoNuPz_mcNuPz->Write();
		h_recoNuE_mcNuE->Write();
		h_secondaryVertex->Scale( 100.0 / n_secondaryVertex );
		h_secondaryVertex->GetYaxis()->SetTitle("#SLDecay [%]");
		h_secondaryVertex->Write();
		m_pTFile->Close();
		delete m_pTFile;
	}
}
