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
	m_nSLDecayOfBHadron(0),
	m_nSLDecayOfCHadron(0),
	m_nSLDecayOfTauLepton(0),
	m_nSLDecayTotal(0),
	m_nSLDecayToElectron(0),
	m_nSLDecayToMuon(0),
	m_nSLDecayToTau(0),
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

	registerOutputCollection(	LCIO::VERTEX,
					"SemiLeptonicDecayVertex",
					"Name of Semi-Leptonic Decay Vertices Collection",
					m_SLDVertex,
					std::string("SemiLeptonicDecayVertex")
				);

	registerOutputCollection(	LCIO::RECONSTRUCTEDPARTICLE,
					"SemiLeptonicDecayVertexRP",
					"Name of Semi-Leptonic Decay Vertices Reconstructed Particle Collection",
					m_SLDVertexRP,
					std::string("SemiLeptonicDecayVertexRP")
				);

	registerOutputCollection(	LCIO::RECONSTRUCTEDPARTICLE,
					"ReconstructedNeutrino",
					"Name of Reconstructed Neutrino Collection",
					m_reconstructedNeutrino,
					std::string("ReconstructedNeutrino")
				);

	registerOutputCollection(	LCIO::LCRELATION,
					"JetSLDLinkName",
					"Name of the JetSemiLeptonicDecayLinkName output collection",
					m_JetSLDLinkName,
					std::string("JetSLDLinkName")
				);

	registerOutputCollection(	LCIO::LCRELATION,
					"SLDJetLinkName",
					"Name of the SemiLeptonicDecayJetLinkName output collection",
					m_SLDJetLinkName,
					std::string("SLDJetLinkName")
				);

	registerOutputCollection(	LCIO::LCRELATION,
					"mcNurecoNuLinkName",
					"Name of the trueNeutrino-reconstructedNeutrino output Link collection",
					m_mcNurecoNuLinkName,
					std::string("mcNurecoNuLinkName")
				);

	registerOutputCollection(	LCIO::LCRELATION,
					"recoNumcNuLinkName",
					"Name of the trueNeutrino-reconstructedNeutrino output Link collection",
					m_recoNumcNuLinkName,
					std::string("recoNumcNuLinkName")
				);

	registerOutputCollection(	LCIO::LCRELATION,
					"SLDNuLinkName",
					"Name of the NeutrinoSemiLeptonicDecayLinkName output collection",
					m_SLDNuLinkName,
					std::string("SLDNuLinkName")
				);

	registerOutputCollection(	LCIO::LCRELATION,
					"NuSLDLinkName",
					"Name of the SemiLeptonicDecayNeutrinoLinkName output collection",
					m_NuSLDLinkName,
					std::string("NuSLDLinkName")
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

	registerProcessorParameter(	"cheatPVAcharged",
					"Cheat Particle ID for Charged Decay Products",
					m_cheatPVAcharged,
					bool(true)
				);

	registerProcessorParameter(	"chargedCosAcceptanceAngleSLD4",
					"Acceptance angle for charged PFOs to be assigned to the semi-leptonic decay when lepton is in secondary vertex",
					m_chargedCosAcceptanceAngleSLD4,
					float(0.0)
				);

	registerProcessorParameter(	"chargedCosAcceptanceAngleSLD5",
					"Acceptance angle for charged PFOs to be assigned to the semi-leptonic decay when lepton is with third vertex",
					m_chargedCosAcceptanceAngleSLD5,
					float(0.0)
				);

	registerProcessorParameter(	"cheatPVAneutral",
					"Cheat Particle ID for Neutral Decay Products",
					m_cheatPVAneutral,
					bool(true)
				);

	registerProcessorParameter(	"neutralCosAcceptanceAngle",
					"Acceptance angle for neutral PFOs to be assigned to the semi-leptonic decay",
					m_neutralCosAcceptanceAngle,
					float(0.0)
				);

	registerProcessorParameter(	"BSLDChargedSLD4InvMassCut",
					"Cut on Invariant mass of visible charged decay products for semi-leptonic decay of B-Hadron when lepton is in secondary vertex",
					m_BSLDChargedSLD4InvMassCut,
					float(0.0)
				);

	registerProcessorParameter(	"BSLDChargedSLD5InvMassCut",
					"Cut on Invariant mass of visible charged decay products for semi-leptonic decay of B-Hadron when lepton is with third vertex",
					m_BSLDChargedSLD5InvMassCut,
					float(0.0)
				);

	registerProcessorParameter(	"BSLDNeutralSLD4InvMassCut",
					"Cut on Invariant mass of visible neutral decay products for semi-leptonic decay of B-Hadron when lepton is in secondary vertex",
					m_BSLDNeutralSLD4InvMassCut,
					float(0.0)
				);

	registerProcessorParameter(	"BSLDNeutralSLD5InvMassCut",
					"Cut on Invariant mass of visible neutral decay products for semi-leptonic decay of B-Hadron when lepton is with third vertex",
					m_BSLDNeutralSLD5InvMassCut,
					float(0.0)
				);

	registerProcessorParameter(	"CSLDChargedSLD4InvMassCut",
					"Cut on Invariant mass of visible charged decay products for semi-leptonic decay of C-Hadron when lepton is in secondary vertex",
					m_CSLDChargedSLD4InvMassCut,
					float(0.0)
				);

	registerProcessorParameter(	"CSLDChargedSLD5InvMassCut",
					"Cut on Invariant mass of visible charged decay products for semi-leptonic decay of C-Hadron when lepton is with third vertex",
					m_CSLDChargedSLD5InvMassCut,
					float(0.0)
				);

	registerProcessorParameter(	"CSLDNeutralSLD4InvMassCut",
					"Cut on Invariant mass of visible neutral decay products for semi-leptonic decay of C-Hadron when lepton is in secondary vertex",
					m_CSLDNeutralSLD4InvMassCut,
					float(0.0)
				);

	registerProcessorParameter(	"CSLDNeutralSLD5InvMassCut",
					"Cut on Invariant mass of visible neutral decay products for semi-leptonic decay of C-Hadron when lepton is with third vertex",
					m_CSLDNeutralSLD5InvMassCut,
					float(0.0)
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
		m_pTTree1->Branch("SLDFlavour", &m_SLDFlavour);
		m_pTTree1->Branch("SLDType", &m_SLDType);
		m_pTTree1->Branch("SLDLeptonID", &m_SLDLeptonID);
		m_pTTree1->Branch("leptonE_to_parentE", &m_leptonE_to_parentE);
		m_pTTree1->Branch("otherChargedE_to_parentE", &m_otherChargedE_to_parentE);
		m_pTTree1->Branch("allChargedE_to_parentE", &m_allChargedE_to_parentE);
		m_pTTree1->Branch("neutralE_to_parentE", &m_neutralE_to_parentE);
		m_pTTree1->Branch("neutrino_to_parentE", &m_neutrino_to_parentE);
		m_pTTree1->Branch("nSLDecayOfBHadron",&m_nSLDecayOfBHadron,"nSLDecayOfBHadron/I");
		m_pTTree1->Branch("nSLDecayOfCHadron",&m_nSLDecayOfCHadron,"nSLDecayOfCHadron/I");
		m_pTTree1->Branch("nSLDecayOfTauLepton",&m_nSLDecayOfTauLepton,"nSLDecayOfTauLepton/I");
		m_pTTree1->Branch("nSLDecayTotal",&m_nSLDecayTotal,"nSLDecayTotal/I");
		m_pTTree1->Branch("nSLDecayToElectron",&m_nSLDecayToElectron,"nSLDecayToElectron/I");
		m_pTTree1->Branch("nSLDecayToMuon",&m_nSLDecayToMuon,"nSLDecayToMuon/I");
		m_pTTree1->Branch("nSLDecayToTau",&m_nSLDecayToTau,"nSLDecayToTau/I");
		m_pTTree1->Branch("SLDStatus", &m_SLDStatus );
		m_pTTree1->Branch("nTauNeutrino",&m_nTauNeutrino,"nTauNeutrino/I");
		m_pTTree1->Branch("nNeutrino",&m_nNeutrino,"nNeutrino/I");
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
		m_pTTree1->Branch("true_E_vis", &m_true_E_vis);
		m_pTTree1->Branch("true_E_vis_prime", &m_true_E_vis_prime);
		m_pTTree1->Branch("true_P_vis_par", &m_true_P_vis_par);
		m_pTTree1->Branch("true_P_vis_par_prime", &m_true_P_vis_par_prime);
		m_pTTree1->Branch("true_P_vis_nor", &m_true_P_vis_nor);
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
		m_pTTree1->Branch("flightDirectionErrorDeltaTheta" , &m_FlightDirectionErrorDeltaTheta);
		m_pTTree1->Branch("flightDirectionErrorDeltaPhi" , &m_FlightDirectionErrorDeltaPhi);
		m_pTTree1->Branch("distRecoLeptonToDownStreamVertex", &m_distRecoLeptonToDownStreamVertex);
		m_pTTree1->Branch("dsVertexResidualX", &m_dsVertexResidualX);
		m_pTTree1->Branch("dsVertexResidualY", &m_dsVertexResidualY);
		m_pTTree1->Branch("dsVertexResidualZ", &m_dsVertexResidualZ);
		m_pTTree1->Branch("secVertexResidualX", &m_SecVertexResidualX);
		m_pTTree1->Branch("secVertexResidualY", &m_SecVertexResidualY);
		m_pTTree1->Branch("secVertexResidualZ", &m_SecVertexResidualZ);
		m_pTTree1->Branch("parentHadronMass", &m_parentHadronMass );
		m_pTTree1->Branch("parentHadronPDG", &m_parentHadronPDG );
		m_pTTree1->Branch("trueParentHadronFlightDistance", &m_trueParentHadronFlightDistance );
		m_pTTree1->Branch("recoParentHadronFlightDistance", &m_recoParentHadronFlightDistance );
		m_pTTree1->Branch("daughterHadronMass", &m_daughterHadronMass );
		m_pTTree1->Branch("daughterHadronPDG", &m_daughterHadronPDG );
		m_pTTree1->Branch("daughterHadronFlightDistance", &m_daughterHadronFlightDistance );
		m_pTTree1->Branch("jetEnergy", &m_jetEnergy );
		m_pTTree1->Branch("jetEnergyFractionCharged", &m_jetEnergyFractionCharged );
		m_pTTree1->Branch("jetEnergyFractionNeutralHadron", &m_jetEnergyFractionNeutralHadron );
		m_pTTree1->Branch("jetEnergyFractionPhoton", &m_jetEnergyFractionPhoton );
		m_pTTree1->Branch("jetEnergyFractionNeutrals", &m_jetEnergyFractionNeutrals );
		m_pTTree1->Branch("nTrueNeutralDecayProducts",&m_nTrueNeutralDecayProducts);
		m_pTTree1->Branch("nTrueAloneChargedDecayProducts",&m_nTrueAloneChargedDecayProducts);
		m_pTTree1->Branch("nTrueVertices",&m_nTrueVertices);
		m_pTTree1->Branch("nJetsVerticesDistributedIn",&m_nJetsVerticesDistributedIn);
		m_pTTree1->Branch("nRecoVerticesInJet",&m_nRecoVerticesInJet);
		m_pTTree1->Branch("nAloneChargedPFOs",&m_nAloneChargedPFOs);
		m_pTTree1->Branch("nAloneChargedPFOsFromSLD",&m_nAloneChargedPFOsFromSLD);
		m_pTTree1->Branch("distLeptonAlonePFOsNotFromSLD",&m_distLeptonAlonePFOsNotFromSLD);
		m_pTTree1->Branch("distLeptonAlonePFOsFromSLD",&m_distLeptonAlonePFOsFromSLD);
		m_pTTree1->Branch("alphaTrueNeutrals", &m_alphaTrueNeutrals );
		m_pTTree1->Branch("alphaTrueCharged", &m_alphaTrueCharged );
		m_pTTree1->Branch("cosAlphaNeutrals", &m_cosAlphaNeutrals );
		m_pTTree1->Branch("cosAlphaAloneCharged", &m_cosAlphaCharged );
		m_pTTree1->Branch("alphaChargedPFOsFromSLDwrtLepton", &m_alphaChargedPFOsFromSLDwrtLepton );
		m_pTTree1->Branch("cosAlphaChargedPFOsFromSLDwrtLepton", &m_cosAlphaChargedPFOsFromSLDwrtLepton );
		m_pTTree1->Branch("alphaChargedPFOsFromSLDwrtFD", &m_alphaChargedPFOsFromSLDwrtFD );
		m_pTTree1->Branch("cosAlphaChargedPFOsFromSLDwrtFD", &m_cosAlphaChargedPFOsFromSLDwrtFD );
		m_pTTree1->Branch("alphaChargedPFOsFromSLDwrtJet", &m_alphaChargedPFOsFromSLDwrtJet );
		m_pTTree1->Branch("cosAlphaChargedPFOsFromSLDwrtJet", &m_cosAlphaChargedPFOsFromSLDwrtJet );
		m_pTTree1->Branch("alphaChargedPFOsNotFromSLDwrtLepton", &m_alphaChargedPFOsNotFromSLDwrtLepton );
		m_pTTree1->Branch("cosAlphaChargedPFOsNotFromSLDwrtLepton", &m_cosAlphaChargedPFOsNotFromSLDwrtLepton );
		m_pTTree1->Branch("alphaChargedPFOsNotFromSLDwrtFD", &m_alphaChargedPFOsNotFromSLDwrtFD );
		m_pTTree1->Branch("cosAlphaChargedPFOsNotFromSLDwrtFD", &m_cosAlphaChargedPFOsNotFromSLDwrtFD );
		m_pTTree1->Branch("alphaChargedPFOsNotFromSLDwrtJet", &m_alphaChargedPFOsNotFromSLDwrtJet );
		m_pTTree1->Branch("cosAlphaChargedPFOsNotFromSLDwrtJet", &m_cosAlphaChargedPFOsNotFromSLDwrtJet );
		m_pTTree1->Branch("alphaAloneChargedPFOsFromSLDwrtLepton", &m_alphaAloneChargedPFOsFromSLDwrtLepton );
		m_pTTree1->Branch("cosAlphaAloneChargedPFOsFromSLDwrtLepton", &m_cosAlphaAloneChargedPFOsFromSLDwrtLepton );
		m_pTTree1->Branch("alphaAloneChargedPFOsFromSLDwrtFD", &m_alphaAloneChargedPFOsFromSLDwrtFD );
		m_pTTree1->Branch("cosAlphaAloneChargedPFOsFromSLDwrtFD", &m_cosAlphaAloneChargedPFOsFromSLDwrtFD );
		m_pTTree1->Branch("alphaAloneChargedPFOsFromSLDwrtJet", &m_alphaAloneChargedPFOsFromSLDwrtJet );
		m_pTTree1->Branch("cosAlphaAloneChargedPFOsFromSLDwrtJet", &m_cosAlphaAloneChargedPFOsFromSLDwrtJet );
		m_pTTree1->Branch("alphaAloneChargedPFOsNotFromSLDwrtLepton", &m_alphaAloneChargedPFOsNotFromSLDwrtLepton );
		m_pTTree1->Branch("cosAlphaAloneChargedPFOsNotFromSLDwrtLepton", &m_cosAlphaAloneChargedPFOsNotFromSLDwrtLepton );
		m_pTTree1->Branch("alphaAloneChargedPFOsNotFromSLDwrtFD", &m_alphaAloneChargedPFOsNotFromSLDwrtFD );
		m_pTTree1->Branch("cosAlphaAloneChargedPFOsNotFromSLDwrtFD", &m_cosAlphaAloneChargedPFOsNotFromSLDwrtFD );
		m_pTTree1->Branch("alphaAloneChargedPFOsNotFromSLDwrtJet", &m_alphaAloneChargedPFOsNotFromSLDwrtJet );
		m_pTTree1->Branch("cosAlphaAloneChargedPFOsNotFromSLDwrtJet", &m_cosAlphaAloneChargedPFOsNotFromSLDwrtJet );
		m_pTTree1->Branch("alphaNeutralPFOsFromSLDwrtLepton", &m_alphaNeutralPFOsFromSLDwrtLepton );
		m_pTTree1->Branch("cosAlphaNeutralPFOsFromSLDwrtLepton", &m_cosAlphaNeutralPFOsFromSLDwrtLepton );
		m_pTTree1->Branch("alphaNeutralPFOsFromSLDwrtFD", &m_alphaNeutralPFOsFromSLDwrtFD );
		m_pTTree1->Branch("cosAlphaNeutralPFOsFromSLDwrtFD", &m_cosAlphaNeutralPFOsFromSLDwrtFD );
		m_pTTree1->Branch("alphaNeutralPFOsFromSLDwrtJet", &m_alphaNeutralPFOsFromSLDwrtJet );
		m_pTTree1->Branch("cosAlphaNeutralPFOsFromSLDwrtJet", &m_cosAlphaNeutralPFOsFromSLDwrtJet );
		m_pTTree1->Branch("alphaNeutralPFOsNotFromSLDwrtLepton", &m_alphaNeutralPFOsNotFromSLDwrtLepton );
		m_pTTree1->Branch("cosAlphaNeutralPFOsNotFromSLDwrtLepton", &m_cosAlphaNeutralPFOsNotFromSLDwrtLepton );
		m_pTTree1->Branch("alphaNeutralPFOsNotFromSLDwrtFD", &m_alphaNeutralPFOsNotFromSLDwrtFD );
		m_pTTree1->Branch("cosAlphaNeutralPFOsNotFromSLDwrtFD", &m_cosAlphaNeutralPFOsNotFromSLDwrtFD );
		m_pTTree1->Branch("alphaNeutralPFOsNotFromSLDwrtJet", &m_alphaNeutralPFOsNotFromSLDwrtJet );
		m_pTTree1->Branch("cosAlphaNeutralPFOsNotFromSLDwrtJet", &m_cosAlphaNeutralPFOsNotFromSLDwrtJet );
		m_pTTree1->Branch("weightPFOtoMCP_Lepton", &m_weightPFOtoMCP_Lepton );
		m_pTTree1->Branch("weightMCPtoPFO_Lepton", &m_weightMCPtoPFO_Lepton );
		m_pTTree1->Branch("weightPFOtoMCP_Neutral", &m_weightPFOtoMCP_Neutral );
		m_pTTree1->Branch("weightMCPtoPFO_Neutral", &m_weightMCPtoPFO_Neutral );
		m_pTTree1->Branch("weightPFOtoMCP_Charged", &m_weightPFOtoMCP_Charged );
		m_pTTree1->Branch("weightMCPtoPFO_Charged", &m_weightMCPtoPFO_Charged );

		m_pTTree2 = new TTree("FourMomentums", "FourMomentums");
		m_pTTree2->Branch("SLDStatus", &m_SLDStatus );
		m_pTTree2->Branch("SLDFlavour", &m_SLDFlavour);
		m_pTTree2->Branch("SLDType", &m_SLDType);
		m_pTTree2->Branch("SLDLeptonID", &m_SLDLeptonID);
		m_pTTree2->Branch("leptonE_to_parentE", &m_leptonE_to_parentE);
		m_pTTree2->Branch("otherChargedE_to_parentE", &m_otherChargedE_to_parentE);
		m_pTTree2->Branch("allChargedE_to_parentE", &m_allChargedE_to_parentE);
		m_pTTree2->Branch("neutralE_to_parentE", &m_neutralE_to_parentE);
		m_pTTree2->Branch("neutrino_to_parentE", &m_neutrino_to_parentE);
		m_pTTree2->Branch("visibleChargedInvMassCut", &m_visibleChargedInvMassCut );
		m_pTTree2->Branch("visibleNeutralInvMassCut", &m_visibleNeutralInvMassCut );
		m_pTTree2->Branch("trueNeutralPx", &m_trueNeutralPx );
		m_pTTree2->Branch("trueNeutralPy", &m_trueNeutralPy );
		m_pTTree2->Branch("trueNeutralPz", &m_trueNeutralPz );
		m_pTTree2->Branch("trueNeutralE", &m_trueNeutralE );
		m_pTTree2->Branch("trueNeutralM", &m_trueNeutralM );
		m_pTTree2->Branch("trueChargedPx", &m_trueChargedPx );
		m_pTTree2->Branch("trueChargedPy", &m_trueChargedPy );
		m_pTTree2->Branch("trueChargedPz", &m_trueChargedPz );
		m_pTTree2->Branch("trueChargedE", &m_trueChargedE );
		m_pTTree2->Branch("trueChargedM", &m_trueChargedM );
		m_pTTree2->Branch("trueLeptonPx", &m_trueLeptonPx );
		m_pTTree2->Branch("trueLeptonPy", &m_trueLeptonPy );
		m_pTTree2->Branch("trueLeptonPz", &m_trueLeptonPz );
		m_pTTree2->Branch("trueLeptonE", &m_trueLeptonE );
		m_pTTree2->Branch("trueVisiblePx", &m_trueVisiblePx );
		m_pTTree2->Branch("trueVisiblePy", &m_trueVisiblePy );
		m_pTTree2->Branch("trueVisiblePz", &m_trueVisiblePz );
		m_pTTree2->Branch("trueVisibleE", &m_trueVisibleE );
		m_pTTree2->Branch("trueVisibleM", &m_trueVisibleM );
		m_pTTree2->Branch("trueNeutrinoPx", &m_trueNeutrinoPx );
		m_pTTree2->Branch("trueNeutrinoPy", &m_trueNeutrinoPy );
		m_pTTree2->Branch("trueNeutrinoPz", &m_trueNeutrinoPz );
		m_pTTree2->Branch("trueNeutrinoE", &m_trueNeutrinoE );
		m_pTTree2->Branch("trueHadronPx", &m_trueHadronPx );
		m_pTTree2->Branch("trueHadronPy", &m_trueHadronPy );
		m_pTTree2->Branch("trueHadronPz", &m_trueHadronPz );
		m_pTTree2->Branch("trueHadronE", &m_trueHadronE );
		m_pTTree2->Branch("trueHadronM", &m_trueHadronM );
		m_pTTree2->Branch("recoNeutralPx", &m_recoNeutralPx );
		m_pTTree2->Branch("recoNeutralPy", &m_recoNeutralPy );
		m_pTTree2->Branch("recoNeutralPz", &m_recoNeutralPz );
		m_pTTree2->Branch("recoNeutralE", &m_recoNeutralE );
		m_pTTree2->Branch("recoNeutralM", &m_recoNeutralM );
		m_pTTree2->Branch("recoChargedPx", &m_recoChargedPx );
		m_pTTree2->Branch("recoChargedPy", &m_recoChargedPy );
		m_pTTree2->Branch("recoChargedPz", &m_recoChargedPz );
		m_pTTree2->Branch("recoChargedE", &m_recoChargedE );
		m_pTTree2->Branch("recoChargedM", &m_recoChargedM );
		m_pTTree2->Branch("recoLeptonPx", &m_recoLeptonPx );
		m_pTTree2->Branch("recoLeptonPy", &m_recoLeptonPy );
		m_pTTree2->Branch("recoLeptonPz", &m_recoLeptonPz );
		m_pTTree2->Branch("recoLeptonE", &m_recoLeptonE );
		m_pTTree2->Branch("recoVisiblePx", &m_recoVisiblePx );
		m_pTTree2->Branch("recoVisiblePy", &m_recoVisiblePy );
		m_pTTree2->Branch("recoVisiblePz", &m_recoVisiblePz );
		m_pTTree2->Branch("recoVisibleE", &m_recoVisibleE );
		m_pTTree2->Branch("recoVisibleM", &m_recoVisibleM );
		m_pTTree2->Branch("recoNeutrinoPx", &m_recoNeutrinoPx );
		m_pTTree2->Branch("recoNeutrinoPy", &m_recoNeutrinoPy );
		m_pTTree2->Branch("recoNeutrinoPz", &m_recoNeutrinoPz );
		m_pTTree2->Branch("recoNeutrinoE", &m_recoNeutrinoE );
		m_pTTree2->Branch("recoHadronPx", &m_recoHadronPx );
		m_pTTree2->Branch("recoHadronPy", &m_recoHadronPy );
		m_pTTree2->Branch("recoHadronPz", &m_recoHadronPz );
		m_pTTree2->Branch("recoHadronE", &m_recoHadronE );
		m_pTTree2->Branch("recoHadronM", &m_recoHadronM );
		m_pTTree2->Branch("usedNeutralPx", &m_usedNeutralPx );
		m_pTTree2->Branch("usedNeutralPy", &m_usedNeutralPy );
		m_pTTree2->Branch("usedNeutralPz", &m_usedNeutralPz );
		m_pTTree2->Branch("usedNeutralE", &m_usedNeutralE );
		m_pTTree2->Branch("usedNeutralM", &m_usedNeutralM );
		m_pTTree2->Branch("usedChargedPx", &m_usedChargedPx );
		m_pTTree2->Branch("usedChargedPy", &m_usedChargedPy );
		m_pTTree2->Branch("usedChargedPz", &m_usedChargedPz );
		m_pTTree2->Branch("usedChargedE", &m_usedChargedE );
		m_pTTree2->Branch("usedChargedM", &m_usedChargedM );
		m_pTTree2->Branch("usedLeptonPx", &m_usedLeptonPx );
		m_pTTree2->Branch("usedLeptonPy", &m_usedLeptonPy );
		m_pTTree2->Branch("usedLeptonPz", &m_usedLeptonPz );
		m_pTTree2->Branch("usedLeptonE", &m_usedLeptonE );
		m_pTTree2->Branch("usedVisiblePx", &m_usedVisiblePx );
		m_pTTree2->Branch("usedVisiblePy", &m_usedVisiblePy );
		m_pTTree2->Branch("usedVisiblePz", &m_usedVisiblePz );
		m_pTTree2->Branch("usedVisibleE", &m_usedVisibleE );
		m_pTTree2->Branch("usedVisibleM", &m_usedVisibleM );
		m_pTTree2->Branch("cheatedPVARecoNeutralPx", &m_cheatedPVARecoNeutralPx );
		m_pTTree2->Branch("cheatedPVARecoNeutralPy", &m_cheatedPVARecoNeutralPy );
		m_pTTree2->Branch("cheatedPVARecoNeutralPz", &m_cheatedPVARecoNeutralPz );
		m_pTTree2->Branch("cheatedPVARecoNeutralE", &m_cheatedPVARecoNeutralE );
		m_pTTree2->Branch("cheatedPVARecoNeutralM", &m_cheatedPVARecoNeutralM );
		m_pTTree2->Branch("cheatedPVARecoChargedPx", &m_cheatedPVARecoChargedPx );
		m_pTTree2->Branch("cheatedPVARecoChargedPy", &m_cheatedPVARecoChargedPy );
		m_pTTree2->Branch("cheatedPVARecoChargedPz", &m_cheatedPVARecoChargedPz );
		m_pTTree2->Branch("cheatedPVARecoChargedE", &m_cheatedPVARecoChargedE );
		m_pTTree2->Branch("cheatedPVARecoChargedM", &m_cheatedPVARecoChargedM );
		m_pTTree2->Branch("neutralEnergy", &m_neutralEnergy );
		m_pTTree2->Branch("neutralEnergyFromVertexing", &m_neutralEnergyFromVertexing );
		m_pTTree2->Branch("neutralEnergyFromPVA", &m_neutralEnergyFromPVA );
		m_pTTree2->Branch("neutralMomentum", &m_neutralMomentum );
		m_pTTree2->Branch("neutralMomentumFromVertexing", &m_neutralMomentumFromVertexing );
		m_pTTree2->Branch("neutralMomentumFromPVA", &m_neutralMomentumFromPVA );
		m_pTTree2->Branch("chargedEnergy", &m_chargedEnergy );
		m_pTTree2->Branch("chargedEnergyFromVertexing", &m_chargedEnergyFromVertexing );
		m_pTTree2->Branch("chargedEnergyFromPVA", &m_chargedEnergyFromPVA );
		m_pTTree2->Branch("chargedMomentum", &m_chargedMomentum );
		m_pTTree2->Branch("chargedMomentumFromVertexing", &m_chargedMomentumFromVertexing );
		m_pTTree2->Branch("chargedMomentumFromPVA", &m_chargedMomentumFromPVA );
		m_pTTree2->Branch("expectedNeutralEnergy", &m_expectedNeutralEnergy );
		m_pTTree2->Branch("expectedNeutralEnergyFromVertexing", &m_expectedNeutralEnergyFromVertexing );
		m_pTTree2->Branch("expectedNeutralEnergyFromPVA", &m_expectedNeutralEnergyFromPVA );
		m_pTTree2->Branch("expectedNeutralMomentum", &m_expectedNeutralMomentum );
		m_pTTree2->Branch("expectedNeutralMomentumFromVertexing", &m_expectedNeutralMomentumFromVertexing );
		m_pTTree2->Branch("expectedNeutralMomentumFromPVA", &m_expectedNeutralMomentumFromPVA );
		m_pTTree2->Branch("expectedChargedEnergy", &m_expectedChargedEnergy );
		m_pTTree2->Branch("expectedChargedEnergyFromVertexing", &m_expectedChargedEnergyFromVertexing );
		m_pTTree2->Branch("expectedChargedEnergyFromPVA", &m_expectedChargedEnergyFromPVA );
		m_pTTree2->Branch("expectedChargedMomentum", &m_expectedChargedMomentum );
		m_pTTree2->Branch("expectedChargedMomentumFromVertexing", &m_expectedChargedMomentumFromVertexing );
		m_pTTree2->Branch("expectedChargedMomentumFromPVA", &m_expectedChargedMomentumFromPVA );
		m_pTTree2->Branch("NuPxResidual", &m_NuPxResidual);
		m_pTTree2->Branch("NuPyResidual", &m_NuPyResidual);
		m_pTTree2->Branch("NuPzResidual", &m_NuPzResidual);
		m_pTTree2->Branch("NuEResidual", &m_NuEResidual);



		h_SLDStatus = new TH1I( "SLDStatus" , ";" , 7 , 0 , 7 ); n_SLDStatus = 0;
		h_SLDStatus->GetXaxis()->SetBinLabel(1,"No l^{REC}");
		h_SLDStatus->GetXaxis()->SetBinLabel(2,"lep not in jet");
		h_SLDStatus->GetXaxis()->SetBinLabel(3,"lep in Prim. Vtx");
		h_SLDStatus->GetXaxis()->SetBinLabel(4,"lep in Sec. Vtx");
		h_SLDStatus->GetXaxis()->SetBinLabel(5,"lep + 3^{rd} Vtx");
		h_SLDStatus->GetXaxis()->SetBinLabel(6,"lep + alone track");
		h_SLDStatus->GetXaxis()->SetBinLabel(7,"other");
//		h_SLDStatus->GetXaxis()->SetBinLabel(1,"lep not found");
//		h_SLDStatus->GetXaxis()->SetBinLabel(2,"lep not in jet");
//		h_SLDStatus->GetXaxis()->SetBinLabel(3,"lep in Prim. Vtx");
//		h_SLDStatus->GetXaxis()->SetBinLabel(4,"lep in Sec. Vtx");
//		h_SLDStatus->GetXaxis()->SetBinLabel(5,"lep + 3^{rd} Vtx");
//		h_SLDStatus->GetXaxis()->SetBinLabel(6,"lep + alone track");
//		h_SLDStatus->GetXaxis()->SetBinLabel(7,"other");

		h_BHadronType = new TH1F( "BHadronType" , ";" , 8 , -0.5 , 7.5 );
		h_BHadronType->GetXaxis()->SetBinLabel(1,"B^{0}");//PDG = 511
		h_BHadronType->GetXaxis()->SetBinLabel(2,"B^{#pm}");//PDG = 521
		h_BHadronType->GetXaxis()->SetBinLabel(3,"B^{0}_{s}");//PDG = 531
		h_BHadronType->GetXaxis()->SetBinLabel(4,"B^{#pm}_{c}");//PDG = 541
		h_BHadronType->GetXaxis()->SetBinLabel(5,"#Lambda^{0}_{b}");//PDG = 5122
		h_BHadronType->GetXaxis()->SetBinLabel(6,"#Xi^{#pm}_{b}");//PDG = 5132
		h_BHadronType->GetXaxis()->SetBinLabel(7,"#Xi^{0}_{b}");//PDG = 5232
		h_BHadronType->GetXaxis()->SetBinLabel(8,"#Omega^{-}_{b}");//PDG = 5332

		h_CHadronType = new TH1F( "CHadronType" , ";" , 16 , -0.5 , 15.5 );
		h_CHadronType->GetXaxis()->SetBinLabel(1,"D^{#pm}");//PDG = 411
		h_CHadronType->GetXaxis()->SetBinLabel(2,"D^{*}(2010)^{#pm}");//PDG = 413
		h_CHadronType->GetXaxis()->SetBinLabel(3,"D^{*}_{2}(2460)^{#pm}");//PDG = 415
		h_CHadronType->GetXaxis()->SetBinLabel(4,"D^{0}");//PDG = 421
		h_CHadronType->GetXaxis()->SetBinLabel(5,"D^{*}(2007)^{0}");//PDG = 423
		h_CHadronType->GetXaxis()->SetBinLabel(6,"D^{*}_{2}(2460)^{0}");//PDG = 425
		h_CHadronType->GetXaxis()->SetBinLabel(7,"D^{#pm}_{s}");//PDG = 431
		h_CHadronType->GetXaxis()->SetBinLabel(8,"D^{*#pm}_{s}");//PDG = 433
		h_CHadronType->GetXaxis()->SetBinLabel(9,"D^{*}_{s2}");//PDG = 435
		h_CHadronType->GetXaxis()->SetBinLabel(10,"#eta_{c}");//PDG = 441
		h_CHadronType->GetXaxis()->SetBinLabel(11,"J/#psi");//PDG = 443
		h_CHadronType->GetXaxis()->SetBinLabel(12,"#Lambda^{#pm}_{c}");//PDG = 4122
		h_CHadronType->GetXaxis()->SetBinLabel(13,"#Xi^{0}_{c}");//PDG = 4132
		h_CHadronType->GetXaxis()->SetBinLabel(14,"#Xi^{#pm}_{c}");//PDG = 4232
		h_CHadronType->GetXaxis()->SetBinLabel(15,"#Omega^{0}_{c}");//PDG = 4332
		h_CHadronType->GetXaxis()->SetBinLabel(16,"#Omega^{*0}_{c}");//PDG = 4334

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
		h_SLDecayFlavour = new TH1I( "SLDecayFlavour" , "; SLDecay Type" , 3 , 0 , 3 );
		h_SLDecayFlavour->GetXaxis()->SetBinLabel(1,"SLD of B-Hadron");
		h_SLDecayFlavour->GetXaxis()->SetBinLabel(2,"SLD of C-Hadron");
		h_SLDecayFlavour->GetXaxis()->SetBinLabel(3,"LD of Tau-Lepton");
		h_SLDecayModeB = new TH1I( "SLDMode B-Hadron" , "; SLDecay Mode" , 3 , 0 , 3 );
		h_SLDecayModeB->GetXaxis()->SetBinLabel(1,"X#rightarrow e#nu_{e}Y");
		h_SLDecayModeB->GetXaxis()->SetBinLabel(2,"X#rightarrow #mu#nu_{#mu}Y");
		h_SLDecayModeB->GetXaxis()->SetBinLabel(3,"X#rightarrow #tau#nu_{#tau}Y");
		h_SLDecayModeC = new TH1I( "SLDMode C-Hadron" , "; SLDecay Mode" , 3 , 0 , 3 );
		h_SLDecayModeC->GetXaxis()->SetBinLabel(1,"X#rightarrow e#nu_{e}Y");
		h_SLDecayModeC->GetXaxis()->SetBinLabel(2,"X#rightarrow #mu#nu_{#mu}Y");
		h_SLDecayModeC->GetXaxis()->SetBinLabel(3,"X#rightarrow #tau#nu_{#tau}Y");
		h_SLDecayOrder = new TH1I( "SLDecayOrder" , "; SLDecay Type" , 2 , 0 , 2 );
		h_SLDecayOrder->GetXaxis()->SetBinLabel(1,"with UpStream/DownStream SLD");
		h_SLDecayOrder->GetXaxis()->SetBinLabel(2,"without UpStream/DownStream SLD");
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

	BHadPDGs = { 511 , 521 , 531 , 541 , 5122 , 5132 , 5232 , 5332 };
	CHadPDGs = { 411 , 413 , 415 , 421 , 423 , 425 , 431 , 433 , 435 , 441 , 443 , 4122 , 4132 , 4232 , 4332 , 4334 };
}

void SLDCorrection::Clear()
{
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
	m_SLDFlavour.clear();
	m_SLDType.clear();
	m_SLDLeptonID.clear();
	m_leptonE_to_parentE.clear();
	m_otherChargedE_to_parentE.clear();
	m_allChargedE_to_parentE.clear();
	m_neutralE_to_parentE.clear();
	m_neutrino_to_parentE.clear();
	m_nSLDecayOfBHadron = 0;
	m_nSLDecayOfCHadron = 0;
	m_nSLDecayOfTauLepton = 0;
	m_nSLDecayTotal = 0;
	m_nSLDecayToElectron = 0;
	m_nSLDecayToMuon = 0;
	m_nSLDecayToTau = 0;
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
	m_true_E_vis.clear();
	m_true_E_vis_prime.clear();
	m_true_P_vis_par.clear();
	m_true_P_vis_par_prime.clear();
	m_true_P_vis_nor.clear();
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
	m_FlightDirectionErrorDeltaTheta.clear();
	m_FlightDirectionErrorDeltaPhi.clear();
	m_dsVertexResidualX.clear();
	m_dsVertexResidualY.clear();
	m_dsVertexResidualZ.clear();
	m_SecVertexResidualX.clear();
	m_SecVertexResidualY.clear();
	m_SecVertexResidualZ.clear();
	m_parentHadronMass.clear();
	m_parentHadronPDG.clear();
	m_trueParentHadronFlightDistance.clear();
	m_recoParentHadronFlightDistance.clear();
	m_daughterHadronMass.clear();
	m_daughterHadronPDG.clear();
	m_daughterHadronFlightDistance.clear();
	m_jetEnergy.clear();
	m_jetEnergyFractionCharged.clear();
	m_jetEnergyFractionNeutralHadron.clear();
	m_jetEnergyFractionPhoton.clear();
	m_jetEnergyFractionNeutrals.clear();
	m_SLDStatus.clear();
	m_alphaTrueNeutrals.clear();
	m_alphaTrueCharged.clear();
	m_cosAlphaNeutrals.clear();
	m_cosAlphaCharged.clear();
	m_alphaChargedPFOsFromSLDwrtLepton.clear();
	m_cosAlphaChargedPFOsFromSLDwrtLepton.clear();
	m_alphaChargedPFOsFromSLDwrtFD.clear();
	m_cosAlphaChargedPFOsFromSLDwrtFD.clear();
	m_alphaChargedPFOsFromSLDwrtJet.clear();
	m_cosAlphaChargedPFOsFromSLDwrtJet.clear();
	m_alphaChargedPFOsNotFromSLDwrtLepton.clear();
	m_cosAlphaChargedPFOsNotFromSLDwrtLepton.clear();
	m_alphaChargedPFOsNotFromSLDwrtFD.clear();
	m_cosAlphaChargedPFOsNotFromSLDwrtFD.clear();
	m_alphaChargedPFOsNotFromSLDwrtJet.clear();
	m_cosAlphaChargedPFOsNotFromSLDwrtJet.clear();
	m_alphaAloneChargedPFOsFromSLDwrtLepton.clear();
	m_cosAlphaAloneChargedPFOsFromSLDwrtLepton.clear();
	m_alphaAloneChargedPFOsFromSLDwrtFD.clear();
	m_cosAlphaAloneChargedPFOsFromSLDwrtFD.clear();
	m_alphaAloneChargedPFOsFromSLDwrtJet.clear();
	m_cosAlphaAloneChargedPFOsFromSLDwrtJet.clear();
	m_alphaAloneChargedPFOsNotFromSLDwrtLepton.clear();
	m_cosAlphaAloneChargedPFOsNotFromSLDwrtLepton.clear();
	m_alphaAloneChargedPFOsNotFromSLDwrtFD.clear();
	m_cosAlphaAloneChargedPFOsNotFromSLDwrtFD.clear();
	m_alphaAloneChargedPFOsNotFromSLDwrtJet.clear();
	m_cosAlphaAloneChargedPFOsNotFromSLDwrtJet.clear();
	m_alphaNeutralPFOsFromSLDwrtLepton.clear();
	m_cosAlphaNeutralPFOsFromSLDwrtLepton.clear();
	m_alphaNeutralPFOsFromSLDwrtFD.clear();
	m_cosAlphaNeutralPFOsFromSLDwrtFD.clear();
	m_alphaNeutralPFOsFromSLDwrtJet.clear();
	m_cosAlphaNeutralPFOsFromSLDwrtJet.clear();
	m_alphaNeutralPFOsNotFromSLDwrtLepton.clear();
	m_cosAlphaNeutralPFOsNotFromSLDwrtLepton.clear();
	m_alphaNeutralPFOsNotFromSLDwrtFD.clear();
	m_cosAlphaNeutralPFOsNotFromSLDwrtFD.clear();
	m_alphaNeutralPFOsNotFromSLDwrtJet.clear();
	m_cosAlphaNeutralPFOsNotFromSLDwrtJet.clear();
	m_weightPFOtoMCP_Lepton.clear();
	m_weightMCPtoPFO_Lepton.clear();
	m_weightPFOtoMCP_Neutral.clear();
	m_weightMCPtoPFO_Neutral.clear();
	m_weightPFOtoMCP_Charged.clear();
	m_weightMCPtoPFO_Charged.clear();
	m_distRecoLeptonToDownStreamVertex.clear();

	m_visibleChargedInvMassCut.clear();
	m_visibleNeutralInvMassCut.clear();
	m_trueNeutralPx.clear();
	m_trueNeutralPy.clear();
	m_trueNeutralPz.clear();
	m_trueNeutralE.clear();
	m_trueNeutralM.clear();
	m_trueChargedPx.clear();
	m_trueChargedPy.clear();
	m_trueChargedPz.clear();
	m_trueChargedE.clear();
	m_trueChargedM.clear();
	m_trueLeptonPx.clear();
	m_trueLeptonPy.clear();
	m_trueLeptonPz.clear();
	m_trueLeptonE.clear();
	m_trueVisiblePx.clear();
	m_trueVisiblePy.clear();
	m_trueVisiblePz.clear();
	m_trueVisibleE.clear();
	m_trueVisibleM.clear();
	m_trueNeutrinoPx.clear();
	m_trueNeutrinoPy.clear();
	m_trueNeutrinoPz.clear();
	m_trueNeutrinoE.clear();
	m_trueHadronPx.clear();
	m_trueHadronPy.clear();
	m_trueHadronPz.clear();
	m_trueHadronE.clear();
	m_trueHadronM.clear();
	m_recoNeutralPx.clear();
	m_recoNeutralPy.clear();
	m_recoNeutralPz.clear();
	m_recoNeutralE.clear();
	m_recoNeutralM.clear();
	m_recoChargedPx.clear();
	m_recoChargedPy.clear();
	m_recoChargedPz.clear();
	m_recoChargedE.clear();
	m_recoChargedM.clear();
	m_recoLeptonPx.clear();
	m_recoLeptonPy.clear();
	m_recoLeptonPz.clear();
	m_recoLeptonE.clear();
	m_recoVisiblePx.clear();
	m_recoVisiblePy.clear();
	m_recoVisiblePz.clear();
	m_recoVisibleE.clear();
	m_recoVisibleM.clear();
	m_recoNeutrinoPx.clear();
	m_recoNeutrinoPy.clear();
	m_recoNeutrinoPz.clear();
	m_recoNeutrinoE.clear();
	m_recoHadronPx.clear();
	m_recoHadronPy.clear();
	m_recoHadronPz.clear();
	m_recoHadronE.clear();
	m_recoHadronM.clear();
	m_usedNeutralPx.clear();
	m_usedNeutralPy.clear();
	m_usedNeutralPz.clear();
	m_usedNeutralE.clear();
	m_usedNeutralM.clear();
	m_usedChargedPx.clear();
	m_usedChargedPy.clear();
	m_usedChargedPz.clear();
	m_usedChargedE.clear();
	m_usedChargedM.clear();
	m_usedLeptonPx.clear();
	m_usedLeptonPy.clear();
	m_usedLeptonPz.clear();
	m_usedLeptonE.clear();
	m_usedVisiblePx.clear();
	m_usedVisiblePy.clear();
	m_usedVisiblePz.clear();
	m_usedVisibleE.clear();
	m_usedVisibleM.clear();
	m_cheatedPVARecoNeutralPx.clear();
	m_cheatedPVARecoNeutralPy.clear();
	m_cheatedPVARecoNeutralPz.clear();
	m_cheatedPVARecoNeutralE.clear();
	m_cheatedPVARecoNeutralM.clear();
	m_cheatedPVARecoChargedPx.clear();
	m_cheatedPVARecoChargedPy.clear();
	m_cheatedPVARecoChargedPz.clear();
	m_cheatedPVARecoChargedE.clear();
	m_cheatedPVARecoChargedM.clear();
	m_neutralEnergy.clear();
	m_neutralEnergyFromVertexing.clear();
	m_neutralEnergyFromPVA.clear();
	m_neutralMomentum.clear();
	m_neutralMomentumFromVertexing.clear();
	m_neutralMomentumFromPVA.clear();
	m_chargedEnergy.clear();
	m_chargedEnergyFromVertexing.clear();
	m_chargedEnergyFromPVA.clear();
	m_chargedMomentum.clear();
	m_chargedMomentumFromVertexing.clear();
	m_chargedMomentumFromPVA.clear();
	m_expectedNeutralEnergy.clear();
	m_expectedNeutralEnergyFromVertexing.clear();
	m_expectedNeutralEnergyFromPVA.clear();
	m_expectedNeutralMomentum.clear();
	m_expectedNeutralMomentumFromVertexing.clear();
	m_expectedNeutralMomentumFromPVA.clear();
	m_expectedChargedEnergy.clear();
	m_expectedChargedEnergyFromVertexing.clear();
	m_expectedChargedEnergyFromPVA.clear();
	m_expectedChargedMomentum.clear();
	m_expectedChargedMomentumFromVertexing.clear();
	m_expectedChargedMomentumFromPVA.clear();
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
	int nTauNeutrino = 0;
	++m_nEvtSum;
	this->Clear();
	streamlog_out(MESSAGE) << "" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////	Processing event 	" << m_nEvt << "	////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;

	IMPL::LCCollectionVec* semiLeptonicVertex(NULL);
	semiLeptonicVertex = new IMPL::LCCollectionVec( LCIO::VERTEX );
//	semiLeptonicVertex->setSubset( true );
	IMPL::LCCollectionVec* semiLeptonicVertexRP(NULL);
	semiLeptonicVertexRP = new IMPL::LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
//	semiLeptonicVertexRP->setSubset( true );
	IMPL::LCCollectionVec* Neutrinos(NULL);
	Neutrinos = new IMPL::LCCollectionVec( LCIO::RECONSTRUCTEDPARTICLE );
	EVENT::LCCollection* JetSLDLink(NULL);
//	JetSLDLink = new EVENT::LCCollection( LCIO::LCRELATION );
	EVENT::LCCollection* SLDJetLink(NULL);
//	SLDJetLink = new EVENT::LCCollection( LCIO::LCRELATION );
	EVENT::LCCollection* mcNurecoNuLink(NULL);
	EVENT::LCCollection* recoNumcNuLink(NULL);
	EVENT::LCCollection* NuSLDLink(NULL);
	EVENT::LCCollection* SLDNuLink(NULL);

	LCRelationNavigator JetSLDRelNav(LCIO::RECONSTRUCTEDPARTICLE , LCIO::VERTEX  );
	LCRelationNavigator SLDJetRelNav(LCIO::VERTEX , LCIO::RECONSTRUCTEDPARTICLE  );
	LCRelationNavigator NeutrinoSLDRelNav(LCIO::RECONSTRUCTEDPARTICLE , LCIO::VERTEX  );
	LCRelationNavigator SLDNeutrinoRelNav(LCIO::VERTEX , LCIO::RECONSTRUCTEDPARTICLE  );
	LCRelationNavigator MCNuRecoNuRelNav(LCIO::MCPARTICLE , LCIO::RECONSTRUCTEDPARTICLE  );
	LCRelationNavigator RecoNuMCNuRelNav(LCIO::RECONSTRUCTEDPARTICLE , LCIO::MCPARTICLE  );

	try
	{
		vtxVector semiLeptonicVertices{};
		pfoVector semiLeptonicVertexRecoParticles{};
		pfoVector neutrinos{};
		pfoVector jetsOfSemiLeptonicDecays{};
		mcpVector mcNeutrinos{};

		vtxVector tempSemiLeptonicVertices{};
		pfoVector tempSemiLeptonicVertexRecoParticles{};
		pfoVector tempNeutrinos{};
		pfoVector tempJetsOfSemiLeptonicDecays{};

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
			if ( ( abs( testLepton->getPDG() ) == 11 || abs( testLepton->getPDG() ) == 13 || abs( testLepton->getPDG() ) == 15 ) )
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
					if ( isBHadronSLDecay )
					{
						++m_nSLDecayOfBHadron;
						m_SLDFlavour.push_back( 5 );
						if ( m_fillRootTree ) h_SLDecayFlavour->Fill( 0.5 );
					}
					else if ( isCHadronSLDecay )
					{
						++m_nSLDecayOfCHadron;
						m_SLDFlavour.push_back( 4 );
						if ( m_fillRootTree ) h_SLDecayFlavour->Fill( 1.5 );
					}
					else if ( isTauLeptonSLDecay )
					{
						++m_nSLDecayOfTauLepton;
						m_SLDFlavour.push_back( 15 );
						if ( m_fillRootTree ) h_SLDecayFlavour->Fill( 2.5 );
					}
					else
					{
						m_SLDFlavour.push_back( 0 );
					}
					if ( downStreamSLDecay || upStreamSLDecay )
					{
						m_SLDType.push_back( 0 );
					}
					else
					{
						m_SLDType.push_back( 1 );
					}
					m_SLDLeptonID.push_back( testLepton->getPDG() );
					if ( abs( testLepton->getPDG() ) == 11 )
					{
						++m_nSLDecayToElectron;
						if ( isBHadronSLDecay && m_fillRootTree ) h_SLDecayModeB->Fill( 0.5 );
						if ( isCHadronSLDecay && m_fillRootTree ) h_SLDecayModeC->Fill( 0.5 );
					}
					else if ( abs( testLepton->getPDG() ) == 13 )
					{
						++m_nSLDecayToMuon;
						if ( isBHadronSLDecay && m_fillRootTree ) h_SLDecayModeB->Fill( 1.5 );
						if ( isCHadronSLDecay && m_fillRootTree ) h_SLDecayModeC->Fill( 1.5 );
					}
					else if ( abs( testLepton->getPDG() ) == 15 )
					{
						++m_nSLDecayToTau;
						if ( isBHadronSLDecay && m_fillRootTree ) h_SLDecayModeB->Fill( 2.5 );
						if ( isCHadronSLDecay && m_fillRootTree ) h_SLDecayModeC->Fill( 2.5 );
					}
					m_SLDecayXi.push_back( testLepton->getParents()[ 0 ]->getVertex()[ 0 ] );
					m_SLDecayYi.push_back( testLepton->getParents()[ 0 ]->getVertex()[ 1 ] );
					m_SLDecayZi.push_back( testLepton->getParents()[ 0 ]->getVertex()[ 2 ] );
					m_SLDecayRi.push_back( sqrt( pow( testLepton->getParents()[ 0 ]->getVertex()[ 0 ] , 2 ) + pow( testLepton->getParents()[ 0 ]->getVertex()[ 1 ] , 2 ) + pow( testLepton->getParents()[ 0 ]->getVertex()[ 2 ] , 2 ) ) );
					m_SLDecayXf.push_back( testLepton->getParents()[ 0 ]->getEndpoint()[ 0 ] );
					m_SLDecayYf.push_back( testLepton->getParents()[ 0 ]->getEndpoint()[ 1 ] );
					m_SLDecayZf.push_back( testLepton->getParents()[ 0 ]->getEndpoint()[ 2 ] );
					m_SLDecayRf.push_back( sqrt( pow( testLepton->getParents()[ 0 ]->getEndpoint()[ 0 ] , 2 ) + pow( testLepton->getParents()[ 0 ]->getEndpoint()[ 1 ] , 2 ) + pow( testLepton->getParents()[ 0 ]->getEndpoint()[ 2 ] , 2 ) ) );
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
					}
					if ( upStreamSLDecay )
					{
						streamlog_out(DEBUG3) << "	There is/are upstream semi-leptonic(s) decay in primary semi-leptonic decay products" << std::endl;
					}
					if ( !downStreamSLDecay && !upStreamSLDecay )
					{
						if ( m_fillRootTree ) h_SLDecayOrder->Fill( 1.5 );
						streamlog_out(DEBUG3) << "	<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
						streamlog_out(DEBUG3) << "	<<<<<<<<<<<<<<<< There are no upstream and downstream semi-leptonic decay >>>>>>>>>>>>>>>>>" << std::endl;
						streamlog_out(DEBUG3) << "	<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
						doSLDCorrection( pLCEvent , testLepton , semiLeptonicVertices , semiLeptonicVertexRecoParticles , jetsOfSemiLeptonicDecays , neutrinos );
						m_parentHadronMass.push_back( ( testLepton->getParents()[ 0 ] )->getMass() );
						m_parentHadronPDG.push_back( ( testLepton->getParents()[ 0 ] )->getPDG() );
						for ( unsigned int i_Btype = 0 ; i_Btype < BHadPDGs.size() ; ++i_Btype )
						{
							if ( abs( ( testLepton->getParents()[ 0 ] )->getPDG() ) == BHadPDGs[ i_Btype ] && m_fillRootTree ) h_BHadronType->Fill( i_Btype );
						}
						m_trueParentHadronFlightDistance.push_back( std::sqrt( std::pow( testLepton->getParents()[ 0 ]->getEndpoint()[ 0 ] - testLepton->getParents()[ 0 ]->getVertex()[ 0 ] , 2 ) + std::pow( testLepton->getParents()[ 0 ]->getEndpoint()[ 1 ] - testLepton->getParents()[ 0 ]->getVertex()[ 1 ] , 2 ) + std::pow( testLepton->getParents()[ 0 ]->getEndpoint()[ 2 ] - testLepton->getParents()[ 0 ]->getVertex()[ 2 ] , 2 ) ) );
						for ( unsigned int i_d = 0 ; i_d < ( testLepton->getParents()[ 0 ] )->getDaughters().size() ; ++i_d )
						{
							MCParticle *mcDaughter = ( testLepton->getParents()[ 0 ] )->getDaughters()[ i_d ];
							int daughterPDG = std::abs( mcDaughter->getPDG() );
							if ( daughterPDG < 11 || daughterPDG > 16 )
							{
								m_daughterHadronMass.push_back( mcDaughter->getMass() );
								m_daughterHadronPDG.push_back( mcDaughter->getPDG() );
								for ( unsigned int i_Ctype = 0 ; i_Ctype < CHadPDGs.size() ; ++i_Ctype )
								{
									if ( abs( mcDaughter->getPDG() ) == CHadPDGs[ i_Ctype ] && m_fillRootTree ) h_CHadronType->Fill( i_Ctype );
								}
							}
							if ( daughterPDG == std::abs( testLepton->getPDG() ) + 1 ) mcNeutrinos.push_back( mcDaughter );
						}
					}
					else
					{
						if ( m_fillRootTree ) h_SLDecayOrder->Fill( 0.5 );
					}
				}
			}
		}
		m_nTauNeutrino = nTauNeutrino;
		m_nSLDecayTotal = m_nSLDecayOfBHadron + m_nSLDecayOfCHadron;

		if ( m_fillRootTree )
		{
			m_pTTree1->Fill();
			m_pTTree2->Fill();
		}
		semiLeptonicVertex->parameters().setValue( "nBHadronSLD_found" , ( int )m_nSLDecayOfBHadron );
		semiLeptonicVertex->parameters().setValue( "nCHadronSLD_found" , ( int )m_nSLDecayOfCHadron );
		semiLeptonicVertex->parameters().setValue( "nTauLeptonSLD_found" , ( int )m_nSLDecayOfTauLepton );
		semiLeptonicVertex->parameters().setValue( "nTotalSLD_found" , ( int )m_nSLDecayTotal );
		for ( unsigned int i_sld = 0 ; i_sld < semiLeptonicVertices.size() ; ++i_sld )
		{
			semiLeptonicVertex->addElement( semiLeptonicVertices[ i_sld ] );
			semiLeptonicVertexRP->addElement( semiLeptonicVertexRecoParticles[ i_sld ] );
			for ( unsigned int i_nu = 0 ; i_nu < 3 ; ++i_nu )
			{
				Neutrinos->addElement( neutrinos[ i_sld * 3 + i_nu ] );
				NeutrinoSLDRelNav.addRelation( neutrinos[ i_sld * 3 + i_nu ] , semiLeptonicVertices[ i_sld ] , 1.0 );
				SLDNeutrinoRelNav.addRelation( semiLeptonicVertices[ i_sld ] , neutrinos[ i_sld * 3 + i_nu ] , 1.0 );
			}
			JetSLDRelNav.addRelation( jetsOfSemiLeptonicDecays[ i_sld ] , semiLeptonicVertices[ i_sld ] , 1.0 );
			SLDJetRelNav.addRelation( semiLeptonicVertices[ i_sld ] , jetsOfSemiLeptonicDecays[ i_sld ] , 1.0 );
		}
		if ( mcNeutrinos.size() == neutrinos.size() / 3 )
		{
			for ( unsigned int i_nu = 0 ; i_nu < mcNeutrinos.size() ; ++i_nu )
			{
				MCNuRecoNuRelNav.addRelation( mcNeutrinos[ i_nu ] , neutrinos[ i_nu * 3 ] , 1.0 );
				RecoNuMCNuRelNav.addRelation( neutrinos[ i_nu * 3 ] , mcNeutrinos[ i_nu ] , 1.0 );
			}
			mcNurecoNuLink = MCNuRecoNuRelNav.createLCCollection();
			recoNumcNuLink = RecoNuMCNuRelNav.createLCCollection();
		}
		JetSLDLink = JetSLDRelNav.createLCCollection();
		SLDJetLink = SLDJetRelNav.createLCCollection();
		NuSLDLink = NeutrinoSLDRelNav.createLCCollection();
		SLDNuLink = SLDNeutrinoRelNav.createLCCollection();
		pLCEvent->addCollection( semiLeptonicVertex , m_SLDVertex );
		pLCEvent->addCollection( semiLeptonicVertexRP , m_SLDVertexRP );
		pLCEvent->addCollection( Neutrinos , m_reconstructedNeutrino );
		pLCEvent->addCollection( JetSLDLink , m_JetSLDLinkName );
		pLCEvent->addCollection( SLDJetLink , m_SLDJetLinkName );
		pLCEvent->addCollection( NuSLDLink , m_NuSLDLinkName );
		pLCEvent->addCollection( SLDNuLink , m_SLDNuLinkName );
		if ( mcNeutrinos.size() == neutrinos.size() / 3 )
		{
			pLCEvent->addCollection( mcNurecoNuLink , m_mcNurecoNuLinkName );
			pLCEvent->addCollection( recoNumcNuLink , m_recoNumcNuLinkName );
		}
//		if ( semiLeptonicVertices.size() != 0 )
//		{
//		}

		streamlog_out(MESSAGE) << "	Found " << semiLeptonicVertices.size() << " semi-leptonic decays with " << semiLeptonicVertexRecoParticles.size() << " associated reconstructed particles in " << jetsOfSemiLeptonicDecays.size() << " jets" << std::endl;
//		semiLeptonicVertices.clear();
//		semiLeptonicVertexRecoParticles.clear();
//		jetsOfSemiLeptonicDecays.clear();

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
			if ( daughter->isOverlay() ) continue;
			if ( daughter->getGeneratorStatus() == 1 )
			{
				if ( abs( daughter->getPDG() ) == 11 || abs( daughter->getPDG() ) == 13 || abs( daughter->getPDG() ) == 15 )
				{
					int leptonCharge = ( int ) daughter->getCharge();
					int expectedNeutrinoPDG = -1 * ( daughter->getPDG() - leptonCharge );
					for ( long unsigned int i_NuCandidate = 0 ; i_NuCandidate < ( parentHadron->getDaughters() ).size() ; ++i_NuCandidate )
					{
						MCParticle *secondDaughter = parentHadron->getDaughters()[ i_NuCandidate ];
						if ( secondDaughter->getPDG() == expectedNeutrinoPDG && secondDaughter->getGeneratorStatus() == 1 ) hasSLDecay = true;
					}
				}
			}
			else if( abs( daughter->getPDG() ) == 15 )
			{
				int leptonCharge = ( int ) daughter->getCharge();
				int expectedNeutrinoPDG = -1 * ( daughter->getPDG() - leptonCharge );
				for ( long unsigned int i_NuCandidate = 0 ; i_NuCandidate < ( parentHadron->getDaughters() ).size() ; ++i_NuCandidate )
				{
					MCParticle *secondDaughter = parentHadron->getDaughters()[ i_NuCandidate ];
					if ( secondDaughter->getPDG() == expectedNeutrinoPDG && secondDaughter->getGeneratorStatus() == 1 ) hasSLDecay = true;
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

void SLDCorrection::doSLDCorrection( EVENT::LCEvent *pLCEvent , MCParticle *SLDLepton , vtxVector& semiLeptonicVertices , pfoVector& semiLeptonicVertexRecoParticles , pfoVector& jetsOfSemiLeptonicDecays , pfoVector& neutrinos )
{
	++n_SLDStatus;
	showTrueParameters( SLDLepton );

	VertexImpl* semiLeptonicVertex = new VertexImpl;
	ReconstructedParticleImpl* semiLeptonicVertexRecoParticle = new ReconstructedParticleImpl;
	ReconstructedParticleImpl* recoNeutrinoZero = new ReconstructedParticleImpl;
	ReconstructedParticleImpl* recoNeutrinoPos = new ReconstructedParticleImpl;
	ReconstructedParticleImpl* recoNeutrinoNeg = new ReconstructedParticleImpl;

	vtxVector SLDVertices{};
	pfoVector SLDVerticesRP{};

	mcpVector trueNeutralDecayProducts{};
	mcpVector trueChargedDecayProducts{};

	pfoVector recoNeutralDecayProducts{};
	pfoVector recoChargedDecayProducts{};

	pfoVector jetVector{};
	pfoVector allPFOsInJet{};
	pfoVector allPFOsInJetFromSLD{};
	pfoVector allPFOsInJetNotFromSLD{};
	pfoVector chargedPFOsInJet{};
	pfoVector chargedPFOsInJetFromSLD{};
	pfoVector chargedPFOsInJetNotFromSLD{};
	pfoVector neutralPFOsInJet{};
	pfoVector neutralPFOsInJetFromSLD{};
	pfoVector neutralPFOsInJetNotFromSLD{};
	pfoVector aloneChargedPFOsInJet{};
	pfoVector aloneChargedPFOsInJetFromSLD{};
	pfoVector aloneChargedPFOsInJetNotFromSLD{};

	pfoVector associatedParticles{};
	pfoVector associatedNeutralParticle{};
	pfoVector associatedChargedParticle{};

	pfoVector decayProducts{};
	pfoVector neutralDecayProducts{};
	pfoVector chargedDecayProducts{};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////									    ////
////			Get Primary Information from event		    ////
////									    ////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

	LCRelationNavigator RecoMCParticleNav( pLCEvent->getCollection( m_RecoMCTruthLinkCollection ) );
	LCRelationNavigator MCParticleRecoNav( pLCEvent->getCollection( m_MCTruthRecoLinkCollection ) );
	LCCollection *primaryVertexCollection = pLCEvent->getCollection( m_inputPrimaryVertex );
	Vertex* primaryVertex = dynamic_cast<Vertex*>( primaryVertexCollection->getElementAt( 0 ) );
	Vertex* startVertex = dynamic_cast<Vertex*>( primaryVertexCollection->getElementAt( 0 ) );
	LCCollection *jetCollection = pLCEvent->getCollection( m_inputJetCollection );
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

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////									    ////
////		Draw semi-leptonic decay with MCParticles		    ////
////									    ////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

	if ( m_displayEvent )
	{
		DDMarlinCED::newEvent( this ); // refresh
		DDMarlinCED::drawDD4hepDetector( this->_theDetector , 0 , std::vector<std::string>{} ); // draw geometry
		DDCEDPickingHandler& pHandler = DDCEDPickingHandler::getInstance();
		pHandler.update(pLCEvent);
		drawMCParticles( SLDLepton->getParents()[ 0 ] );
	}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////									    ////
////	Initialize Input/Output variables for neutrino correction	    ////
////									    ////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

	MCParticle* trueNeutrino = getTrueNeutrino( SLDLepton );
	MCParticle *parentHadron = SLDLepton->getParents()[ 0 ];
	TLorentzVector trueVisible4Mom( 0.0 , 0.0 , 0.0 , 0.0 );
	for ( unsigned int i_d = 0 ; i_d < ( SLDLepton->getParents()[ 0 ] )->getDaughters().size() ; ++i_d )
	{
		if ( ( SLDLepton->getParents()[ 0 ] )->getDaughters()[ i_d ] != trueNeutrino ) trueVisible4Mom += TLorentzVector( ( SLDLepton->getParents()[ 0 ] )->getDaughters()[ i_d ]->getMomentum() , ( SLDLepton->getParents()[ 0 ] )->getDaughters()[ i_d ]->getEnergy() );
	}

	TLorentzVector trueNeutrinoFourMomentum( trueNeutrino->getMomentum() , trueNeutrino->getEnergy() );
	TLorentzVector trueHadronFourMomentum( parentHadron->getMomentum() , parentHadron->getEnergy() );
	TLorentzVector trueLeptonFourMomentum( SLDLepton->getMomentum() , SLDLepton->getEnergy() );
	TLorentzVector trueChargedFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector trueNeutralFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector trueVisibleFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );

	TLorentzVector cheatedPVARecoChargedFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector cheatedPVARecoNeutralFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );

	TLorentzVector recoLeptonFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector recoChargedFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector recoNeutralFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector recoVisibleFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector recoNeutrinoFourMomentumPos( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector recoNeutrinoFourMomentumNeg( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector recoNeutrinoFourMomentumClose( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector recoHadronFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );

	std::vector< float > NeutrinoCovMat( 10 , 0.0 );
	std::vector< float > NeutrinoCovMatPos( 10 , 0.0 );
	std::vector< float > NeutrinoCovMatNeg( 10 , 0.0 );
	std::vector< float > NeutrinoCovMatZero( 10 , 0.0 );
	std::vector< float > CovMatrixPVA( 10, 0.0 );
	std::vector< float > CovMatrixFlightDirection( 6, 0.0 );
	std::vector< float > CovMatrixDetPar( 10 , 0.0 );
	std::vector< float > CovMatrixDetNor( 10 , 0.0 );
	std::vector< float > CovMatrixDetector( 10 , 0.0 );

	TLorentzVector leptonFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector chargedFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector neutralFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector visibleFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );

	std::vector<float> truePrimaryVertex{};
	truePrimaryVertex.push_back( parentHadron->getVertex()[ 0 ] );
	truePrimaryVertex.push_back( parentHadron->getVertex()[ 1 ] );
	truePrimaryVertex.push_back( parentHadron->getVertex()[ 2 ] );

	TVector3 trueFlightDirection( 0.0 , 0.0 , 0.0 );
	TVector3 recoFlightDirection( 0.0 , 0.0 , 0.0 );
	TVector3 flightDirection( 0.0 , 0.0 , 0.0 );
	double parentHadronMass = parentHadron->getMass();// Cheated for the time being
	float helicesDistance = 0.0;
/*
	float expectedNeutralEnergy = 0.0;
	float expectedNeutralEnergyFromVertexing = 0.0;
	float expectedNeutralEnergyFromPVA = 0.0;
	float expectedChargedEnergy = 0.0;
	float expectedChargedEnergyFromVertexing = 0.0;
	float expectedChargedEnergyFromPVA = 0.0;

	float assignedNeutralEnergy = 0.0;
	float assignedChargedEnergy = 0.0;
	float correctAssignedNeutralEnergy = 0.0;
	float correctAssignedChargedEnergy = 0.0;
	float wrongAssignedNeutralEnergy = 0.0;
	float wrongAssignedChargedEnergy = 0.0;

	float neutralEnergy = 0.0;
	float chargedEnergy = 0.0;
	float associatedNeutralEnergy = 0.0;
	float associatedChargedEnergy = 0.0;
	float missedNeutralEnergy = 0.0;
	float missedChargedEnergy = 0.0;
*/
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////									    ////
////	Get True Visible decay products, four-momentum & flight direction   ////
////									    ////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

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
	int nTrueChargedMCPs = getChargedMCPs( SLDLepton , parentHadron , trueChargedDecayProducts );
	streamlog_out(DEBUG2) << "----------------------------------------------------------------------" << std::endl;
	streamlog_out(DEBUG2) << "----------------- " << nTrueChargedMCPs << " Stable Charged MCParticles Found------------------" << std::endl;
	streamlog_out(DEBUG2) << "----------------------------------------------------------------------" << std::endl;
	streamlog_out(DEBUG2) << "" << std::endl;
	for ( unsigned int i_par = 0 ; i_par < trueNeutralDecayProducts.size() ; ++i_par )
	{
		trueNeutralFourMomentum += TLorentzVector( ( trueNeutralDecayProducts[ i_par ] )->getMomentum() , ( trueNeutralDecayProducts[ i_par ] )->getEnergy() );
	}
	for ( unsigned int i_par = 0 ; i_par < trueChargedDecayProducts.size() ; ++i_par )
	{
		trueChargedFourMomentum += TLorentzVector( ( trueChargedDecayProducts[ i_par ] )->getMomentum() , ( trueChargedDecayProducts[ i_par ] )->getEnergy() );
	}
	trueVisibleFourMomentum = trueNeutralFourMomentum + trueChargedFourMomentum + trueLeptonFourMomentum;

	std::vector<double> trueStartVertex{};
	std::vector<double> trueSLDVertex{};
	getTrueFlightDirection( SLDLepton , trueFlightDirection , trueStartVertex , trueSLDVertex );

	m_true_E_vis.push_back( trueVisible4Mom.E() );
	m_true_E_vis_prime.push_back( ( pow( parentHadron->getMass() , 2 ) + pow( trueVisible4Mom.M() , 2 ) ) / ( 2 * parentHadron->getMass() ) );
	TVector3 p_vis = TVector3( trueVisible4Mom.Px() , trueVisible4Mom.Py() , trueVisible4Mom.Pz() );
	m_true_P_vis_par.push_back( p_vis.Dot( trueFlightDirection ) );
	m_true_P_vis_nor.push_back( sqrt( pow( p_vis.Mag() , 2 ) - pow( p_vis.Dot( trueFlightDirection ) , 2 ) ) );
	m_true_P_vis_par_prime.push_back( sqrt( pow( ( pow( parentHadron->getMass() , 2 ) - pow( trueVisible4Mom.M() , 2 ) ) / ( 2 * parentHadron->getMass() ) , 2 ) - pow( sqrt( pow( p_vis.Mag() , 2 ) - pow( p_vis.Dot( trueFlightDirection ) , 2 ) ) , 2 ) ) );
	streamlog_out(DEBUG8) << "		 Visible Energy:		" << trueVisible4Mom.E() << std::endl;
	streamlog_out(DEBUG8) << "		 Visible Energy prime:		" << ( pow( parentHadron->getMass() , 2 ) + pow( trueVisible4Mom.M() , 2 ) ) / ( 2 * parentHadron->getMass() ) << std::endl;
	streamlog_out(DEBUG8) << "		 Visible Momentum:		" << p_vis.Px() << " , " << p_vis.Py() << " , " << p_vis.Pz() << std::endl;
	streamlog_out(DEBUG8) << "		|Visible Momentum (par)|:	" << p_vis.Dot( trueFlightDirection ) << std::endl;
	streamlog_out(DEBUG8) << "		|Visible Momentum (nor)|:	" << sqrt( pow( p_vis.Mag() , 2 ) - pow( p_vis.Dot( trueFlightDirection ) , 2 ) ) << std::endl;
	streamlog_out(DEBUG8) << "		|Visible Momentum (par-prime)|:	" << sqrt( pow( ( pow( parentHadron->getMass() , 2 ) - pow( trueVisible4Mom.M() , 2 ) ) / ( 2 * parentHadron->getMass() ) , 2 ) - pow( sqrt( pow( p_vis.Mag() , 2 ) - pow( p_vis.Dot( trueFlightDirection ) , 2 ) ) , 2 ) ) << std::endl;


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////									    ////
////	Identify Reconstrcuted Lepton and the jet/vertex assigned to SLD    ////
////									    ////
////////////////////////////////////////////////////////////////////////////////
////			If there is no Reconstrcuted Lepton		    ////
////					OR				    ////
////			Reconstrcuted Lepton is not in a jet		    ////
////					OR				    ////
////		     Reconstrcuted Lepton is in Primary Vertex		    ////
////									    ////
////			     Neutrino Correction FAILS			    ////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

	float weightPFOtoMCP = 0.0;
	float weightMCPtoPFO = 0.0;
	ReconstructedParticle* linkedRecoLepton = getLinkedPFO( SLDLepton , RecoMCParticleNav , MCParticleRecoNav , true , false , weightPFOtoMCP , weightMCPtoPFO );
	m_weightPFOtoMCP_Lepton.push_back( weightPFOtoMCP );
	m_weightMCPtoPFO_Lepton.push_back( weightMCPtoPFO );
	if ( linkedRecoLepton == NULL )
	{
		m_SLDStatus.push_back( 1 );
		streamlog_out(WARNING) << "	||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << std::endl;
		streamlog_out(WARNING) << "	||||||||||||||||| Reconstructed Lepton is not found ||||||||||||||||||" << std::endl;
		streamlog_out(WARNING) << "	||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << std::endl;
		if ( m_fillRootTree ) h_SLDStatus->Fill( 0.5 );
		return;
	}
	recoLeptonFourMomentum = TLorentzVector( linkedRecoLepton->getMomentum() , linkedRecoLepton->getEnergy() );

	bool recoLeptonIsInJet = false;
	ReconstructedParticle *assignedJet = NULL;
	assignedJet = getJetAssignedToParticle( linkedRecoLepton , jetVector , recoLeptonIsInJet );
	if ( !recoLeptonIsInJet )
	{
		m_SLDStatus.push_back( 2 );
		streamlog_out(WARNING) << "	||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << std::endl;
		streamlog_out(WARNING) << "	||||||||||| Reconstructed Lepton doesn't belong to any jet |||||||||||" << std::endl;
		streamlog_out(WARNING) << "	||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << std::endl;
		if ( m_fillRootTree ) h_SLDStatus->Fill( 1.5 );
		return;
	}

	if ( isParticleInVertex( linkedRecoLepton , primaryVertex ) )
	{
		m_SLDStatus.push_back( 3 );
		streamlog_out(WARNING) << "	||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << std::endl;
		streamlog_out(WARNING) << "	||||||||||||| Reconstructed Lepton is in primary vertex ||||||||||||||" << std::endl;
		streamlog_out(WARNING) << "	||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << std::endl;
		if ( m_fillRootTree ) h_SLDStatus->Fill( 2.5 );
		return;
	}
	investigateJetEnergyContent( assignedJet );
	jetsOfSemiLeptonicDecays.push_back( assignedJet );

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////									    ////
////	   Investigate Alone charged Particles in the assigned jet	    ////
////									    ////
////		     FromSLD and NotFromSLD is Cheated			    ////
////									    ////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

	aloneChargedPFOsInJet = getParticlesWithAloneTracks( linkedRecoLepton , assignedJet , primaryVertex , buildUpVertexVector );
	for ( unsigned int i_par = 0 ; i_par < aloneChargedPFOsInJet.size() ; ++i_par )
	{
		ReconstructedParticle* testParticle = aloneChargedPFOsInJet[ i_par ];
		weightPFOtoMCP = 0.0;
		weightMCPtoPFO = 0.0;
		MCParticle* linkedChargedMCP = getLinkedMCP( testParticle , RecoMCParticleNav , MCParticleRecoNav , true , false , weightPFOtoMCP , weightMCPtoPFO );
		if ( linkedChargedMCP != NULL )
		{
			bool MCPisFromSLD = false;
			isMCParticleFromSLD( parentHadron , linkedChargedMCP , MCPisFromSLD );
			if ( MCPisFromSLD )
			{
				aloneChargedPFOsInJetFromSLD.push_back( testParticle );
			}
			else
			{
				aloneChargedPFOsInJetNotFromSLD.push_back( testParticle );
			}
		}
	}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////									    ////
////	   Reconstruction of Flight Direction of Parent Hadron		    ////
////									    ////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

	TVector3 leptonDirection = TVector3( SLDLepton->getMomentum() ); leptonDirection.SetMag( 1.0 );
	TVector3 jetAxis = TVector3( assignedJet->getMomentum() ); jetAxis.SetMag( 1.0 );

	vtxVector verticesInJet = getVerticesInJet( assignedJet , buildUpVertexVector );
	m_nRecoVerticesInJet.push_back( verticesInJet.size() );

	int vertexingScenario = m_vertexingScenario;
	pfoVector sortedChargedPFOs;
	if ( aloneChargedPFOsInJet.size() != 0 ) sortParticles( sortedChargedPFOs , aloneChargedPFOsInJet , jetAxis );
	if ( sortedChargedPFOs.size() == 0 && m_vertexingScenario == 4 ) vertexingScenario = 1;
	m_flightDirectionStatus.push_back( vertexingScenario );
	double hadronFlightLength;
	TVector3 daughterHadronFlightDirection;
	double daughterHadronFlightDistance = 0.0;
	std::vector<float> sldVertexPosition{};
	int SLDStatus = getRecoFlightDirection( linkedRecoLepton , recoFlightDirection , hadronFlightLength , primaryVertex , startVertex , SLDVertices , SLDVerticesRP , assignedJet , verticesInJet , sortedChargedPFOs , helicesDistance , vertexingScenario , daughterHadronFlightDirection , daughterHadronFlightDistance , sldVertexPosition );
	m_daughterHadronFlightDistance.push_back( daughterHadronFlightDistance );
	m_recoParentHadronFlightDistance.push_back( hadronFlightLength );
	m_distRecoLeptonToDownStreamVertex.push_back( helicesDistance );
	m_FlightDirectionErrorCosAlpha.push_back( trueFlightDirection.Dot( recoFlightDirection ) );
	m_FlightDirectionErrorSinAlpha.push_back( std::sin( acos( trueFlightDirection.Dot( recoFlightDirection ) ) ) );
	m_FlightDirectionErrorAlpha.push_back( acos( trueFlightDirection.Dot( recoFlightDirection ) ) * 180.0 / 3.14159265 );
	m_FlightDirectionErrorDeltaTheta.push_back( recoFlightDirection.Theta() - trueFlightDirection.Theta() );
	m_FlightDirectionErrorDeltaPhi.push_back( recoFlightDirection.Phi() - trueFlightDirection.Phi() );

	m_SLDStatus.push_back( SLDStatus );
	if ( m_fillRootTree ) h_SLDStatus->Fill( SLDStatus - 0.5 );

	if ( m_cheatFlightDirection )
	{
		flightDirection = trueFlightDirection;
	}
	else
	{
		flightDirection = recoFlightDirection;
	}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////									    ////
////		Particle Association to the semi-leptonic decay		    ////
////									    ////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

	for ( unsigned int i_par = 0 ; i_par < trueNeutralDecayProducts.size() ; ++i_par )
	{
		ReconstructedParticle *linkedPFO = getLinkedPFO( trueNeutralDecayProducts[ i_par ] , RecoMCParticleNav , MCParticleRecoNav , false , true , weightPFOtoMCP , weightMCPtoPFO );
		m_weightPFOtoMCP_Neutral.push_back( weightPFOtoMCP );
		m_weightMCPtoPFO_Neutral.push_back( weightMCPtoPFO );
		if ( linkedPFO != NULL ) recoNeutralDecayProducts.push_back( linkedPFO );
	}
	for ( unsigned int i_par = 0 ; i_par < trueChargedDecayProducts.size() ; ++i_par )
	{
		ReconstructedParticle *linkedPFO = getLinkedPFO( trueChargedDecayProducts[ i_par ] , RecoMCParticleNav , MCParticleRecoNav , true , false , weightPFOtoMCP , weightMCPtoPFO );
		m_weightPFOtoMCP_Charged.push_back( weightPFOtoMCP );
		m_weightMCPtoPFO_Charged.push_back( weightMCPtoPFO );
		if ( linkedPFO != NULL ) recoChargedDecayProducts.push_back( linkedPFO );
	}

	float chargedCosAcceptanceAngle = 0.0;
	if ( SLDStatus == 4 )
	{
		chargedCosAcceptanceAngle = m_chargedCosAcceptanceAngleSLD4;
	}
	else if ( SLDStatus == 5 )
	{
		chargedCosAcceptanceAngle = m_chargedCosAcceptanceAngleSLD5;
	}


	for ( unsigned int i_par = 0 ; i_par < assignedJet->getParticles().size() ; ++i_par )
	{
		ReconstructedParticle* jetParticle = assignedJet->getParticles()[ i_par ];
		TVector3 jetParticleMomentum = TVector3( jetParticle->getMomentum() );
		jetParticleMomentum.SetMag( 1.0 );
		weightPFOtoMCP = 0.0;
		weightMCPtoPFO = 0.0;
		MCParticle* linkedNeutralMCP = getLinkedMCP( jetParticle , RecoMCParticleNav , MCParticleRecoNav , false , true , weightPFOtoMCP , weightMCPtoPFO );
		MCParticle* linkedChargedMCP = getLinkedMCP( jetParticle , RecoMCParticleNav , MCParticleRecoNav , true , false , weightPFOtoMCP , weightMCPtoPFO );
		bool MCPisFromSLD = false;
		if ( jetParticle->getTracks().size() == 0 && jetParticleMomentum.Dot( flightDirection ) >= m_neutralCosAcceptanceAngle )
		{
			streamlog_out(DEBUG2) << "----------------------------------------------------------------------" << std::endl;
			streamlog_out(DEBUG2) << "-------- Added One Neutral PFO to SLDecay products candidates --------" << std::endl;
			streamlog_out(DEBUG2) << "----------------------------------------------------------------------" << std::endl;
			streamlog_out(DEBUG2) << *jetParticle << std::endl;
			allPFOsInJet.push_back( jetParticle );
			neutralPFOsInJet.push_back( jetParticle );
			if ( linkedNeutralMCP != NULL )
			{
				isMCParticleFromSLD( parentHadron , linkedNeutralMCP , MCPisFromSLD );
				if ( MCPisFromSLD )
				{
					neutralPFOsInJetFromSLD.push_back( jetParticle );
					allPFOsInJetFromSLD.push_back( jetParticle );
				}
				else
				{
					neutralPFOsInJetNotFromSLD.push_back( jetParticle );
					allPFOsInJetNotFromSLD.push_back( jetParticle );
				}
			}
		}
		else if ( jetParticle->getTracks().size() != 0 && jetParticleMomentum.Dot( flightDirection ) >= chargedCosAcceptanceAngle )
		{
			streamlog_out(DEBUG2) << "----------------------------------------------------------------------" << std::endl;
			streamlog_out(DEBUG2) << "-------- Added One Charged PFO to SLDecay products candidates --------" << std::endl;
			streamlog_out(DEBUG2) << "----------------------------------------------------------------------" << std::endl;
			streamlog_out(DEBUG2) << *jetParticle << std::endl;
			streamlog_out(DEBUG2) << *jetParticle << std::endl;
			allPFOsInJet.push_back( jetParticle );
			chargedPFOsInJet.push_back( jetParticle );
			if ( linkedChargedMCP != NULL )
			{
				isMCParticleFromSLD( parentHadron , linkedChargedMCP , MCPisFromSLD );
				if ( MCPisFromSLD )
				{
					chargedPFOsInJetFromSLD.push_back( jetParticle );
					allPFOsInJetFromSLD.push_back( jetParticle );
				}
				else
				{
					chargedPFOsInJetNotFromSLD.push_back( jetParticle );
					allPFOsInJetNotFromSLD.push_back( jetParticle );
				}
			}
		}
	}


	evaluatePFOsAngle( aloneChargedPFOsInJetFromSLD , aloneChargedPFOsInJetNotFromSLD , chargedPFOsInJetFromSLD , chargedPFOsInJetNotFromSLD , neutralPFOsInJetFromSLD , neutralPFOsInJetNotFromSLD , leptonDirection , jetAxis , recoFlightDirection , SLDStatus );

//	associatedParticles.push_back( linkedRecoLepton );
//	double InvMassCut = 2.1;
	double InvMassCutCharged = 0.0;
	double InvMassCutNeutral = 0.0;
	if ( SLDStatus == 4 )
	{
		InvMassCutCharged = m_BSLDChargedSLD4InvMassCut;
		InvMassCutNeutral = m_BSLDNeutralSLD4InvMassCut;
		ReconstructedParticle* sldVertexRP = SLDVerticesRP[ 0 ];
		for ( unsigned int i_par = 0 ; i_par < sldVertexRP->getParticles().size() ; ++i_par )
		{
			ReconstructedParticle* chargedDecayProduct = sldVertexRP->getParticles()[ i_par ];
			if ( chargedDecayProduct != linkedRecoLepton )
			{
				associatedChargedParticle.push_back( chargedDecayProduct );
				associatedParticles.push_back( chargedDecayProduct );
//				if ( m_displayEvent ) drawReconstructedParticle( chargedDecayProduct , primaryVertex , 0x0075df , 0x000000 );
			}
		}
	}
	else if ( SLDStatus == 5 )
	{
		InvMassCutCharged = m_BSLDChargedSLD5InvMassCut;
		InvMassCutNeutral = m_BSLDNeutralSLD5InvMassCut;
		ReconstructedParticle* sldVertexRP = SLDVerticesRP[ 0 ];
		for ( unsigned int i_par = 0 ; i_par < sldVertexRP->getParticles().size() ; ++i_par )
		{
			ReconstructedParticle* chargedDecayProduct = sldVertexRP->getParticles()[ i_par ];
			associatedChargedParticle.push_back( chargedDecayProduct );
			associatedParticles.push_back( chargedDecayProduct );
//			if ( m_displayEvent ) drawReconstructedParticle( chargedDecayProduct , primaryVertex , 0x0075df , 0x2e8e04 );
		}
	}

	vtxVector availableVerticesInJet{};
	for ( unsigned int i_vtx = 0 ; i_vtx < verticesInJet.size() ; ++i_vtx )
	{
		if ( verticesInJet[ i_vtx ] != SLDVertices[ 0 ] ) availableVerticesInJet.push_back( verticesInJet[ i_vtx ] );
	}

	if ( m_cheatPVAneutral )
	{
		for ( unsigned int i_par = 0 ; i_par < recoNeutralDecayProducts.size() ; ++i_par )
		{
			decayProducts.push_back( recoNeutralDecayProducts[ i_par ] );
//			if ( m_displayEvent ) drawReconstructedParticle( recoNeutralDecayProducts[ i_par ] , primaryVertex , 0x0075df , 0x2e8e04 );
		}
		if ( m_cheatPVAcharged )
		{
			for ( unsigned int i_par = 0 ; i_par < recoChargedDecayProducts.size() ; ++i_par )
			{
				decayProducts.push_back( recoChargedDecayProducts[ i_par ] );
//				if ( m_displayEvent ) drawReconstructedParticle( recoChargedDecayProducts[ i_par ] , primaryVertex , 0x0075df , 0x000000 );
			}
		}
		else
		{
			for ( unsigned int i_par = 0 ; i_par < associatedParticles.size() ; ++i_par )
			{
				if ( associatedParticles[ i_par ]->getTracks().size() != 0 && associatedParticles[ i_par ] != linkedRecoLepton )
				{
					decayProducts.push_back( associatedParticles[ i_par ] );
				}
			}
			assignVerticesToSemiLeptonicDecay( decayProducts , availableVerticesInJet , InvMassCutCharged , daughterHadronFlightDirection , startVertex );
			assignParticlesToSemiLeptonicDecay( decayProducts , chargedPFOsInJet , InvMassCutCharged , daughterHadronFlightDirection );
		}
	}
	else
	{
		if ( m_cheatPVAcharged )
		{
			for ( unsigned int i_par = 0 ; i_par < recoChargedDecayProducts.size() ; ++i_par )
			{
				decayProducts.push_back( recoChargedDecayProducts[ i_par ] );
//				if ( m_displayEvent ) drawReconstructedParticle( recoChargedDecayProducts[ i_par ] , primaryVertex , 0x0075df , 0x000000 );
			}
		}
		else
		{
			for ( unsigned int i_par = 0 ; i_par < associatedParticles.size() ; ++i_par )
			{
				if ( associatedParticles[ i_par ]->getTracks().size() != 0 && associatedParticles[ i_par ] != linkedRecoLepton )
				{
					decayProducts.push_back( associatedParticles[ i_par ] );
				}
			}
			assignVerticesToSemiLeptonicDecay( decayProducts , availableVerticesInJet , InvMassCutCharged , daughterHadronFlightDirection , startVertex );
			assignParticlesToSemiLeptonicDecay( decayProducts , chargedPFOsInJet , InvMassCutCharged , daughterHadronFlightDirection );
		}
		assignParticlesToSemiLeptonicDecay( decayProducts , neutralPFOsInJet , InvMassCutNeutral , daughterHadronFlightDirection );
	}
	for ( unsigned int i_par = 0 ; i_par < decayProducts.size() ; ++i_par )
	{
		if ( m_displayEvent ) drawReconstructedParticle( decayProducts[ i_par ] , primaryVertex , 0x0075df , 0x2e8e04 );
		if ( decayProducts[ i_par ]->getTracks().size() == 0 )
		{
			neutralDecayProducts.push_back( decayProducts[ i_par ] );
		}
		else
		{
			chargedDecayProducts.push_back( decayProducts[ i_par ] );
		}
	}
	if ( m_displayEvent ) drawReconstructedParticle( linkedRecoLepton , primaryVertex , 0xfc0000 , 0xfc0000 );

	m_visibleChargedInvMassCut.push_back( InvMassCutCharged );
	m_visibleNeutralInvMassCut.push_back( InvMassCutNeutral );


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////									    ////
////		Form visible four-momentum for neutrino correction	    ////
////									    ////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

	for ( unsigned int i_par = 0 ; i_par < chargedDecayProducts.size() ; ++i_par )
	{
		if ( chargedDecayProducts[ i_par ]->getTracks().size() != 0 && chargedDecayProducts[ i_par ] != linkedRecoLepton ) recoChargedFourMomentum += TLorentzVector( chargedDecayProducts[ i_par ]->getMomentum() , chargedDecayProducts[ i_par ]->getEnergy() );
	}

	for ( unsigned int i_par = 0 ; i_par < neutralDecayProducts.size() ; ++i_par )
	{
		if ( neutralDecayProducts[ i_par ]->getTracks().size() == 0 ) recoNeutralFourMomentum += TLorentzVector( neutralDecayProducts[ i_par ]->getMomentum() , neutralDecayProducts[ i_par ]->getEnergy() );
	}

	for ( unsigned int i_par = 0 ; i_par < recoChargedDecayProducts.size() ; ++i_par )
	{
		if ( recoChargedDecayProducts[ i_par ]->getTracks().size() != 0 && recoChargedDecayProducts[ i_par ] != linkedRecoLepton ) cheatedPVARecoChargedFourMomentum += TLorentzVector( recoChargedDecayProducts[ i_par ]->getMomentum() , recoChargedDecayProducts[ i_par ]->getEnergy() );
	}

	for ( unsigned int i_par = 0 ; i_par < recoNeutralDecayProducts.size() ; ++i_par )
	{
		if ( recoNeutralDecayProducts[ i_par ]->getTracks().size() == 0 ) cheatedPVARecoNeutralFourMomentum += TLorentzVector( recoNeutralDecayProducts[ i_par ]->getMomentum() , recoNeutralDecayProducts[ i_par ]->getEnergy() );
	}

	if ( m_cheatLepton4momentum )
	{
		leptonFourMomentum = trueLeptonFourMomentum;
	}
	else
	{
		leptonFourMomentum = recoLeptonFourMomentum;
	}
	if ( m_cheatNeutral4momentum )
	{
		neutralFourMomentum = trueNeutralFourMomentum;
	}
	else
	{
		neutralFourMomentum = recoNeutralFourMomentum;
	}
	if ( m_cheatCharged4momentum )
	{
		chargedFourMomentum = trueChargedFourMomentum;
	}
	else
	{
		chargedFourMomentum = recoChargedFourMomentum;
	}
	recoVisibleFourMomentum = recoLeptonFourMomentum + recoNeutralFourMomentum + recoChargedFourMomentum;
	visibleFourMomentum = leptonFourMomentum + neutralFourMomentum + chargedFourMomentum;

	streamlog_out(DEBUG8) << "			     (  PDG	, Mass		, Px		, Py		, Pz		, E		, Charge	)" << std::endl;
	streamlog_out(DEBUG8) << "		Neutrino" << std::endl;
	streamlog_out(DEBUG8) << "			True:(	" << trueNeutrino->getPDG() << "	, " << trueNeutrino->getMass() << "	, " << trueNeutrino->getMomentum()[ 0 ] << "	, " << trueNeutrino->getMomentum()[ 1 ] << "	, " << trueNeutrino->getMomentum()[ 2 ] << "	, " << trueNeutrino->getEnergy() << "	, " << trueNeutrino->getCharge() << "	)" << std::endl;
	streamlog_out(DEBUG8) << "" << std::endl;
	streamlog_out(DEBUG8) << "		Hadron" << std::endl;
	streamlog_out(DEBUG8) << "			True:(	" << parentHadron->getPDG() << "	, " << parentHadron->getMass() << "	, " << parentHadron->getMomentum()[ 0 ] << "	, " << parentHadron->getMomentum()[ 1 ] << "	, " << parentHadron->getMomentum()[ 2 ] << "	, " << parentHadron->getEnergy() << "	, " << parentHadron->getCharge() << "	)" << std::endl;
	streamlog_out(DEBUG8) << "" << std::endl;
	streamlog_out(DEBUG8) << "		Lepton" << std::endl;
	streamlog_out(DEBUG8) << "			True:(	" << SLDLepton->getPDG() << "	, " << SLDLepton->getMass() << "	, " << SLDLepton->getMomentum()[ 0 ] << "	, " << SLDLepton->getMomentum()[ 1 ] << "	, " << SLDLepton->getMomentum()[ 2 ] << "	, " << SLDLepton->getEnergy() << "	, " << SLDLepton->getCharge() << "	)" << std::endl;
	streamlog_out(DEBUG8) << "			Reco:(	" << linkedRecoLepton->getType() << "	, " << linkedRecoLepton->getMass() << "	, " << linkedRecoLepton->getMomentum()[ 0 ] << "	, " << linkedRecoLepton->getMomentum()[ 1 ] << "	, " << linkedRecoLepton->getMomentum()[ 2 ] << "	, " << linkedRecoLepton->getEnergy() << "	, " << linkedRecoLepton->getCharge() << "	)" << std::endl;
	streamlog_out(DEBUG8) << "			Used:(	" << "xxx" << "	, " << leptonFourMomentum.M() << "	, " << leptonFourMomentum.Px() << "	, " << leptonFourMomentum.Py() << "	, " << leptonFourMomentum.Pz() << "	, " << leptonFourMomentum.E() << "	, " << (int)SLDLepton->getCharge() << "		)" << std::endl;
	streamlog_out(DEBUG8) << "" << std::endl;
	streamlog_out(DEBUG8) << "		Charged" << std::endl;
	streamlog_out(DEBUG8) << "			True:(	" << "qqq" << "	, " << trueChargedFourMomentum.M() << "	, " << trueChargedFourMomentum.Px() << "	, " << trueChargedFourMomentum.Py() << "	, " << trueChargedFourMomentum.Pz() << "	, " << trueChargedFourMomentum.E() << "	, " << "qqq" << "		)" << std::endl;
	streamlog_out(DEBUG8) << "			Reco:(	" << "qqq" << "	, " << recoChargedFourMomentum.M() << "	, " << recoChargedFourMomentum.Px() << "	, " << recoChargedFourMomentum.Py() << "	, " << recoChargedFourMomentum.Pz() << "	, " << recoChargedFourMomentum.E() << "	, " << "qqq" << "		)" << std::endl;
	streamlog_out(DEBUG8) << "			Used:(	" << "qqq" << "	, " << chargedFourMomentum.M() << "	, " << chargedFourMomentum.Px() << "	, " << chargedFourMomentum.Py() << "	, " << chargedFourMomentum.Pz() << "	, " << chargedFourMomentum.E() << "	, " << "nnn" << "		)" << std::endl;
	streamlog_out(DEBUG8) << "" << std::endl;
	streamlog_out(DEBUG8) << "		Neutral" << std::endl;
	streamlog_out(DEBUG8) << "			True:(	" << "nnn" << "	, " << trueNeutralFourMomentum.M() << "	, " << trueNeutralFourMomentum.Px() << "	, " << trueNeutralFourMomentum.Py() << "	, " << trueNeutralFourMomentum.Pz() << "	, " << trueNeutralFourMomentum.E() << "	, " << "nnn" << "		)" << std::endl;
	streamlog_out(DEBUG8) << "			Reco:(	" << "nnn" << "	, " << recoNeutralFourMomentum.M() << "	, " << recoNeutralFourMomentum.Px() << "	, " << recoNeutralFourMomentum.Py() << "	, " << recoNeutralFourMomentum.Pz() << "	, " << recoNeutralFourMomentum.E() << "	, " << "nnn" << "		)" << std::endl;
	streamlog_out(DEBUG8) << "			Used:(	" << "nnn" << "	, " << neutralFourMomentum.M() << "	, " << neutralFourMomentum.Px() << "	, " << neutralFourMomentum.Py() << "	, " << neutralFourMomentum.Pz() << "	, " << neutralFourMomentum.E() << "	, " << "nnn" << "		)" << std::endl;
	streamlog_out(DEBUG8) << "" << std::endl;
	streamlog_out(DEBUG8) << "		Visible" << std::endl;
	streamlog_out(DEBUG8) << "			True:(	" << "nnn" << "	, " << trueVisibleFourMomentum.M() << "	, " << trueVisibleFourMomentum.Px() << "	, " << trueVisibleFourMomentum.Py() << "	, " << trueVisibleFourMomentum.Pz() << "	, " << trueVisibleFourMomentum.E() << "	, " << "nnn" << "		)" << std::endl;
	streamlog_out(DEBUG8) << "			Reco:(	" << "nnn" << "	, " << recoVisibleFourMomentum.M() << "	, " << recoVisibleFourMomentum.Px() << "	, " << recoVisibleFourMomentum.Py() << "	, " << recoVisibleFourMomentum.Pz() << "	, " << recoVisibleFourMomentum.E() << "	, " << "nnn" << "		)" << std::endl;
	streamlog_out(DEBUG8) << "			Used:(	" << "nnn" << "	, " << visibleFourMomentum.M() << "	, " << visibleFourMomentum.Px() << "	, " << visibleFourMomentum.Py() << "	, " << visibleFourMomentum.Pz() << "	, " << visibleFourMomentum.E() << "	, " << "nnn" << "		)" << std::endl;
	streamlog_out(DEBUG8) << "" << std::endl;
	streamlog_out(DEBUG8) << "			     	(  X		, Y		, Z	)" << std::endl;
	streamlog_out(DEBUG8) << "		Flight Direction" << std::endl;
	streamlog_out(DEBUG8) << "			True:	(  " << trueFlightDirection.X() << "		, " << trueFlightDirection.Y() << "		, " << trueFlightDirection.Z() << "	)" << std::endl;
	streamlog_out(DEBUG8) << "			Reco:	(  " << recoFlightDirection.X() << "		, " << recoFlightDirection.Y() << "		, " << recoFlightDirection.Z() << "	)" << std::endl;
	streamlog_out(DEBUG8) << "			Used:	(  " << flightDirection.X() << "		, " << flightDirection.Y() << "		, " << flightDirection.Z() << "	)" << std::endl;

//	recoNeutrinoFourMomentumPos = getNeutrinoFourMomentum( flightDirection , visibleFourMomentum , parentHadronMass , +1.0 );
//	recoNeutrinoFourMomentumNeg = getNeutrinoFourMomentum( flightDirection , visibleFourMomentum , parentHadronMass , -1.0 );

	recoNeutrinoFourMomentumPos = getNeutrinoFourMomentumModified( flightDirection , visibleFourMomentum , parentHadronMass , +1.0 );
	recoNeutrinoFourMomentumNeg = getNeutrinoFourMomentumModified( flightDirection , visibleFourMomentum , parentHadronMass , -1.0 );

//	recoNeutrinoFourMomentumPos = getNeutrinoFourMomentumStandardMethod( flightDirection , visibleFourMomentum , parentHadronMass , +1.0 );
//	recoNeutrinoFourMomentumNeg = getNeutrinoFourMomentumStandardMethod( flightDirection , visibleFourMomentum , parentHadronMass , -1.0 );

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////									    ////
////		Estimate Errors due to Particle to Vertex Association	    ////
////									    ////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

	getCovMatPVA( decayProducts , associatedParticles , SLDStatus , CovMatrixPVA );
	float sigmaTheta = 0.0;
	float sigmaPhi = 0.0;
	if ( m_cheatFlightDirection )
	{
		sigmaTheta = 0.0;
		sigmaPhi = 0.0;
	}
	else
	{
		if ( SLDStatus == 4 )
		{
			sigmaTheta = 0.015;
			sigmaPhi = 0.018;
		}
		else if ( SLDStatus == 5 )
		{
			sigmaTheta = 0.030;
			sigmaPhi = 0.030;
		}
	}
	getCovMatFlightDirection( flightDirection , sigmaTheta , sigmaPhi , CovMatrixFlightDirection );
	getCovMatDetFlightDirection( decayProducts , flightDirection , CovMatrixFlightDirection , CovMatrixDetector , linkedRecoLepton , CovMatrixDetPar , CovMatrixDetNor );
	getNeutrinoCovMat( recoNeutrinoFourMomentumPos , visibleFourMomentum , flightDirection , parentHadronMass , CovMatrixPVA , CovMatrixDetector , CovMatrixDetPar , CovMatrixDetNor , NeutrinoCovMatPos );
	getNeutrinoCovMat( recoNeutrinoFourMomentumNeg , visibleFourMomentum , flightDirection , parentHadronMass , CovMatrixPVA , CovMatrixDetector , CovMatrixDetPar , CovMatrixDetNor , NeutrinoCovMatNeg );
	streamlog_out(DEBUG9) << "	CovMatNeutrino(+) :	" << NeutrinoCovMatPos[ 0 ] << std::endl;
	streamlog_out(DEBUG9) << "				" << NeutrinoCovMatPos[ 1 ] << "	,	" << NeutrinoCovMatPos[ 2 ] << std::endl;
	streamlog_out(DEBUG9) << "				" << NeutrinoCovMatPos[ 3 ] << "	,	" << NeutrinoCovMatPos[ 4 ] << "	,	" << NeutrinoCovMatPos[ 5 ] << std::endl;
	streamlog_out(DEBUG9) << "				" << NeutrinoCovMatPos[ 6 ] << "	,	" << NeutrinoCovMatPos[ 7 ] << "	,	" << NeutrinoCovMatPos[ 8 ] << "	,	" << NeutrinoCovMatPos[ 9 ] << std::endl;
	streamlog_out(DEBUG9) << "" << std::endl;
	streamlog_out(DEBUG9) << "	CovMatNeutrino(-) :	" << NeutrinoCovMatNeg[ 0 ] << std::endl;
	streamlog_out(DEBUG9) << "				" << NeutrinoCovMatNeg[ 1 ] << "	,	" << NeutrinoCovMatNeg[ 2 ] << std::endl;
	streamlog_out(DEBUG9) << "				" << NeutrinoCovMatNeg[ 3 ] << "	,	" << NeutrinoCovMatNeg[ 4 ] << "	,	" << NeutrinoCovMatNeg[ 5 ] << std::endl;
	streamlog_out(DEBUG9) << "				" << NeutrinoCovMatNeg[ 6 ] << "	,	" << NeutrinoCovMatNeg[ 7 ] << "	,	" << NeutrinoCovMatNeg[ 8 ] << "	,	" << NeutrinoCovMatNeg[ 9 ] << std::endl;
	streamlog_out(DEBUG9) << "" << std::endl;

	TLorentzVector expectedNeutralFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector expectedNeutralFourMomentumFromVertexing( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector expectedNeutralFourMomentumFromPVA( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector expectedChargedFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector expectedChargedFourMomentumFromVertexing( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector expectedChargedFourMomentumFromPVA( 0.0 , 0.0 , 0.0 , 0.0 );
	for ( unsigned int i_par = 0 ; i_par < recoNeutralDecayProducts.size() ; ++i_par )
	{
		bool perfectAssociation = false;
		TLorentzVector particleFourMomentum( ( recoNeutralDecayProducts[ i_par ] )->getMomentum() , ( recoNeutralDecayProducts[ i_par ] )->getEnergy() );
		expectedNeutralFourMomentum += particleFourMomentum;
		for ( unsigned int i_p = 0 ; i_p < associatedParticles.size() ; ++i_p )
		{
			if ( recoNeutralDecayProducts[ i_par ] == associatedParticles[ i_p ] ) perfectAssociation = true;
		}
		if ( perfectAssociation )
		{
			expectedNeutralFourMomentumFromVertexing += particleFourMomentum;
		}
		else
		{
			expectedNeutralFourMomentumFromPVA += particleFourMomentum;
		}
	}
	for ( unsigned int i_par = 0 ; i_par < recoChargedDecayProducts.size() ; ++i_par )
	{
		bool perfectAssociation = false;
		TLorentzVector particleFourMomentum( ( recoChargedDecayProducts[ i_par ] )->getMomentum() , ( recoChargedDecayProducts[ i_par ] )->getEnergy() );
		expectedChargedFourMomentum += particleFourMomentum;
		for ( unsigned int i_p = 0 ; i_p < associatedParticles.size() ; ++i_p )
		{
			if ( recoChargedDecayProducts[ i_par ] == associatedParticles[ i_p ] ) perfectAssociation = true;
		}
		if ( perfectAssociation )
		{
			expectedChargedFourMomentumFromVertexing += particleFourMomentum;
		}
		else
		{
			expectedChargedFourMomentumFromPVA += particleFourMomentum;
		}
	}
	m_expectedNeutralEnergy.push_back( expectedNeutralFourMomentum.E() );
	m_expectedNeutralEnergyFromVertexing.push_back( expectedNeutralFourMomentumFromVertexing.E() );
	m_expectedNeutralEnergyFromPVA.push_back( expectedNeutralFourMomentumFromPVA.E() );
	m_expectedChargedEnergy.push_back( expectedChargedFourMomentum.E() );
	m_expectedChargedEnergyFromVertexing.push_back( expectedChargedFourMomentumFromVertexing.E() );
	m_expectedChargedEnergyFromPVA.push_back( expectedChargedFourMomentumFromPVA.E() );
	m_expectedNeutralMomentum.push_back( ( expectedNeutralFourMomentum.Vect() ).Mag() );
	m_expectedNeutralMomentumFromVertexing.push_back( ( expectedNeutralFourMomentumFromVertexing.Vect() ).Mag() );
	m_expectedNeutralMomentumFromPVA.push_back( ( expectedNeutralFourMomentumFromPVA.Vect() ).Mag() );
	m_expectedChargedMomentum.push_back( ( expectedChargedFourMomentum.Vect() ).Mag() );
	m_expectedChargedMomentumFromVertexing.push_back( ( expectedChargedFourMomentumFromVertexing.Vect() ).Mag() );
	m_expectedChargedMomentumFromPVA.push_back( ( expectedChargedFourMomentumFromPVA.Vect() ).Mag() );
/*
	TVector3 assignedNeutralP( 0.0 , 0.0 , 0.0 );
	TVector3 correctAssignedNeutralP( 0.0 , 0.0 , 0.0 );
	TVector3 wrongAssignedNeutralP( 0.0 , 0.0 , 0.0 );
	TVector3 assignedChargedP( 0.0 , 0.0 , 0.0 );
	TVector3 correctAssignedChargedP( 0.0 , 0.0 , 0.0 );
	TVector3 wrongAssignedChargedP( 0.0 , 0.0 , 0.0 );
	TVector3 neutralP( 0.0 , 0.0 , 0.0 );
	TVector3 associatedNeutralP( 0.0 , 0.0 , 0.0 );
	TVector3 missedNeutralP( 0.0 , 0.0 , 0.0 );
	TVector3 chargedP( 0.0 , 0.0 , 0.0 );
	TVector3 associatedChargedP( 0.0 , 0.0 , 0.0 );
	TVector3 missedChargedP( 0.0 , 0.0 , 0.0 );
	for ( unsigned int i_par = 0 ; i_par < decayProducts.size() ; ++i_par )
	{
		bool perfectAssociation = false;
		TVector3 assignedP = TVector3( decayProducts[ i_par ]->getMomentum() );
		if ( decayProducts[ i_par ]->getTracks().size() == 0 )
		{
			assignedNeutralEnergy += decayProducts[ i_par ]->getEnergy();
			assignedNeutralP += assignedP;
			for ( unsigned int i_p = 0 ; i_p < recoNeutralDecayProducts.size() ; ++i_p )
			{
				if ( decayProducts[ i_par ] == recoNeutralDecayProducts[ i_p ] ) perfectAssociation = true;
			}
			if ( perfectAssociation )
			{
				correctAssignedNeutralEnergy += decayProducts[ i_par ]->getEnergy();
				correctAssignedNeutralP += assignedP;
			}
			else
			{
				wrongAssignedNeutralEnergy += decayProducts[ i_par ]->getEnergy();
				wrongAssignedNeutralP += assignedP;
			}
		}
		else
		{
			assignedChargedEnergy += decayProducts[ i_par ]->getEnergy();
			assignedChargedP += assignedP;
			for ( unsigned int i_p = 0 ; i_p < recoChargedDecayProducts.size() ; ++i_p )
			{
				if ( decayProducts[ i_par ] == recoChargedDecayProducts[ i_p ] ) perfectAssociation = true;
			}
			if ( perfectAssociation )
			{
				correctAssignedChargedEnergy += decayProducts[ i_par ]->getEnergy();
				correctAssignedChargedP += assignedP;
			}
			else
			{
				wrongAssignedChargedEnergy += decayProducts[ i_par ]->getEnergy();
				wrongAssignedChargedP += assignedP;
			}
		}
	}
	for ( unsigned int i_p = 0 ; i_p < recoNeutralDecayProducts.size() ; ++i_p )
	{
		bool particleAssigned = false;
		TVector3 assignedP = TVector3( recoNeutralDecayProducts[ i_p ]->getMomentum() );
		neutralEnergy += recoNeutralDecayProducts[ i_p ]->getEnergy();
		neutralP += assignedP;
		for ( unsigned int i_par = 0 ; i_par < decayProducts.size() ; ++i_par )
		{
			if ( decayProducts[ i_par ] == recoNeutralDecayProducts[ i_p ] ) particleAssigned = true;
		}
		if ( particleAssigned )
		{
			associatedNeutralEnergy += recoNeutralDecayProducts[ i_p ]->getEnergy();
			associatedNeutralP += assignedP;
		}
		else
		{
			missedNeutralEnergy += recoNeutralDecayProducts[ i_p ]->getEnergy();
			missedNeutralP += assignedP;
		}
	}
	for ( unsigned int i_p = 0 ; i_p < recoChargedDecayProducts.size() ; ++i_p )
	{
		bool particleAssigned = false;
		TVector3 assignedP = TVector3( recoChargedDecayProducts[ i_p ]->getMomentum() );
		chargedEnergy += recoChargedDecayProducts[ i_p ]->getEnergy();
		chargedP += assignedP;
		for ( unsigned int i_par = 0 ; i_par < decayProducts.size() ; ++i_par )
		{
			if ( decayProducts[ i_par ] == recoChargedDecayProducts[ i_p ] ) particleAssigned = true;
		}
		if ( particleAssigned )
		{
			associatedChargedEnergy += recoChargedDecayProducts[ i_p ]->getEnergy();
			associatedChargedP += assignedP;
		}
		else
		{
			missedChargedEnergy += recoChargedDecayProducts[ i_p ]->getEnergy();
			missedChargedP += assignedP;
		}
	}

	m_assignedNeutralEnergy.push_back( assignedNeutralEnergy );
	m_assignedNeutralMomentum.push_back( assignedNeutralP.Mag() );
	m_assignedChargedEnergy.push_back( assignedChargedEnergy );
	m_assignedChargedMomentum.push_back( assignedChargedP.Mag() );
	m_correctAssignedNeutralEnergy.push_back( correctAssignedNeutralEnergy );
	m_correctAssignedNeutralMomentum.push_back( correctAssignedNeutralP.Mag() );
	m_correctAssignedChargedEnergy.push_back( correctAssignedChargedEnergy );
	m_correctAssignedChargedMomentum.push_back( correctAssignedChargedP.Mag() );
	m_wrongAssignedNeutralEnergy.push_back( wrongAssignedNeutralEnergy );
	m_wrongAssignedNeutralMomentum.push_back( wrongAssignedNeutralP.Mag() );
	m_wrongAssignedChargedEnergy.push_back( correctAssignedChargedEnergy );
	m_wrongAssignedChargedMomentum.push_back( wrongAssignedChargedP.Mag() );

	m_neutralEnergy.push_back( neutralEnergy );
	m_neutralMomentum.push_back( neutralP.Mag() );
	m_chargedEnergy.push_back( chargedEnergy );
	m_chargedMomentum.push_back( chargedP.Mag() );
	m_associatedNeutralEnergy.push_back( associatedNeutralEnergy );
	m_associatedNeutralMomentum.push_back( associatedNeutralP.Mag() );
	m_associatedChargedEnergy.push_back( associatedChargedEnergy );
	m_associatedChargedMomentum.push_back( associatedChargedP.Mag() );
	m_missedNeutralEnergy.push_back( missedNeutralEnergy );
	m_missedNeutralMomentum.push_back( missedNeutralP.Mag() );
	m_missedChargedEnergy.push_back( missedChargedEnergy );
	m_missedChargedMomentum.push_back( missedChargedP.Mag() );
*/


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////									    ////
////		Checking results and preparing solutions in LCIO format	    ////
////									    ////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


	if ( fabs( recoNeutrinoFourMomentumPos.E() - trueNeutrinoFourMomentum.E() ) < fabs( recoNeutrinoFourMomentumNeg.E() - trueNeutrinoFourMomentum.E() ) )
	{
		recoNeutrinoFourMomentumClose = recoNeutrinoFourMomentumPos;
		for ( int i_Element = 0 ; i_Element < 10 ; ++i_Element )
		{
			NeutrinoCovMat[ i_Element ] = NeutrinoCovMatPos[ i_Element ];
		}
	}
	else
	{
		recoNeutrinoFourMomentumClose = recoNeutrinoFourMomentumNeg;
		for ( int i_Element = 0 ; i_Element < 10 ; ++i_Element )
		{
			NeutrinoCovMat[ i_Element ] = NeutrinoCovMatNeg[ i_Element ];
		}
	}
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

	if ( m_displayEvent )
	{
		DDMarlinCED::draw( this , 1); // draw everything
	}
	recoHadronFourMomentum = recoVisibleFourMomentum + recoNeutrinoFourMomentumClose;
	fillTrueRecoFourMomentum( trueNeutralFourMomentum , trueChargedFourMomentum , trueLeptonFourMomentum , trueVisibleFourMomentum , trueNeutrinoFourMomentum , trueHadronFourMomentum , recoNeutralFourMomentum , recoChargedFourMomentum , recoLeptonFourMomentum , recoVisibleFourMomentum , recoNeutrinoFourMomentumClose , recoHadronFourMomentum , neutralFourMomentum , chargedFourMomentum , leptonFourMomentum , visibleFourMomentum , cheatedPVARecoChargedFourMomentum , cheatedPVARecoNeutralFourMomentum );
/*
	VertexImpl* semiLeptonicVertex = new VertexImpl;
	ReconstructedParticleImpl* semiLeptonicVertexRecoParticle = new ReconstructedParticleImpl;
	ReconstructedParticleImpl* recoNeutrinoZero = new ReconstructedParticleImpl;
	ReconstructedParticleImpl* recoNeutrinoPos = new ReconstructedParticleImpl;
	ReconstructedParticleImpl* recoNeutrinoNeg = new ReconstructedParticleImpl;
*/

	streamlog_out(DEBUG8) << "		Semi-leptonic decay type: " << SLDStatus << std::endl;

	double MomentumPos[3]{ recoNeutrinoFourMomentumPos.Px() , recoNeutrinoFourMomentumPos.Py() , recoNeutrinoFourMomentumPos.Pz() };
	recoNeutrinoPos->setType( -1 * SLDLepton->getPDG() + SLDLepton->getCharge() );
	recoNeutrinoPos->setMomentum( MomentumPos );
	recoNeutrinoPos->setEnergy( recoNeutrinoFourMomentumPos.E() );
	recoNeutrinoPos->setCovMatrix( NeutrinoCovMatPos );
	recoNeutrinoPos->setMass( 0.0 );
	recoNeutrinoPos->setCharge( 0.0 );
	recoNeutrinoPos->setReferencePoint( linkedRecoLepton->getReferencePoint() );
	for ( unsigned int j = 0 ; j < linkedRecoLepton->getParticleIDs().size() ; ++j )
	{
		ParticleIDImpl* inPID = dynamic_cast<ParticleIDImpl*>( linkedRecoLepton->getParticleIDs()[ j ] );
		ParticleIDImpl* outPID = new ParticleIDImpl;
		outPID->setType( -1 * inPID->getType() + linkedRecoLepton->getCharge() );
		outPID->setPDG( -1 * inPID->getType() + linkedRecoLepton->getCharge() );
		outPID->setLikelihood( inPID->getLikelihood() );
		outPID->setAlgorithmType( inPID->getAlgorithmType() ) ;
		for ( unsigned int k = 0 ; k < inPID->getParameters().size() ; ++k ) outPID->addParameter( inPID->getParameters()[ k ] );
		recoNeutrinoPos->addParticleID( outPID );
	}
	recoNeutrinoPos->setParticleIDUsed( linkedRecoLepton->getParticleIDUsed() );
	recoNeutrinoPos->setGoodnessOfPID( linkedRecoLepton->getGoodnessOfPID() );
	recoNeutrinoPos->setStartVertex( semiLeptonicVertex );

	double MomentumNeg[3]{ recoNeutrinoFourMomentumNeg.Px() , recoNeutrinoFourMomentumNeg.Py() , recoNeutrinoFourMomentumNeg.Pz() };
	recoNeutrinoNeg->setType( -1 * SLDLepton->getPDG() + SLDLepton->getCharge() );
	recoNeutrinoNeg->setMomentum( MomentumNeg );
	recoNeutrinoNeg->setEnergy( recoNeutrinoFourMomentumNeg.E() );
	recoNeutrinoNeg->setCovMatrix( NeutrinoCovMatNeg );
	recoNeutrinoNeg->setMass( 0.0 );
	recoNeutrinoNeg->setCharge( 0.0 );
	recoNeutrinoNeg->setReferencePoint( linkedRecoLepton->getReferencePoint() );
	for ( unsigned int j = 0 ; j < linkedRecoLepton->getParticleIDs().size() ; ++j )
	{
		ParticleIDImpl* inPID = dynamic_cast<ParticleIDImpl*>( linkedRecoLepton->getParticleIDs()[ j ] );
		ParticleIDImpl* outPID = new ParticleIDImpl;
		outPID->setType( -1 * inPID->getType() + linkedRecoLepton->getCharge() );
		outPID->setPDG( -1 * inPID->getType() + linkedRecoLepton->getCharge() );
		outPID->setLikelihood( inPID->getLikelihood() );
		outPID->setAlgorithmType( inPID->getAlgorithmType() ) ;
		for ( unsigned int k = 0 ; k < inPID->getParameters().size() ; ++k ) outPID->addParameter( inPID->getParameters()[ k ] );
		recoNeutrinoNeg->addParticleID( outPID );
	}
	recoNeutrinoNeg->setParticleIDUsed( linkedRecoLepton->getParticleIDUsed() );
	recoNeutrinoNeg->setGoodnessOfPID( linkedRecoLepton->getGoodnessOfPID() );
	recoNeutrinoNeg->setStartVertex( semiLeptonicVertex );

	double MomentumZero[3]{ 0.0 , 0.0 , 0.0 };
	recoNeutrinoZero->setType( -1 * SLDLepton->getPDG() + SLDLepton->getCharge() );
	recoNeutrinoZero->setMomentum( MomentumZero );
	recoNeutrinoZero->setEnergy( 0.0 );
	recoNeutrinoZero->setCovMatrix( NeutrinoCovMatZero );
	recoNeutrinoZero->setMass( 0.0 );
	recoNeutrinoZero->setCharge( 0.0 );
	recoNeutrinoZero->setReferencePoint( linkedRecoLepton->getReferencePoint() );
	for ( unsigned int j = 0 ; j < linkedRecoLepton->getParticleIDs().size() ; ++j )
	{
		ParticleIDImpl* inPID = dynamic_cast<ParticleIDImpl*>( linkedRecoLepton->getParticleIDs()[ j ] );
		ParticleIDImpl* outPID = new ParticleIDImpl;
		outPID->setType( -1 * inPID->getType() + linkedRecoLepton->getCharge() );
		outPID->setPDG( -1 * inPID->getType() + linkedRecoLepton->getCharge() );
		outPID->setLikelihood( inPID->getLikelihood() );
		outPID->setAlgorithmType( inPID->getAlgorithmType() ) ;
		for ( unsigned int k = 0 ; k < inPID->getParameters().size() ; ++k ) outPID->addParameter( inPID->getParameters()[ k ] );
		recoNeutrinoZero->addParticleID( outPID );
	}
	recoNeutrinoZero->setParticleIDUsed( linkedRecoLepton->getParticleIDUsed() );
	recoNeutrinoZero->setGoodnessOfPID( linkedRecoLepton->getGoodnessOfPID() );
	recoNeutrinoZero->addParticle( recoNeutrinoPos );
	recoNeutrinoZero->addParticle( recoNeutrinoNeg );
	recoNeutrinoZero->setStartVertex( semiLeptonicVertex );

	streamlog_out(DEBUG8) << "	Creating semi-leptonic vertex LCObject " << std::endl;
	float vertexPosition[ 3 ]{};
	if ( SLDStatus == 4 || SLDStatus == 5 )
	{
		vertexPosition[ 0 ] = sldVertexPosition[ 0 ];
		vertexPosition[ 1 ] = sldVertexPosition[ 1 ];
		vertexPosition[ 2 ] = sldVertexPosition[ 2 ];
	}
	else
	{
		vertexPosition[ 0 ] = primaryVertex->getPosition()[ 0 ];
		vertexPosition[ 1 ] = primaryVertex->getPosition()[ 1 ];
		vertexPosition[ 2 ] = primaryVertex->getPosition()[ 2 ];
	}
	streamlog_out(DEBUG8) << "		Position: ( " << vertexPosition[ 0 ] << " , " << vertexPosition[ 1 ] << " , " << vertexPosition[ 2 ] << " )" << std::endl;
	semiLeptonicVertex->setPrimary( false );
	streamlog_out(DEBUG8) << "		IsPrimary? " << semiLeptonicVertex->isPrimary() << std::endl;
	if ( SLDStatus == 4 )
	{
		semiLeptonicVertex->setAlgorithmType( "LepIn2ndVertex" );
	}
	else if ( SLDStatus == 5 )
	{
		semiLeptonicVertex->setAlgorithmType( "Lep+3rdVertex" );
	}
	streamlog_out(DEBUG8) << "		Algorithm Type : " << semiLeptonicVertex->getAlgorithmType() << std::endl;
	semiLeptonicVertex->setChi2( 0.0 );
	streamlog_out(DEBUG8) << "		Chi2 : " << semiLeptonicVertex->getChi2() << std::endl;
	semiLeptonicVertex->setProbability( 0.0 );
	streamlog_out(DEBUG8) << "		Probability : " << semiLeptonicVertex->getProbability() << std::endl;
	semiLeptonicVertex->setPosition( vertexPosition );
	streamlog_out(DEBUG8) << "		Position : " << semiLeptonicVertex->getPosition()[ 0 ] << " , " << semiLeptonicVertex->getPosition()[ 1 ] << " , " << semiLeptonicVertex->getPosition()[ 2 ] << std::endl;
	semiLeptonicVertex->setAssociatedParticle( semiLeptonicVertexRecoParticle );
	streamlog_out(DEBUG8) << "		AssociatedParticle : " << semiLeptonicVertex->getAssociatedParticle() << std::endl;

	double VisMomentumSLD[3]{ visibleFourMomentum.Px() , visibleFourMomentum.Py() , visibleFourMomentum.Pz() };
//	semiLeptonicVertexRecoParticle->setType( -1 * SLDLepton->getPDG() + SLDLepton->getCharge() );
	semiLeptonicVertexRecoParticle->setMomentum( VisMomentumSLD );
	semiLeptonicVertexRecoParticle->setEnergy( visibleFourMomentum.E() );
	semiLeptonicVertexRecoParticle->setCovMatrix( NeutrinoCovMat );
	semiLeptonicVertexRecoParticle->setMass( visibleFourMomentum.M() );
	semiLeptonicVertexRecoParticle->setCharge( 0.0 );
	semiLeptonicVertexRecoParticle->setReferencePoint( linkedRecoLepton->getReferencePoint() );
/*	for ( unsigned int j = 0 ; j < linkedRecoLepton->getParticleIDs().size() ; ++j )
	{
		ParticleIDImpl* inPID = dynamic_cast<ParticleIDImpl*>( linkedRecoLepton->getParticleIDs()[ j ] );
		ParticleIDImpl* outPID = new ParticleIDImpl;
		outPID->setType( -1 * inPID->getType() + linkedRecoLepton->getCharge() );
		outPID->setPDG( -1 * inPID->getType() + linkedRecoLepton->getCharge() );
		outPID->setLikelihood( inPID->getLikelihood() );
		outPID->setAlgorithmType( inPID->getAlgorithmType() ) ;
		for ( unsigned int k = 0 ; k < inPID->getParameters().size() ; ++k ) outPID->addParameter( inPID->getParameters()[ k ] );
		semiLeptonicVertexRecoParticle->addParticleID( outPID );
	}
*/
	semiLeptonicVertexRecoParticle->addParticle( recoNeutrinoZero );
	semiLeptonicVertexRecoParticle->addParticle( linkedRecoLepton );
	for ( unsigned int i_par = 0 ; i_par < decayProducts.size() ; ++i_par )
	{
		if ( decayProducts[ i_par ] != linkedRecoLepton ) semiLeptonicVertexRecoParticle->addParticle( decayProducts[ i_par ] );
	}

//	semiLeptonicVertexRecoParticle->setParticleIDUsed( linkedRecoLepton->getParticleIDUsed() );
//	semiLeptonicVertexRecoParticle->setGoodnessOfPID( linkedRecoLepton->getGoodnessOfPID() );
	semiLeptonicVertexRecoParticle->setStartVertex( startVertex );

	semiLeptonicVertices.push_back( semiLeptonicVertex );
	semiLeptonicVertexRecoParticles.push_back( semiLeptonicVertexRecoParticle );
	neutrinos.push_back( recoNeutrinoZero );
	neutrinos.push_back( recoNeutrinoPos );
	neutrinos.push_back( recoNeutrinoNeg );
}

void SLDCorrection::getCovMatPVA(	std::vector<EVENT::ReconstructedParticle*> decayProducts ,
						std::vector<EVENT::ReconstructedParticle*> associatedParticles ,
						int SLDStatus , std::vector< float > &CovMatrixPVA )
{
	float sigmaE_NeutralPVA = 0.0;
	float sigmaE_ChargedPVA = 0.0;
	std::vector< float > CovMatrixChargedPVA( 10, 0.0 );
	std::vector< float > CovMatrixNeutralPVA( 10, 0.0 );
	TLorentzVector chargedTLV( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector chargedTLVFromVertexing( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector chargedTLVFromPVA( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector neutralTLV( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector neutralTLVFromVertexing( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector neutralTLVFromPVA( 0.0 , 0.0 , 0.0 , 0.0 );

	float chargedEnergyFromPVA = 0.0;
	float chargedMomentumFromPVA = 0.0;
	float neutralEnergyFromPVA = 0.0;
	float neutralMomentumFromPVA = 0.0;
	TVector3 neutralPfromPVA = TVector3( 0.0 , 0.0 , 0.0 );
	TVector3 neutralPfromVertexing = TVector3( 0.0 , 0.0 , 0.0 );
	TVector3 chargedPfromPVA = TVector3( 0.0 , 0.0 , 0.0 );
	TVector3 chargedPfromVertexing = TVector3( 0.0 , 0.0 , 0.0 );

	for ( unsigned int i_par = 0 ; i_par < decayProducts.size() ; ++i_par )
	{
		bool perfectAssociation = false;
		TLorentzVector particleTLV( decayProducts[ i_par ]->getMomentum() , decayProducts[ i_par ]->getEnergy() );
		if ( decayProducts[ i_par ]->getTracks().size() == 0 )
		{
			neutralTLV += particleTLV;
			for ( unsigned int i_p = 0 ; i_p < associatedParticles.size() ; ++i_p )
			{
				if ( decayProducts[ i_par ] == associatedParticles[ i_p ] ) perfectAssociation = true;
			}
			if ( perfectAssociation )
			{
				neutralTLVFromVertexing += particleTLV;
			}
			else
			{
				neutralTLVFromPVA += particleTLV;
			}
		}
		else
		{
			chargedTLV += particleTLV;
			for ( unsigned int i_p = 0 ; i_p < associatedParticles.size() ; ++i_p )
			{
				if ( decayProducts[ i_par ] == associatedParticles[ i_p ] ) perfectAssociation = true;
			}
			if ( perfectAssociation )
			{
				chargedTLVFromVertexing += particleTLV;
			}
			else
			{
				chargedTLVFromPVA += particleTLV;
			}
		}
	}
	m_neutralEnergy.push_back( neutralTLV.E() );
	m_neutralEnergyFromVertexing.push_back( neutralTLVFromVertexing.E() );
	m_neutralEnergyFromPVA.push_back( neutralTLVFromPVA.E() );
	m_neutralMomentum.push_back( ( neutralTLV.Vect() ).Mag() );
	m_neutralMomentumFromVertexing.push_back( ( neutralTLVFromVertexing.Vect() ).Mag() );
	m_neutralMomentumFromPVA.push_back( ( neutralTLVFromPVA.Vect() ).Mag() );

	m_chargedEnergy.push_back( chargedTLV.E() );
	m_chargedEnergyFromVertexing.push_back( chargedTLVFromVertexing.E() );
	m_chargedEnergyFromPVA.push_back( chargedTLVFromPVA.E() );
	m_chargedMomentum.push_back( ( chargedTLV.Vect() ).Mag() );
	m_chargedMomentumFromVertexing.push_back( ( chargedTLVFromVertexing.Vect() ).Mag() );
	m_chargedMomentumFromPVA.push_back( ( chargedTLVFromPVA.Vect() ).Mag() );

	float chargedTheta = chargedPfromPVA.Theta();
	float chargedPhi = chargedPfromPVA.Phi();
	float neutralTheta = neutralPfromPVA.Theta();
	float neutralPhi = neutralPfromPVA.Phi();
	if ( SLDStatus == 4 )
	{
		sigmaE_NeutralPVA = ( m_cheatPVAneutral ? 0.0 : 6.7 );
		sigmaE_ChargedPVA = ( m_cheatPVAcharged ? 0.0 : 9.2 );
	}
	else if ( SLDStatus == 5 )
	{
		sigmaE_NeutralPVA = ( m_cheatPVAneutral ? 0.0 : 6.5 );
		sigmaE_ChargedPVA = ( m_cheatPVAcharged ? 0.0 : 9.5 );
	}
	CovMatrixChargedPVA[ 0 ] = pow( sin( chargedTheta ) , 2 ) * pow( cos( chargedPhi ) , 2 );
	CovMatrixChargedPVA[ 1 ] = pow( sin( chargedTheta ) , 2 ) * sin( chargedPhi ) * cos( chargedPhi );
	CovMatrixChargedPVA[ 2 ] = pow( sin( chargedTheta ) , 2 ) * pow( sin( chargedPhi ) , 2 );
	CovMatrixChargedPVA[ 3 ] = sin( chargedTheta ) * cos( chargedTheta ) * cos( chargedPhi );
	CovMatrixChargedPVA[ 4 ] = sin( chargedTheta ) * cos( chargedTheta ) * sin( chargedPhi );
	CovMatrixChargedPVA[ 5 ] = pow( cos( chargedTheta ) , 2 );
	CovMatrixChargedPVA[ 6 ] = ( chargedMomentumFromPVA / chargedEnergyFromPVA ) * sin( chargedTheta ) * cos( chargedPhi );
	CovMatrixChargedPVA[ 7 ] = ( chargedMomentumFromPVA / chargedEnergyFromPVA ) * sin( chargedTheta ) * sin( chargedPhi );
	CovMatrixChargedPVA[ 8 ] = ( chargedMomentumFromPVA / chargedEnergyFromPVA ) * cos( chargedTheta );
	CovMatrixChargedPVA[ 9 ] = pow( chargedMomentumFromPVA , 2 ) / pow( chargedEnergyFromPVA , 2 );
	CovMatrixNeutralPVA[ 0 ] = pow( sin( neutralTheta ) , 2 ) * pow( cos( neutralPhi ) , 2 );
	CovMatrixNeutralPVA[ 1 ] = pow( sin( neutralTheta ) , 2 ) * sin( neutralPhi ) * cos( neutralPhi );
	CovMatrixNeutralPVA[ 2 ] = pow( sin( neutralTheta ) , 2 ) * pow( sin( neutralPhi ) , 2 );
	CovMatrixNeutralPVA[ 3 ] = sin( neutralTheta ) * cos( neutralTheta ) * cos( neutralPhi );
	CovMatrixNeutralPVA[ 4 ] = sin( neutralTheta ) * cos( neutralTheta ) * sin( neutralPhi );
	CovMatrixNeutralPVA[ 5 ] = pow( cos( neutralTheta ) , 2 );
	CovMatrixNeutralPVA[ 6 ] = ( neutralMomentumFromPVA / neutralEnergyFromPVA ) * sin( neutralTheta ) * cos( neutralPhi );
	CovMatrixNeutralPVA[ 7 ] = ( neutralMomentumFromPVA / neutralEnergyFromPVA ) * sin( neutralTheta ) * sin( neutralPhi );
	CovMatrixNeutralPVA[ 8 ] = ( neutralMomentumFromPVA / neutralEnergyFromPVA ) * cos( neutralTheta );
	CovMatrixNeutralPVA[ 9 ] = pow( neutralMomentumFromPVA , 2 ) / pow( neutralEnergyFromPVA , 2 );
	streamlog_out(DEBUG9) << "	sigmaE_NeutralPVA = " << sigmaE_NeutralPVA << std::endl;
	streamlog_out(DEBUG9) << "	neutralMomentumFromPVA = " << neutralMomentumFromPVA << std::endl;
	streamlog_out(DEBUG9) << "	neutralEnergyFromPVA = " << neutralEnergyFromPVA << std::endl;
	streamlog_out(DEBUG9) << "	sigmaE_ChargedPVA = " << sigmaE_ChargedPVA << std::endl;
	streamlog_out(DEBUG9) << "	chargedMomentumFromPVA = " << chargedMomentumFromPVA << std::endl;
	streamlog_out(DEBUG9) << "	chargedEnergyFromPVA = " << chargedEnergyFromPVA << std::endl;
	for ( int i_Element = 0 ; i_Element < 10 ; ++i_Element )
	{
		CovMatrixPVA.push_back( pow( neutralEnergyFromPVA / neutralMomentumFromPVA , 2 ) * sigmaE_NeutralPVA * CovMatrixNeutralPVA[ i_Element ] + pow( chargedEnergyFromPVA / chargedMomentumFromPVA , 2 ) * sigmaE_ChargedPVA * CovMatrixChargedPVA[ i_Element ] );
	}
	streamlog_out(DEBUG9) << "	CovMatPVA :	" << CovMatrixPVA[ 0 ] << std::endl;
	streamlog_out(DEBUG9) << "			" << CovMatrixPVA[ 1 ] << "	,	" << CovMatrixPVA[ 2 ] << std::endl;
	streamlog_out(DEBUG9) << "			" << CovMatrixPVA[ 3 ] << "	,	" << CovMatrixPVA[ 4 ] << "	,	" << CovMatrixPVA[ 5 ] << std::endl;
	streamlog_out(DEBUG9) << "			" << CovMatrixPVA[ 6 ] << "	,	" << CovMatrixPVA[ 7 ] << "	,	" << CovMatrixPVA[ 8 ] << "	,	" << CovMatrixPVA[ 9 ] << std::endl;
	streamlog_out(DEBUG9) << "" << std::endl;

}

void SLDCorrection::getCovMatFlightDirection(	TVector3 flightDirection , float sigmaTheta , float sigmaPhi , std::vector< float > &CovMatrixFlightDirection )
{
	//	Obtain covariance matrix on flight direction (ux,uy,uz) from the
	//	flight direction errors (angular errors: sigmaTheta and sigmaPhi).
	//
	//	define the jacobian as the 2x3 matrix:
	//
	//
	//
	//			Dux/DTheta		Duy/DTheta		Duz/DTheta
	//	 J =
	//			Dux/DPhi		Duy/DPhi		Duz/DPhi
	//
	//
	//
	//			cosTheta.cosPhi		cosTheta.sinPhi		-sinPhi
	//	 J =
	//			-sinTheta.sinPhi	sinTheta.cosPhi		0
	//
	//
	//
	//	Order in the covariance matrix on helix parameters:
	//
	//			Theta.Theta		Theta.phi
	//	Cov =
	//			phi.Theta		phi.phi
	//
	//
	//
	CovMatrixFlightDirection.clear();
	const int rows			= 2; // n rows jacobian
	const int columns		= 3; // n columns jacobian
	const int kspace_dim		= 3;

	TMatrixD covMatrixFlightDir( kspace_dim , kspace_dim );
	double jacobian_by_rows[rows*columns] =
	{
		cos( flightDirection.Theta() ) * cos( flightDirection.Phi() )		,	cos( flightDirection.Theta() ) * sin( flightDirection.Phi() )	,	-1.0 * sin( flightDirection.Phi() )	,
		-1.0 * sin( flightDirection.Theta() ) * sin( flightDirection.Phi() )	,	sin( flightDirection.Theta() ) * cos( flightDirection.Phi() )	,	0.0
	};
	TMatrixD jacobian(rows,columns, jacobian_by_rows, "C");
	streamlog_out(DEBUG9) << "	JacobianFlightDir :	" << jacobian( 0 , 0 ) << "	,	" << jacobian( 0 , 1 ) << "	,	" << jacobian( 0 , 2 ) << std::endl;
	streamlog_out(DEBUG9) << "				" << jacobian( 1 , 0 ) << "	,	" << jacobian( 1 , 1 ) << "	,	" << jacobian( 1 , 2 ) << std::endl;
	double angular_cov_matrix_by_rows[rows*rows] =
	{
		pow( sigmaTheta , 2 )	,	0.0			,
		0.0			,	pow( sigmaPhi , 2 )
	};
	streamlog_out(DEBUG9) << "	angular_cov_matrix_by_rows :	" << angular_cov_matrix_by_rows[ 0 ] << "	,	" << angular_cov_matrix_by_rows[ 1 ] << std::endl;
	streamlog_out(DEBUG9) << "					" << angular_cov_matrix_by_rows[ 2 ] << "	,	" << angular_cov_matrix_by_rows[ 3 ] << std::endl;
	TMatrixD covMatrix_FD(rows,rows, angular_cov_matrix_by_rows, "C");
	streamlog_out(DEBUG9) << "	covMatrix_FD :	" << covMatrix_FD( 0 , 0 ) << "	,	" << covMatrix_FD( 0 , 1 ) << std::endl;
	streamlog_out(DEBUG9) << "			" << covMatrix_FD( 1 , 0 ) << "	,	" << covMatrix_FD( 1 , 1 ) << std::endl;
	covMatrixFlightDir.Mult( TMatrixD( jacobian , TMatrixD::kTransposeMult , covMatrix_FD ) , jacobian );
	CovMatrixFlightDirection.push_back( covMatrixFlightDir( 0 , 0 ) );
	CovMatrixFlightDirection.push_back( covMatrixFlightDir( 1 , 0 ) );
	CovMatrixFlightDirection.push_back( covMatrixFlightDir( 1 , 1 ) );
	CovMatrixFlightDirection.push_back( covMatrixFlightDir( 2 , 0 ) );
	CovMatrixFlightDirection.push_back( covMatrixFlightDir( 2 , 1 ) );
	CovMatrixFlightDirection.push_back( covMatrixFlightDir( 2 , 2 ) );
	streamlog_out(DEBUG9) << "	Theta		= " << flightDirection.Theta() << std::endl;
	streamlog_out(DEBUG9) << "	Phi		= " << flightDirection.Phi() << std::endl;
	streamlog_out(DEBUG9) << "	sigmaTheta	= " << sigmaTheta << std::endl;
	streamlog_out(DEBUG9) << "	sigmaPhi	= " << sigmaPhi << std::endl;
	streamlog_out(DEBUG9) << "	CovMatFlightDir :	" << CovMatrixFlightDirection[ 0 ] << std::endl;
	streamlog_out(DEBUG9) << "				" << CovMatrixFlightDirection[ 1 ] << "	,	" << CovMatrixFlightDirection[ 2 ] << std::endl;
	streamlog_out(DEBUG9) << "				" << CovMatrixFlightDirection[ 3 ] << "	,	" << CovMatrixFlightDirection[ 4 ] << "	,	" << CovMatrixFlightDirection[ 5 ] << std::endl;
	streamlog_out(DEBUG9) << "" << std::endl;
}

void SLDCorrection::getCovMatDetFlightDirection(	std::vector<EVENT::ReconstructedParticle*> decayProducts ,
							TVector3 flightDirection , std::vector< float > CovMatrixFlightDirection ,
							std::vector< float > &CovMatrixDetector , EVENT::ReconstructedParticle* linkedRecoLepton ,
							std::vector< float > &CovMatrixDetPar , std::vector< float > &CovMatrixDetNor )
{
	std::vector< float > initialCovMatrixDetector( 10 , 0.0 );
	TLorentzVector visibleFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	CovMatrixDetector.clear();
	CovMatrixDetPar.clear();
	CovMatrixDetNor.clear();
	float takeCovMat = ( m_cheatLepton4momentum ? 0.0 : 1.0 );
	for ( int i_Element = 0 ; i_Element < 10 ; ++i_Element )
	{
		initialCovMatrixDetector[ i_Element ] = takeCovMat * linkedRecoLepton->getCovMatrix()[ i_Element ];
	}

	for ( unsigned int i_par = 0 ; i_par < decayProducts.size() ; ++i_par )
	{
		if ( decayProducts[ i_par ]->getTracks().size() == 0 )
		{
			takeCovMat = ( m_cheatNeutral4momentum ? 0.0 : 1.0 );
		}
		else
		{
			takeCovMat = ( m_cheatCharged4momentum ? 0.0 : 1.0 );
		}
		visibleFourMomentum += TLorentzVector( decayProducts[ i_par ]->getMomentum() , decayProducts[ i_par ]->getEnergy() );
		for ( int i_Element = 0 ; i_Element < 10 ; ++i_Element )
		{
			initialCovMatrixDetector[ i_Element ] += takeCovMat * decayProducts[ i_par ]->getCovMatrix()[ i_Element ];
		}
	}
	for ( int i_Element = 0 ; i_Element < 10 ; ++i_Element )
	{
		CovMatrixDetector.push_back( initialCovMatrixDetector[ i_Element ] );
	}
	streamlog_out(DEBUG9) << "	CovMatDetector :	" << CovMatrixDetector[ 0 ] << std::endl;
	streamlog_out(DEBUG9) << "				" << CovMatrixDetector[ 1 ] << "	,	" << CovMatrixDetector[ 2 ] << std::endl;
	streamlog_out(DEBUG9) << "				" << CovMatrixDetector[ 3 ] << "	,	" << CovMatrixDetector[ 4 ] << "	,	" << CovMatrixDetector[ 5 ] << std::endl;
	streamlog_out(DEBUG9) << "				" << CovMatrixDetector[ 6 ] << "	,	" << CovMatrixDetector[ 7 ] << "	,	" << CovMatrixDetector[ 8 ] << "	,	" << CovMatrixDetector[ 9 ] << std::endl;
	streamlog_out(DEBUG9) << "" << std::endl;
	getCovMatrixDetPar( flightDirection , visibleFourMomentum , CovMatrixFlightDirection , CovMatrixDetector , CovMatrixDetPar );
	getCovMatrixDetNor( flightDirection , visibleFourMomentum , CovMatrixFlightDirection , CovMatrixDetector , CovMatrixDetNor );
}

void SLDCorrection::getCovMatrixDetPar(		TVector3 flightDirection , TLorentzVector visibleFourMomentum , std::vector< float > CovMatrixFlightDirection ,
						std::vector< float > initialCovMatrixDetector , std::vector< float > &CovMatrixDetPar )
{
	//	Obtain covariance matrix of visible decay products parallel to flight direction from the
	//	flight direction covariance matrix and detector resolution.
	//
	//	P_x_par = ( Px * ux + Py * uy + Pz * uz ) * ux
	//
	//	P_y_par = ( Px * ux + Py * uy + Pz * uz ) * uy
	//
	//	P_z_par = ( Px * ux + Py * uy + Pz * uz ) * uz
	//
	//	M_transverse^2 = E^2 - p_z^2	=> E_transverse^2 = 2 * E^2 - 2 * P_par^2 - M^2
	//
	//	E_par = E - sqrt( 2 * E^2 - 2 * P_par^2 - M^2 )
	//
	//	define the jacobian as the 7x4 matrix:
	//
	//
	//
	//			DPxpar/DPx		DPypar/DPx		DPzpar/DPx		DEpar/DPx
	//
	//			DPxpar/DPy		DPypar/DPy		DPzpar/DPy		DEpar/DPy
	//
	//			DPxpar/DPz		DPypar/DPz		DPzpar/DPz		DEpar/DPz
	//
	//	J = 		DPxpar/DE		DPypar/DE		DPzpar/DE		DEpar/DE
	//
	//			DPxpar/Dux		DPypar/Dux		DPzpar/Dux		DEpar/Dux
	//
	//			DPxpar/Duy		DPypar/Duy		DPzpar/Duy		DEpar/Duy
	//
	//			DPxpar/Duz		DPypar/Duz		DPzpar/Duz		DEpar/Duz
	//
	//
	//
	//	Order in the covariance matrix on detector resolution and flight direction error:
	//
	//			Px.Px		Px.Py		Px.Pz		Px.E		Px.ux		Px.uy		Px.uz
	//
	//			Py.Px		Py.Py		Py.Pz		Py.E		Py.ux		Py.uy		Py.uz
	//
	//			Pz.Px		Pz.Py		Pz.Pz		Pz.E		Pz.ux		Pz.uy		Pz.uz
	//
	//	Cov =		E.Px		E.Py		E.Pz		E.E		E.ux		E.uy		E.uz
	//
	//			ux.Px		ux.Py		ux.Pz		ux.E		ux.ux		ux.uy		ux.uz
	//
	//			uy.Px		uy.Py		uy.Pz		uy.E		uy.ux		uy.uy		uy.uz
	//
	//			uz.Px		uz.Py		uz.Pz		uz.E		uz.ux		uz.uy		uz.uz
	//
	//
	//
	//
	CovMatrixDetPar.clear();
	const int rows			= 7; // n rows jacobian
	const int columns		= 4; // n columns jacobian
	const int SpaceTime_dim		= 4;

	TVector3 Momentum = visibleFourMomentum.Vect();
	double Energy = visibleFourMomentum.E();
	double mass = visibleFourMomentum.M();
	double P_par = Momentum.Dot( flightDirection );
	double P_nor = sqrt( Momentum.Mag2() - pow( P_par , 2 ) );
	double Px = visibleFourMomentum.Px();
	double Py = visibleFourMomentum.Py();
	double Pz = visibleFourMomentum.Pz();

	double uX = flightDirection.X();
	double uY = flightDirection.Y();
	double uZ = flightDirection.Z();
	streamlog_out(DEBUG9) << "	E_vis 	:	" << Energy << std::endl;
	streamlog_out(DEBUG9) << "	mass 	:	" << mass << std::endl;
	streamlog_out(DEBUG9) << "	P_par 	:	" << P_par << std::endl;
	streamlog_out(DEBUG9) << "	P_nor 	:	" << P_nor << std::endl;
	streamlog_out(DEBUG9) << "	Px 	:	" << Px << std::endl;
	streamlog_out(DEBUG9) << "	Py 	:	" << Py << std::endl;
	streamlog_out(DEBUG9) << "	Pz 	:	" << Pz << std::endl;
	streamlog_out(DEBUG9) << "	uX 	:	" << uX << std::endl;
	streamlog_out(DEBUG9) << "	uY 	:	" << uY << std::endl;
	streamlog_out(DEBUG9) << "	uZ 	:	" << uZ << std::endl;

	TMatrixD covMatrixFlightDirDet( SpaceTime_dim , SpaceTime_dim );
	double jacobian_by_rows[rows*columns] =
	{
		uX * uX			,	uX * uY			,	uX * uZ			,	2.0 * uX * P_par / sqrt( pow( mass , 2 ) + 2.0 * pow( P_nor , 2 ) )	,
		uY * uX			,	uY * uY			,	uY * uZ			,	2.0 * uY * P_par / sqrt( pow( mass , 2 ) + 2.0 * pow( P_nor , 2 ) )	,
		uZ * uX			,	uZ * uY			,	uZ * uZ			,	2.0 * uZ * P_par / sqrt( pow( mass , 2 ) + 2.0 * pow( P_nor , 2 ) )	,
		0.0			,	0.0			,	0.0			,	1.0 - 2.0 * Energy / sqrt( pow( mass , 2 ) + 2.0 * pow( P_nor , 2 ) )	,
		P_par + Px * uX		,	Px * uY			,	Px * uZ			,	2.0 * Px * P_par / sqrt( pow( mass , 2 ) + 2.0 * pow( P_nor , 2 ) )	,
		Py * uX			,	P_par + Py * uY		,	Py * uZ			,	2.0 * Py * P_par / sqrt( pow( mass , 2 ) + 2.0 * pow( P_nor , 2 ) )	,
		Pz * uX			,	Pz * uY			,	P_par + Pz * uZ		,	2.0 * Pz * P_par / sqrt( pow( mass , 2 ) + 2.0 * pow( P_nor , 2 ) )
	};
	TMatrixD jacobian(rows,columns, jacobian_by_rows, "C");
	streamlog_out(DEBUG9) << "	Jacobian (FD)->(p,E) :	" << jacobian( 0 , 0 ) << "	,	" << jacobian( 0 , 1 ) << "	,	" << jacobian( 0 , 2 ) << "	,	" << jacobian( 0 , 3 ) << std::endl;
	streamlog_out(DEBUG9) << "				" << jacobian( 1 , 0 ) << "	,	" << jacobian( 1 , 1 ) << "	,	" << jacobian( 1 , 2 ) << "	,	" << jacobian( 1 , 3 ) << std::endl;
	streamlog_out(DEBUG9) << "				" << jacobian( 2 , 0 ) << "	,	" << jacobian( 2 , 1 ) << "	,	" << jacobian( 2 , 2 ) << "	,	" << jacobian( 2 , 3 ) << std::endl;
	streamlog_out(DEBUG9) << "				" << jacobian( 3 , 0 ) << "	,	" << jacobian( 3 , 1 ) << "	,	" << jacobian( 3 , 2 ) << "	,	" << jacobian( 3 , 3 ) << std::endl;
	streamlog_out(DEBUG9) << "				" << jacobian( 4 , 0 ) << "	,	" << jacobian( 4 , 1 ) << "	,	" << jacobian( 4 , 2 ) << "	,	" << jacobian( 4 , 3 ) << std::endl;
	streamlog_out(DEBUG9) << "				" << jacobian( 5 , 0 ) << "	,	" << jacobian( 5 , 1 ) << "	,	" << jacobian( 5 , 2 ) << "	,	" << jacobian( 5 , 3 ) << std::endl;
	streamlog_out(DEBUG9) << "				" << jacobian( 6 , 0 ) << "	,	" << jacobian( 6 , 1 ) << "	,	" << jacobian( 6 , 2 ) << "	,	" << jacobian( 6 , 3 ) << std::endl;
	double input_cov_matrix_by_rows[rows*rows] =
	{
		initialCovMatrixDetector[ 0 ]	,	initialCovMatrixDetector[ 1 ]	,	initialCovMatrixDetector[ 3 ]	,	initialCovMatrixDetector[ 6 ]	,	0.0				,	0.0				,	0.0				,
		initialCovMatrixDetector[ 1 ]	,	initialCovMatrixDetector[ 2 ]	,	initialCovMatrixDetector[ 4 ]	,	initialCovMatrixDetector[ 7 ]	,	0.0				,	0.0				,	0.0				,
		initialCovMatrixDetector[ 3 ]	,	initialCovMatrixDetector[ 4 ]	,	initialCovMatrixDetector[ 5 ]	,	initialCovMatrixDetector[ 8 ]	,	0.0				,	0.0				,	0.0				,
		initialCovMatrixDetector[ 6 ]	,	initialCovMatrixDetector[ 7 ]	,	initialCovMatrixDetector[ 8 ]	,	initialCovMatrixDetector[ 9 ]	,	0.0				,	0.0				,	0.0				,
		0.0				,	0.0				,	0.0				,	0.0				,	CovMatrixFlightDirection[ 0 ]	,	CovMatrixFlightDirection[ 1 ]	,	CovMatrixFlightDirection[ 3 ]	,
		0.0				,	0.0				,	0.0				,	0.0				,	CovMatrixFlightDirection[ 1 ]	,	CovMatrixFlightDirection[ 2 ]	,	CovMatrixFlightDirection[ 4 ]	,
		0.0				,	0.0				,	0.0				,	0.0				,	CovMatrixFlightDirection[ 3 ]	,	CovMatrixFlightDirection[ 4 ]	,	CovMatrixFlightDirection[ 5 ]	,
	};
	TMatrixD covMatrix_DetFD(rows,rows, input_cov_matrix_by_rows, "C");
	streamlog_out(DEBUG9) << "	CovMat (p,E,ux,uy,uz) :	" << covMatrix_DetFD( 0 , 0 ) << "	, " << covMatrix_DetFD( 0 , 1 ) << "	, " << covMatrix_DetFD( 0 , 2 ) << "	, " << covMatrix_DetFD( 0 , 3 ) << "	, " << covMatrix_DetFD( 0 , 4 ) << "	, " << covMatrix_DetFD( 0 , 5 ) << "	, " << covMatrix_DetFD( 0 , 6 ) << std::endl;
	streamlog_out(DEBUG9) << "					" << covMatrix_DetFD( 1 , 0 ) << "	, " << covMatrix_DetFD( 1 , 1 ) << "	, " << covMatrix_DetFD( 1 , 2 ) << "	, " << covMatrix_DetFD( 1 , 3 ) << "	, " << covMatrix_DetFD( 1 , 4 ) << "	, " << covMatrix_DetFD( 1 , 5 ) << "	, " << covMatrix_DetFD( 1 , 6 ) << std::endl;
	streamlog_out(DEBUG9) << "					" << covMatrix_DetFD( 2 , 0 ) << "	, " << covMatrix_DetFD( 2 , 1 ) << "	, " << covMatrix_DetFD( 2 , 2 ) << "	, " << covMatrix_DetFD( 2 , 3 ) << "	, " << covMatrix_DetFD( 2 , 4 ) << "	, " << covMatrix_DetFD( 2 , 5 ) << "	, " << covMatrix_DetFD( 2 , 6 ) << std::endl;
	streamlog_out(DEBUG9) << "					" << covMatrix_DetFD( 3 , 0 ) << "	, " << covMatrix_DetFD( 3 , 1 ) << "	, " << covMatrix_DetFD( 3 , 2 ) << "	, " << covMatrix_DetFD( 3 , 3 ) << "	, " << covMatrix_DetFD( 3 , 4 ) << "	, " << covMatrix_DetFD( 3 , 5 ) << "	, " << covMatrix_DetFD( 3 , 6 ) << std::endl;
	streamlog_out(DEBUG9) << "					" << covMatrix_DetFD( 4 , 0 ) << "	, " << covMatrix_DetFD( 4 , 1 ) << "	, " << covMatrix_DetFD( 4 , 2 ) << "	, " << covMatrix_DetFD( 4 , 3 ) << "	, " << covMatrix_DetFD( 4 , 4 ) << "	, " << covMatrix_DetFD( 4 , 5 ) << "	, " << covMatrix_DetFD( 4 , 6 ) << std::endl;
	streamlog_out(DEBUG9) << "					" << covMatrix_DetFD( 5 , 0 ) << "	, " << covMatrix_DetFD( 5 , 1 ) << "	, " << covMatrix_DetFD( 5 , 2 ) << "	, " << covMatrix_DetFD( 5 , 3 ) << "	, " << covMatrix_DetFD( 5 , 4 ) << "	, " << covMatrix_DetFD( 5 , 5 ) << "	, " << covMatrix_DetFD( 5 , 6 ) << std::endl;
	streamlog_out(DEBUG9) << "					" << covMatrix_DetFD( 6 , 0 ) << "	, " << covMatrix_DetFD( 6 , 1 ) << "	, " << covMatrix_DetFD( 6 , 2 ) << "	, " << covMatrix_DetFD( 6 , 3 ) << "	, " << covMatrix_DetFD( 6 , 4 ) << "	, " << covMatrix_DetFD( 6 , 5 ) << "	, " << covMatrix_DetFD( 6 , 6 ) << std::endl;
	covMatrixFlightDirDet.Mult( TMatrixD( jacobian , TMatrixD::kTransposeMult , covMatrix_DetFD) , jacobian );
	CovMatrixDetPar.push_back( covMatrixFlightDirDet( 0 , 0 ) );
	CovMatrixDetPar.push_back( covMatrixFlightDirDet( 1 , 0 ) );
	CovMatrixDetPar.push_back( covMatrixFlightDirDet( 1 , 1 ) );
	CovMatrixDetPar.push_back( covMatrixFlightDirDet( 2 , 0 ) );
	CovMatrixDetPar.push_back( covMatrixFlightDirDet( 2 , 1 ) );
	CovMatrixDetPar.push_back( covMatrixFlightDirDet( 2 , 2 ) );
	CovMatrixDetPar.push_back( covMatrixFlightDirDet( 3 , 0 ) );
	CovMatrixDetPar.push_back( covMatrixFlightDirDet( 3 , 1 ) );
	CovMatrixDetPar.push_back( covMatrixFlightDirDet( 3 , 2 ) );
	CovMatrixDetPar.push_back( covMatrixFlightDirDet( 3 , 3 ) );
	streamlog_out(DEBUG9) << "	CovMatDetFlightPar :	" << CovMatrixDetPar[ 0 ] << std::endl;
	streamlog_out(DEBUG9) << "				" << CovMatrixDetPar[ 1 ] << "	,	" << CovMatrixDetPar[ 2 ] << std::endl;
	streamlog_out(DEBUG9) << "				" << CovMatrixDetPar[ 3 ] << "	,	" << CovMatrixDetPar[ 4 ] << "	,	" << CovMatrixDetPar[ 5 ] << std::endl;
	streamlog_out(DEBUG9) << "				" << CovMatrixDetPar[ 6 ] << "	,	" << CovMatrixDetPar[ 7 ] << "	,	" << CovMatrixDetPar[ 8 ] << "	,	" << CovMatrixDetPar[ 9 ] << std::endl;
	streamlog_out(DEBUG9) << "" << std::endl;
}

void SLDCorrection::getCovMatrixDetNor(		TVector3 flightDirection , TLorentzVector visibleFourMomentum , std::vector< float > CovMatrixFlightDirection ,
						std::vector< float > initialCovMatrixDetector , std::vector< float > &CovMatrixDetNor )
{
	//	Obtain covariance matrix of visible decay products parallel to flight direction from the
	//	flight direction covariance matrix and detector resolution.
	//
	//	P_x_nor = Px - ( Px * ux + Py * uy + Pz * uz ) * ux
	//
	//	P_y_nor = Py - ( Px * ux + Py * uy + Pz * uz ) * uy
	//
	//	P_z_nor = Pz - ( Px * ux + Py * uy + Pz * uz ) * uz
	//
	//	M_transverse^2 = E^2 - p_z^2	=> E_transverse^2 = 2 * E^2 - 2 * P_par^2 - M^2
	//
	//	define the jacobian as the 7x4 matrix:
	//
	//
	//
	//			DPxnor/DPx		DPynor/DPx		DPznor/DPx		DEnor/DPx
	//
	//			DPxnor/DPy		DPynor/DPy		DPznor/DPy		DEnor/DPy
	//
	//			DPxnor/DPz		DPynor/DPz		DPznor/DPz		DEnor/DPz
	//
	//	J = 		DPxnor/DE		DPynor/DE		DPznor/DE		DEnor/DE
	//
	//			DPxnor/Dux		DPynor/Dux		DPznor/Dux		DEnor/Dux
	//
	//			DPxnor/Duy		DPynor/Duy		DPznor/Duy		DEnor/Duy
	//
	//			DPxnor/Duz		DPynor/Duz		DPznor/Duz		DEnor/Duz
	//
	//
	//
	//	Order in the covariance matrix on detector resolution and flight direction error:
	//
	//			Px.Px		Px.Py		Px.Pz		Px.E		Px.ux		Px.uy		Px.uz
	//
	//			Py.Px		Py.Py		Py.Pz		Py.E		Py.ux		Py.uy		Py.uz
	//
	//			Pz.Px		Pz.Py		Pz.Pz		Pz.E		Pz.ux		Pz.uy		Pz.uz
	//
	//	Cov =		E.Px		E.Py		E.Pz		E.E		E.ux		E.uy		E.uz
	//
	//			ux.Px		ux.Py		ux.Pz		ux.E		ux.ux		ux.uy		ux.uz
	//
	//			uy.Px		uy.Py		uy.Pz		uy.E		uy.ux		uy.uy		uy.uz
	//
	//			uz.Px		uz.Py		uz.Pz		uz.E		uz.ux		uz.uy		uz.uz
	//
	//
	//
	//
	CovMatrixDetNor.clear();
	const int rows			= 7; // n rows jacobian
	const int columns		= 4; // n columns jacobian
	const int SpaceTime_dim		= 4;

	TVector3 Momentum = visibleFourMomentum.Vect();
	double Energy = visibleFourMomentum.E();
	double mass = visibleFourMomentum.M();
	double P_par = Momentum.Dot( flightDirection );
	double P_nor = sqrt( Momentum.Mag2() - pow( P_par , 2 ) );
	double Px = visibleFourMomentum.Px();
	double Py = visibleFourMomentum.Py();
	double Pz = visibleFourMomentum.Pz();

	double uX = flightDirection.X();
	double uY = flightDirection.Y();
	double uZ = flightDirection.Z();

	TMatrixD covMatrixFlightDirDet( SpaceTime_dim , SpaceTime_dim );
	double jacobian_by_rows[rows*columns] =
	{
		1.0 - uX * uX		,	-1.0 * uX * uY		,	-1.0 * uX * uZ		,	-2.0 * uX * P_par / sqrt( pow( mass , 2 ) + 2.0 * pow( P_nor , 2 ) )	,
		-1.0 * uY * uX		,	1.0 - uY * uY		,	-1.0 * uY * uZ		,	-2.0 * uY * P_par / sqrt( pow( mass , 2 ) + 2.0 * pow( P_nor , 2 ) )	,
		-1.0 * uZ * uX		,	-1.0 * uZ * uY		,	1.0 - uZ * uZ		,	-2.0 * uZ * P_par / sqrt( pow( mass , 2 ) + 2.0 * pow( P_nor , 2 ) )	,
		0.0			,	0.0			,	0.0			,	2.0 * Energy / sqrt( pow( mass , 2 ) + 2.0 * pow( P_nor , 2 ) )	,
		-1.0 * P_par - Px * uX	,	-1.0 * Px * uY		,	-1.0 * Px * uZ		,	-2.0 * Px * P_par / sqrt( pow( mass , 2 ) + 2.0 * pow( P_nor , 2 ) )	,
		-1.0 * Py * uX		,	-1.0 * P_par - Py * uY	,	-1.0 * Py * uZ		,	-2.0 * Py * P_par / sqrt( pow( mass , 2 ) + 2.0 * pow( P_nor , 2 ) )	,
		-1.0 * Pz * uX		,	-1.0 * Pz * uY		,	-1.0 * P_par - Pz * uZ		,	-2.0 * Pz * P_par / sqrt( pow( mass , 2 ) + 2.0 * pow( P_nor , 2 ) )
	};
	TMatrixD jacobian(rows,columns, jacobian_by_rows, "C");
	double input_cov_matrix_by_rows[rows*rows] =
	{
		initialCovMatrixDetector[ 0 ]	,	initialCovMatrixDetector[ 1 ]	,	initialCovMatrixDetector[ 3 ]	,	initialCovMatrixDetector[ 6 ]	,	0.0				,	0.0				,	0.0				,
		initialCovMatrixDetector[ 1 ]	,	initialCovMatrixDetector[ 2 ]	,	initialCovMatrixDetector[ 4 ]	,	initialCovMatrixDetector[ 7 ]	,	0.0				,	0.0				,	0.0				,
		initialCovMatrixDetector[ 3 ]	,	initialCovMatrixDetector[ 4 ]	,	initialCovMatrixDetector[ 5 ]	,	initialCovMatrixDetector[ 8 ]	,	0.0				,	0.0				,	0.0				,
		initialCovMatrixDetector[ 6 ]	,	initialCovMatrixDetector[ 7 ]	,	initialCovMatrixDetector[ 8 ]	,	initialCovMatrixDetector[ 9 ]	,	0.0				,	0.0				,	0.0				,
		0.0				,	0.0				,	0.0				,	0.0				,	CovMatrixFlightDirection[ 0 ]	,	CovMatrixFlightDirection[ 1 ]	,	CovMatrixFlightDirection[ 3 ]	,
		0.0				,	0.0				,	0.0				,	0.0				,	CovMatrixFlightDirection[ 1 ]	,	CovMatrixFlightDirection[ 2 ]	,	CovMatrixFlightDirection[ 4 ]	,
		0.0				,	0.0				,	0.0				,	0.0				,	CovMatrixFlightDirection[ 3 ]	,	CovMatrixFlightDirection[ 4 ]	,	CovMatrixFlightDirection[ 5 ]	,
	};
	TMatrixD covMatrix_DetFD(rows,rows, input_cov_matrix_by_rows, "C");
	covMatrixFlightDirDet.Mult( TMatrixD( jacobian , TMatrixD::kTransposeMult , covMatrix_DetFD) , jacobian );
	CovMatrixDetNor.push_back( covMatrixFlightDirDet( 0 , 0 ) );
	CovMatrixDetNor.push_back( covMatrixFlightDirDet( 1 , 0 ) );
	CovMatrixDetNor.push_back( covMatrixFlightDirDet( 1 , 1 ) );
	CovMatrixDetNor.push_back( covMatrixFlightDirDet( 2 , 0 ) );
	CovMatrixDetNor.push_back( covMatrixFlightDirDet( 2 , 1 ) );
	CovMatrixDetNor.push_back( covMatrixFlightDirDet( 2 , 2 ) );
	CovMatrixDetNor.push_back( covMatrixFlightDirDet( 3 , 0 ) );
	CovMatrixDetNor.push_back( covMatrixFlightDirDet( 3 , 1 ) );
	CovMatrixDetNor.push_back( covMatrixFlightDirDet( 3 , 2 ) );
	CovMatrixDetNor.push_back( covMatrixFlightDirDet( 3 , 3 ) );
	streamlog_out(DEBUG9) << "	CovMatDetFlightNor :	" << CovMatrixDetNor[ 0 ] << std::endl;
	streamlog_out(DEBUG9) << "				" << CovMatrixDetNor[ 1 ] << "	,	" << CovMatrixDetNor[ 2 ] << std::endl;
	streamlog_out(DEBUG9) << "				" << CovMatrixDetNor[ 3 ] << "	,	" << CovMatrixDetNor[ 4 ] << "	,	" << CovMatrixDetNor[ 5 ] << std::endl;
	streamlog_out(DEBUG9) << "				" << CovMatrixDetNor[ 6 ] << "	,	" << CovMatrixDetNor[ 7 ] << "	,	" << CovMatrixDetNor[ 8 ] << "	,	" << CovMatrixDetNor[ 9 ] << std::endl;
	streamlog_out(DEBUG9) << "" << std::endl;
}

void SLDCorrection::getNeutrinoCovMat(		TLorentzVector recoNeutrinoFourMomentum , TLorentzVector visibleFourMomentum , TVector3 flightDirection ,
						double parentHadronMass , std::vector< float > CovMatrixPVA , std::vector< float > CovMatrixDetector ,
						std::vector< float > CovMatrixDetPar , std::vector< float > CovMatrixDetNor , std::vector< float > &NeutrinoCovMatrix )
{
	NeutrinoCovMatrix.clear();
	double E_nu = recoNeutrinoFourMomentum.E();
	double P_nu_par = ( recoNeutrinoFourMomentum.Vect() ).Dot( flightDirection );
	double E_vis = visibleFourMomentum.E();
	double P_vis_par = ( visibleFourMomentum.Vect() ).Dot( flightDirection );
	double P_vis_nor = sqrt( ( visibleFourMomentum.Vect() ).Mag2() - pow( P_vis_par , 2 ) );
//	double M_vis = visibleFourMomentum.M();
	streamlog_out(DEBUG9) << "	E_nu =		" << E_nu << std::endl;
	streamlog_out(DEBUG9) << "	P_nu_par =	" << P_nu_par << std::endl;
	streamlog_out(DEBUG9) << "	E_vis =		" << E_vis << std::endl;
	streamlog_out(DEBUG9) << "	P_vis_par =	" << P_vis_par << std::endl;
	streamlog_out(DEBUG9) << "	P_vis_nor =	" << P_vis_nor << std::endl;
	streamlog_out(DEBUG9) << "	M_B =		" << parentHadronMass << std::endl;
	double Coefficient = E_nu / ( E_nu * P_vis_par - E_vis * P_nu_par );
	double coefficient_Evis = Coefficient * ( E_nu + E_vis );
	double coefficient_PvisPar = Coefficient * ( P_nu_par + P_vis_par );
	double coefficient_PvisNor = Coefficient * P_vis_nor * ( 1.0 + E_vis / E_nu );
//	double coefficient_MB = Coefficient * parentHadronMass;
	float power = 1.0;

	NeutrinoCovMatrix.push_back( pow( fabs( coefficient_PvisPar ) , power / 2.0 ) * CovMatrixDetPar[ 0 ] + pow( fabs( coefficient_PvisNor ) , power / 2.0 ) * CovMatrixDetNor[ 0 ] );
	NeutrinoCovMatrix.push_back( pow( fabs( coefficient_PvisPar ) , power / 2.0 ) * CovMatrixDetPar[ 1 ] + pow( fabs( coefficient_PvisNor ) , power / 2.0 ) * CovMatrixDetNor[ 1 ] );
	NeutrinoCovMatrix.push_back( pow( fabs( coefficient_PvisPar ) , power / 2.0 ) * CovMatrixDetPar[ 2 ] + pow( fabs( coefficient_PvisNor ) , power / 2.0 ) * CovMatrixDetNor[ 2 ] );
	NeutrinoCovMatrix.push_back( pow( fabs( coefficient_PvisPar ) , power / 2.0 ) * CovMatrixDetPar[ 3 ] + pow( fabs( coefficient_PvisNor ) , power / 2.0 ) * CovMatrixDetNor[ 3 ] );
	NeutrinoCovMatrix.push_back( pow( fabs( coefficient_PvisPar ) , power / 2.0 ) * CovMatrixDetPar[ 4 ] + pow( fabs( coefficient_PvisNor ) , power / 2.0 ) * CovMatrixDetNor[ 4 ] );
	NeutrinoCovMatrix.push_back( pow( fabs( coefficient_PvisPar ) , power / 2.0 ) * CovMatrixDetPar[ 5 ] + pow( fabs( coefficient_PvisNor ) , power / 2.0 ) * CovMatrixDetNor[ 5 ] );
	NeutrinoCovMatrix.push_back( pow( fabs( coefficient_PvisPar ) , power / 2.0 ) * CovMatrixDetPar[ 6 ] + pow( fabs( coefficient_PvisNor ) , power / 2.0 ) * CovMatrixDetNor[ 6 ] + pow( fabs( coefficient_Evis ) , power / 2.0 ) * ( CovMatrixDetector[ 6 ] + CovMatrixPVA[ 6 ] ) );
	NeutrinoCovMatrix.push_back( pow( fabs( coefficient_PvisPar ) , power / 2.0 ) * CovMatrixDetPar[ 7 ] + pow( fabs( coefficient_PvisNor ) , power / 2.0 ) * CovMatrixDetNor[ 7 ] + pow( fabs( coefficient_Evis ) , power / 2.0 ) * ( CovMatrixDetector[ 7 ] + CovMatrixPVA[ 7 ] ) );
	NeutrinoCovMatrix.push_back( pow( fabs( coefficient_PvisPar ) , power / 2.0 ) * CovMatrixDetPar[ 8 ] + pow( fabs( coefficient_PvisNor ) , power / 2.0 ) * CovMatrixDetNor[ 8 ] + pow( fabs( coefficient_Evis ) , power / 2.0 ) * ( CovMatrixDetector[ 8 ] + CovMatrixPVA[ 8 ] ) );
	NeutrinoCovMatrix.push_back( pow( fabs( coefficient_Evis ) , power / 1.0 ) * ( CovMatrixDetector[ 9 ] + CovMatrixPVA[ 9 ] ) );

/*
	NeutrinoCovMatrix.push_back( pow( E_nu / ( E_vis * P_nu_par - E_nu * P_vis_par ) , 2 ) * ( pow( M_vis , 2 ) * CovMatrixPVA[ 0 ] + pow( P_vis_nor * ( 2 + E_vis / E_nu ) , 2 ) * CovMatrixDetNor[ 0 ] + pow( P_nu_par , 2 ) * CovMatrixDetPar[ 0 ] ) );
	NeutrinoCovMatrix.push_back( pow( E_nu / ( E_vis * P_nu_par - E_nu * P_vis_par ) , 2 ) * ( pow( M_vis , 2 ) * CovMatrixPVA[ 1 ] + pow( P_vis_nor * ( 2 + E_vis / E_nu ) , 2 ) * CovMatrixDetNor[ 1 ] + pow( P_nu_par , 2 ) * CovMatrixDetPar[ 1 ] ) );
	NeutrinoCovMatrix.push_back( pow( E_nu / ( E_vis * P_nu_par - E_nu * P_vis_par ) , 2 ) * ( pow( M_vis , 2 ) * CovMatrixPVA[ 2 ] + pow( P_vis_nor * ( 2 + E_vis / E_nu ) , 2 ) * CovMatrixDetNor[ 2 ] + pow( P_nu_par , 2 ) * CovMatrixDetPar[ 2 ] ) );
	NeutrinoCovMatrix.push_back( pow( E_nu / ( E_vis * P_nu_par - E_nu * P_vis_par ) , 2 ) * ( pow( M_vis , 2 ) * CovMatrixPVA[ 3 ] + pow( P_vis_nor * ( 2 + E_vis / E_nu ) , 2 ) * CovMatrixDetNor[ 3 ] + pow( P_nu_par , 2 ) * CovMatrixDetPar[ 3 ] ) );
	NeutrinoCovMatrix.push_back( pow( E_nu / ( E_vis * P_nu_par - E_nu * P_vis_par ) , 2 ) * ( pow( M_vis , 2 ) * CovMatrixPVA[ 4 ] + pow( P_vis_nor * ( 2 + E_vis / E_nu ) , 2 ) * CovMatrixDetNor[ 4 ] + pow( P_nu_par , 2 ) * CovMatrixDetPar[ 4 ] ) );
	NeutrinoCovMatrix.push_back( pow( E_nu / ( E_vis * P_nu_par - E_nu * P_vis_par ) , 2 ) * ( pow( M_vis , 2 ) * CovMatrixPVA[ 5 ] + pow( P_vis_nor * ( 2 + E_vis / E_nu ) , 2 ) * CovMatrixDetNor[ 5 ] + pow( P_nu_par , 2 ) * CovMatrixDetPar[ 5 ] ) );
	NeutrinoCovMatrix.push_back( pow( E_nu / ( E_vis * P_nu_par - E_nu * P_vis_par ) , 2 ) * ( pow( M_vis , 2 ) * CovMatrixPVA[ 6 ] + E_nu * CovMatrixDetector[ 6 ] ) );
	NeutrinoCovMatrix.push_back( pow( E_nu / ( E_vis * P_nu_par - E_nu * P_vis_par ) , 2 ) * ( pow( M_vis , 2 ) * CovMatrixPVA[ 7 ] + E_nu * CovMatrixDetector[ 7 ] ) );
	NeutrinoCovMatrix.push_back( pow( E_nu / ( E_vis * P_nu_par - E_nu * P_vis_par ) , 2 ) * ( pow( M_vis , 2 ) * CovMatrixPVA[ 8 ] + E_nu * CovMatrixDetector[ 8 ] ) );
	NeutrinoCovMatrix.push_back( pow( E_nu / ( E_vis * P_nu_par - E_nu * P_vis_par ) , 2 ) * ( pow( M_vis , 2 ) * CovMatrixPVA[ 9 ] + pow( E_nu , 2 ) * CovMatrixDetector[ 9 ] ) );
*/
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
	streamlog_out(DEBUG4) << "	TRUE PARENT HADRON MASS =  			" << parentHadron->getMass() << std::endl;
	streamlog_out(DEBUG4) << "	TRUE FLIHT DIRECTION (x,y,z): 			" << trueFliDir.X() << "	,	" << trueFliDir.Y() << "	,	" << trueFliDir.Z() << std::endl;
	streamlog_out(DEBUG4) << "	TRUE VISIBLE FOUR-MOMENTUM (Px,Py,Pz,E):	" << true4mom.Px() << "	,	" << true4mom.Py() << "	,	" << true4mom.Pz() << "	,	" << true4mom.E() << std::endl;
	TVector3 truePvisPar = trueFliDir.Dot( true4mom.Vect() ) * trueFliDir;
	TVector3 truePvisNor = true4mom.Vect() - truePvisPar;
	streamlog_out(DEBUG4) << "	TRUE VISIBLE FOUR-MOMENTUM(par) (Px,Py,Pz):	" << truePvisPar.Px() << "	,	" << truePvisPar.Py() << "	,	" << truePvisPar.Pz() << std::endl;
	streamlog_out(DEBUG4) << "	TRUE VISIBLE FOUR-MOMENTUM(nor) (Px,Py,Pz):	" << truePvisNor.Px() << "	,	" << truePvisNor.Py() << "	,	" << truePvisNor.Pz() << std::endl;

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
	streamlog_out(DEBUG4) << "		Visible Momentum (par-prime):		( " << visible_p_par_prime.Px() << "	, " << visible_p_par_prime.Py() << "	, " << visible_p_par_prime.Pz() << "	)" << std::endl;
	if ( pow( ( pow( ParentHadronMass , 2 ) - pow( visible_mass , 2 ) ) / ( 2 * ParentHadronMass ) , 2 ) < visible_p_nor.Mag2() )
	{
		visible_p_par_prime	= solutionSign * std::numeric_limits<double>::min() * flightDirection;
	}
	m_P_vis_par_prime.push_back( visible_p_par_prime.Mag() );
	streamlog_out(DEBUG4) << "		Visible Momentum (par-prime):		( " << visible_p_par_prime.Px() << "	, " << visible_p_par_prime.Py() << "	, " << visible_p_par_prime.Pz() << "	)" << std::endl;

	double parent_hadron_E		= ( ( visibleFourMomentum.E() * ( pow( ParentHadronMass , 2 ) + pow( visible_mass , 2 ) ) / ( 2 * ParentHadronMass ) ) - visible_p_par.Dot( visible_p_par_prime ) ) * ParentHadronMass / ( pow( visible_mass , 2 ) + visible_p_nor.Mag2() );
	streamlog_out(DEBUG4) << "		Parent Hadron Energy =									" << parent_hadron_E << std::endl;
	TVector3 parent_hadron_p	= sqrt( pow( parent_hadron_E , 2 ) - pow( ParentHadronMass , 2 ) ) * flightDirection;
	streamlog_out(DEBUG4) << "		Parent Hadron Momentum:			( " << parent_hadron_p.Px() << "	, " << parent_hadron_p.Py() << "	, " << parent_hadron_p.Pz() << "	, " << parent_hadron_E << " )" << std::endl;

	double sigma_E_vis = 0.0;
	double sigma_E_vis_prime = 0.0;
	double sigma_p_vis_par = 0.0;
	double sigma_p_vis_par_prime = 0.0;
	double sigma_p_vis_nor = 0.0;

	double sigma_parent_hadron_E2	=	pow( ParentHadronMass , 2 ) * (
						pow( visible_E_prime / ( pow( visible_mass , 2 ) + pow( visible_p_nor.Mag() , 2 ) ) , 2 ) * pow( sigma_E_vis , 2 ) +
						pow( visible_E / ( pow( visible_mass , 2 ) + pow( visible_p_nor.Mag() , 2 ) ) , 2 ) * pow( sigma_E_vis_prime , 2 ) +
						pow( visible_p_par_prime.Mag() / ( pow( visible_mass , 2 ) + pow( visible_p_nor.Mag() , 2 ) ) , 2 ) * pow( sigma_p_vis_par , 2 ) +
						pow( visible_p_par.Mag() / ( pow( visible_mass , 2 ) + pow( visible_p_nor.Mag() , 2 ) ) , 2 ) * pow( sigma_p_vis_par_prime , 2 ) +
						pow( 2 * visible_p_nor.Mag() * ( visible_E_prime * visible_E - visible_p_par_prime.Mag() * visible_p_par.Mag() ) / pow( ( pow( visible_mass , 2 ) + pow( visible_p_nor.Mag() , 2 ) ) , 2 ) , 2 ) * pow( sigma_p_vis_nor , 2 )
					);
					streamlog_out(DEBUG4) << "		Parent Hadron Sigma_E =									" << sigma_parent_hadron_E2 << std::endl;

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

TLorentzVector SLDCorrection::getNeutrinoFourMomentumModified( TVector3 flightDirection , TLorentzVector visibleFourMomentum , double ParentHadronMass , float solutionSign )
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
	streamlog_out(DEBUG8) << "		Visible Energy:									" << visible_E << std::endl;

	TVector3 visible_p		= TVector3( visibleFourMomentum.Px() , visibleFourMomentum.Py() , visibleFourMomentum.Pz() );
	streamlog_out(DEBUG4) << "		Visible Momentum:			( " << visible_p.Px() << "	, " << visible_p.Py() << "	, " << visible_p.Pz() << "	)" << std::endl;
	streamlog_out(DEBUG4) << "		(Visible Momentum).(FlightDirection):							" << visible_p.Dot( flightDirection ) << std::endl;

	double P_vis_par		= visible_p.Dot( flightDirection );
	TVector3 visible_p_par		= P_vis_par * flightDirection;
	m_P_vis_par.push_back( P_vis_par );
	streamlog_out(DEBUG8) << "		|Visible Momentum (par)|:	" << P_vis_par << std::endl;
	streamlog_out(DEBUG4) << "		Visible Momentum (par):			( " << visible_p_par.Px() << "	, " << visible_p_par.Py() << "	, " << visible_p_par.Pz() << "	)" << std::endl;

	TVector3 visible_p_nor		= visible_p - visible_p_par;
	double P_vis_nor		= visible_p_nor.Mag();
	m_P_vis_nor.push_back( P_vis_nor );
	m_P_vis_nor_prime.push_back( P_vis_nor );
	streamlog_out(DEBUG8) << "		|Visible Momentum (nor)|:	" << P_vis_nor << std::endl;
	streamlog_out(DEBUG4) << "		Visible Momentum (nor):			( " << visible_p_nor.Px() << "	, " << visible_p_nor.Py() << "	, " << visible_p_nor.Pz() << "	)" << std::endl;

	double visible_E_prime		= ( pow( ParentHadronMass , 2 ) + pow( visible_mass , 2 ) ) / ( 2 * ParentHadronMass );
	m_E_vis_prime.push_back( visible_E_prime );
	streamlog_out(DEBUG4) << "		Visible Energy (prime):								" << visible_E_prime << std::endl;

	TVector3 visible_p_par_prime	= solutionSign * sqrt( pow( ( pow( ParentHadronMass , 2 ) - pow( visible_mass , 2 ) ) , 2 ) / 4 - pow( ParentHadronMass * P_vis_nor , 2 ) ) / ParentHadronMass * flightDirection;
	double P_vis_par_prime	= solutionSign * sqrt( pow( ( pow( ParentHadronMass , 2 ) - pow( visible_mass , 2 ) ) , 2 ) / 4 - pow( ParentHadronMass * P_vis_nor , 2 ) ) / ParentHadronMass;
	m_P_vis_par_prime.push_back( P_vis_par_prime );
	streamlog_out(DEBUG4) << "		Visible Momentum (par-prime):		( " << visible_p_par_prime.Px() << "	, " << visible_p_par_prime.Py() << "	, " << visible_p_par_prime.Pz() << "	)" << std::endl;
	if ( pow( ( pow( ParentHadronMass , 2 ) - pow( visible_mass , 2 ) ) , 2 ) < 4.0 * pow( P_vis_nor * ParentHadronMass , 2 ) )
	{
		visible_p_par_prime	= solutionSign * std::numeric_limits<double>::min() * flightDirection;
		P_vis_par_prime		= 0.0;
	}
	streamlog_out(DEBUG4) << "		Visible Momentum (par-prime):		( " << visible_p_par_prime.Px() << "	, " << visible_p_par_prime.Py() << "	, " << visible_p_par_prime.Pz() << "	)" << std::endl;
	streamlog_out(DEBUG8) << "		|Visible Momentum (par-prime)|:	" << P_vis_par_prime << std::endl;

	double parent_hadron_E		= ( ( visibleFourMomentum.E() * ( pow( ParentHadronMass , 2 ) + pow( visible_mass , 2 ) ) / ( 2 * ParentHadronMass ) ) - visible_p_par.Dot( visible_p_par_prime ) ) * ParentHadronMass / ( pow( visible_mass , 2 ) + visible_p_nor.Mag2() );
	streamlog_out(DEBUG4) << "		Parent Hadron Energy =									" << parent_hadron_E << std::endl;
	TVector3 parent_hadron_p	= ( ParentHadronMass >= parent_hadron_E ? 0.0 : sqrt( pow( parent_hadron_E , 2 ) - pow( ParentHadronMass , 2 ) ) ) * flightDirection;
	streamlog_out(DEBUG4) << "		Parent Hadron Momentum:			( " << parent_hadron_p.Px() << "	, " << parent_hadron_p.Py() << "	, " << parent_hadron_p.Pz() << "	, " << parent_hadron_E << " )" << std::endl;

	TVector3 Neutrino_p_nor		= -1 * visible_p_nor;
	streamlog_out(DEBUG8) << "		Neutrino Momentum (nor):		( " << Neutrino_p_nor.Px() << "	, " << Neutrino_p_nor.Py() << "	, " << Neutrino_p_nor.Pz() << "	)" << std::endl;

	TVector3 Neutrino_p_par		= parent_hadron_p - visible_p_par;
	streamlog_out(DEBUG8) << "		Neutrino Momentum (par):		( " << Neutrino_p_par.Px() << "	, " << Neutrino_p_par.Py() << "	, " << Neutrino_p_par.Pz() << "	)" << std::endl;

	TVector3 Neutrino_p		= Neutrino_p_nor + Neutrino_p_par;
	streamlog_out(DEBUG8) << "		Neutrino Momentum:			( " << Neutrino_p.Px() << "	, " << Neutrino_p.Py() << "	, " << Neutrino_p.Pz() << "	)" << std::endl;

	double Neutrino_E		= Neutrino_p.Mag();
	streamlog_out(DEBUG8) << "		Neutrino Energy:									" << Neutrino_E << std::endl;

	TLorentzVector Neutrino_tlv( Neutrino_p , Neutrino_E );
	streamlog_out(DEBUG4) << "" << std::endl;
	streamlog_out(DEBUG4) << "	Reconstructed Neutrino Mass:	" << Neutrino_tlv.Mag() << std::endl;
	streamlog_out(DEBUG4) << "	Reconstructed Neutrino 4-Momentum:		( " << Neutrino_tlv.Px() << "	, " << Neutrino_tlv.Py() << "	, " << Neutrino_tlv.Pz() << "	, " << Neutrino_tlv.E() << " )" << std::endl;
	return Neutrino_tlv;
}

TLorentzVector SLDCorrection::getNeutrinoFourMomentumStandardMethod( TVector3 flightDirection , TLorentzVector visibleFourMomentum , double parentHadronMass , float solutionSign )
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
	streamlog_out(DEBUG4) << "		Parent Hadron Mass =	 " << parentHadronMass << std::endl;

	streamlog_out(DEBUG4) << "		Visible 4-Momentum:			( " << visibleFourMomentum.Px() << "	, " << visibleFourMomentum.Py() << "	, " << visibleFourMomentum.Pz() << "	, " << visibleFourMomentum.E() << " )" << std::endl;

	double visible_mass		= visibleFourMomentum.M();
	streamlog_out(DEBUG4) << "		Visible Inv Mass:	" << visible_mass << std::endl;

	double visible_E		= visibleFourMomentum.E();
	m_E_vis.push_back( visible_E );
	streamlog_out(DEBUG8) << "		Visible Energy:									" << visible_E << std::endl;

	TVector3 visible_p		= TVector3( visibleFourMomentum.Px() , visibleFourMomentum.Py() , visibleFourMomentum.Pz() );
	streamlog_out(DEBUG4) << "		Visible Momentum:			( " << visible_p.Px() << "	, " << visible_p.Py() << "	, " << visible_p.Pz() << "	)" << std::endl;
	streamlog_out(DEBUG4) << "		(Visible Momentum).(FlightDirection):							" << visible_p.Dot( flightDirection ) << std::endl;

	double P_vis_par		= visible_p.Dot( flightDirection );
	TVector3 visible_p_par		= P_vis_par * flightDirection;
	m_P_vis_par.push_back( P_vis_par );
	streamlog_out(DEBUG8) << "		|Visible Momentum (par)|:	" << P_vis_par << std::endl;
	streamlog_out(DEBUG4) << "		Visible Momentum (par):			( " << visible_p_par.Px() << "	, " << visible_p_par.Py() << "	, " << visible_p_par.Pz() << "	)" << std::endl;

	double P_vis_nor		= sqrt( visible_p.Mag2() - pow( P_vis_par , 2 ) );
	TVector3 visible_p_nor		= visible_p - visible_p_par;
	m_P_vis_nor.push_back( P_vis_nor );
	m_P_vis_nor_prime.push_back( P_vis_nor );
	streamlog_out(DEBUG8) << "		|Visible Momentum (nor)|:	" << P_vis_nor << std::endl;
	streamlog_out(DEBUG4) << "		Visible Momentum (nor):			( " << visible_p_nor.Px() << "	, " << visible_p_nor.Py() << "	, " << visible_p_nor.Pz() << "	)" << std::endl;
/*
	double A = P_vis_par * ( pow( parentHadronMass , 2 ) - pow( visible_mass , 2 ) - 2.0 * pow( P_vis_nor , 2 ) );
	streamlog_out(DEBUG8) << "		(A) ( M_B^2 - M_vis^2 - 2 * P_vis(nor)^2 ) * P_vis(par) =	" << A << std::endl;
	double B = 4.0 * pow( P_vis_nor , 2 ) * pow( visible_E , 2 ) - pow( pow( parentHadronMass , 2 ) - pow( visible_mass , 2 ) - 2.0 * pow( P_vis_nor , 2 ) , 2 );
	streamlog_out(DEBUG8) << "		(B) 4 * E_vis^2 * P_vis(par)^2 - ( M_B^2 - M_vis^2 - 2 * P_vis(nor)^2 )^2 =	" << B << std::endl;
*/
	double Q = pow( parentHadronMass , 2 ) - pow( visible_mass , 2 ) - 2.0 * pow( P_vis_nor , 2 );
	streamlog_out(DEBUG8) << "		(Q) M_B^2 - M_vis^2 - 2 * P_vis(nor)^2 =	" << Q << std::endl;
	double A = pow( visible_E , 2 ) - pow( P_vis_par , 2 );
	streamlog_out(DEBUG8) << "		(A) E_vis^2 - P_vis(par)^2 =	" << A << std::endl;
	double B = Q * P_vis_par;
	streamlog_out(DEBUG8) << "		(B) ( M_B^2 - M_vis^2 - 2 * P_vis(nor)^2 ) * P_vis(par) =	" << B << std::endl;
	double C = pow( P_vis_nor , 2 ) * pow( visible_E , 2 ) - pow( Q , 2 ) / 4.0;
	streamlog_out(DEBUG8) << "		(C) E_vis^2 * P_vis(nor)^2 - ( M_B^2 - M_vis^2 - 2 * P_vis(nor)^2 )^2 / 4 =	" << C << std::endl;

	double P_nu_par = ( pow( B , 2 ) >= 4.0 * A * C ? ( B + solutionSign * std::sqrt( pow( B , 2 ) - 4.0 * A * C ) ) / ( 2.0 * A ) : B / ( 2.0 * A ) );
	streamlog_out(DEBUG8) << "		|Neutrino Momentum (par)|=		" << P_nu_par << std::endl;
	TVector3 Neutrino_p_par		= P_nu_par * flightDirection;
	streamlog_out(DEBUG8) << "		Neutrino Momentum (par):		( " << Neutrino_p_par.Px() << "	, " << Neutrino_p_par.Py() << "	, " << Neutrino_p_par.Pz() << "	)" << std::endl;

	TVector3 Neutrino_p_nor		= -1.0 * visible_p_nor;
	streamlog_out(DEBUG8) << "		Neutrino Momentum (nor):		( " << Neutrino_p_nor.Px() << "	, " << Neutrino_p_nor.Py() << "	, " << Neutrino_p_nor.Pz() << "	)" << std::endl;
	TVector3 Neutrino_p		= Neutrino_p_nor + Neutrino_p_par;
	streamlog_out(DEBUG8) << "		Neutrino Momentum:			( " << Neutrino_p.Px() << "	, " << Neutrino_p.Py() << "	, " << Neutrino_p.Pz() << "	)" << std::endl;
	double Neutrino_E		= Neutrino_p.Mag();
	streamlog_out(DEBUG8) << "		Neutrino Energy:									" << Neutrino_E << std::endl;

	TLorentzVector Neutrino_tlv( Neutrino_p , Neutrino_E );
	streamlog_out(DEBUG4) << "" << std::endl;
	streamlog_out(DEBUG4) << "	Reconstructed Neutrino Mass:	" << Neutrino_tlv.Mag() << std::endl;
	streamlog_out(DEBUG4) << "	Reconstructed Neutrino 4-Momentum:		( " << Neutrino_tlv.Px() << "	, " << Neutrino_tlv.Py() << "	, " << Neutrino_tlv.Pz() << "	, " << Neutrino_tlv.E() << " )" << std::endl;
	return Neutrino_tlv;

}

MCParticle* SLDCorrection::getTrueNeutrino( MCParticle *SLDLepton )
{
	MCParticle* trueNeutrino{};
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

void SLDCorrection::evaluatePFOsAngle(	std::vector<EVENT::ReconstructedParticle*> aloneChargedPFOsInJetFromSLD ,
					std::vector<EVENT::ReconstructedParticle*> aloneChargedPFOsInJetNotFromSLD ,
					std::vector<EVENT::ReconstructedParticle*> chargedPFOsInJetFromSLD ,
					std::vector<EVENT::ReconstructedParticle*> chargedPFOsInJetNotFromSLD ,
					std::vector<EVENT::ReconstructedParticle*> neutralPFOsInJetFromSLD ,
					std::vector<EVENT::ReconstructedParticle*> neutralPFOsInJetNotFromSLD ,
					TVector3 leptonDirection , TVector3 jetAxis , TVector3 recoFlightDirection , int SLDStatus )
{
	leptonDirection.SetMag( 1.0 );
	jetAxis.SetMag( 1.0 );
	if ( SLDStatus == 4 || SLDStatus == 5 )recoFlightDirection.SetMag( 1.0 );
	for ( unsigned int i_par = 0 ; i_par < aloneChargedPFOsInJetFromSLD.size() ; ++i_par )
	{
		TVector3 pfoMomentum( ( aloneChargedPFOsInJetFromSLD[ i_par ] )->getMomentum() ); pfoMomentum.SetMag( 1.0 );
		m_alphaAloneChargedPFOsFromSLDwrtLepton.push_back( acos( pfoMomentum.Dot( leptonDirection ) ) );
		m_cosAlphaAloneChargedPFOsFromSLDwrtLepton.push_back( pfoMomentum.Dot( leptonDirection ) );
		m_alphaAloneChargedPFOsFromSLDwrtJet.push_back( acos( pfoMomentum.Dot( jetAxis ) ) );
		m_cosAlphaAloneChargedPFOsFromSLDwrtJet.push_back( pfoMomentum.Dot( jetAxis ) );
		if ( SLDStatus == 4 || SLDStatus == 5 )
		{
			m_alphaAloneChargedPFOsFromSLDwrtFD.push_back( acos( pfoMomentum.Dot( recoFlightDirection ) ) );
			m_cosAlphaAloneChargedPFOsFromSLDwrtFD.push_back( pfoMomentum.Dot( recoFlightDirection ) );
		}

	}
	for ( unsigned int i_par = 0 ; i_par < aloneChargedPFOsInJetNotFromSLD.size() ; ++i_par )
	{
		TVector3 pfoMomentum( ( aloneChargedPFOsInJetNotFromSLD[ i_par ] )->getMomentum() ); pfoMomentum.SetMag( 1.0 );
		m_alphaAloneChargedPFOsNotFromSLDwrtLepton.push_back( acos( pfoMomentum.Dot( leptonDirection ) ) );
		m_cosAlphaAloneChargedPFOsNotFromSLDwrtLepton.push_back( pfoMomentum.Dot( leptonDirection ) );
		m_alphaAloneChargedPFOsNotFromSLDwrtJet.push_back( acos( pfoMomentum.Dot( jetAxis ) ) );
		m_cosAlphaAloneChargedPFOsNotFromSLDwrtJet.push_back( pfoMomentum.Dot( jetAxis ) );
		if ( SLDStatus == 4 || SLDStatus == 5 )
		{
			m_alphaAloneChargedPFOsNotFromSLDwrtFD.push_back( acos( pfoMomentum.Dot( recoFlightDirection ) ) );
			m_cosAlphaAloneChargedPFOsNotFromSLDwrtFD.push_back( pfoMomentum.Dot( recoFlightDirection ) );
		}

	}
	for ( unsigned int i_par = 0 ; i_par < chargedPFOsInJetFromSLD.size() ; ++i_par )
	{
		TVector3 pfoMomentum( ( chargedPFOsInJetFromSLD[ i_par ] )->getMomentum() ); pfoMomentum.SetMag( 1.0 );
		m_alphaChargedPFOsFromSLDwrtLepton.push_back( acos( pfoMomentum.Dot( leptonDirection ) ) );
		m_cosAlphaChargedPFOsFromSLDwrtLepton.push_back( pfoMomentum.Dot( leptonDirection ) );
		m_alphaChargedPFOsFromSLDwrtJet.push_back( acos( pfoMomentum.Dot( jetAxis ) ) );
		m_cosAlphaChargedPFOsFromSLDwrtJet.push_back( pfoMomentum.Dot( jetAxis ) );
		if ( SLDStatus == 4 || SLDStatus == 5 )
		{
			m_alphaChargedPFOsFromSLDwrtFD.push_back( acos( pfoMomentum.Dot( recoFlightDirection ) ) );
			m_cosAlphaChargedPFOsFromSLDwrtFD.push_back( pfoMomentum.Dot( recoFlightDirection ) );
		}

	}
	for ( unsigned int i_par = 0 ; i_par < chargedPFOsInJetNotFromSLD.size() ; ++i_par )
	{
		TVector3 pfoMomentum( ( chargedPFOsInJetNotFromSLD[ i_par ] )->getMomentum() ); pfoMomentum.SetMag( 1.0 );
		m_alphaChargedPFOsNotFromSLDwrtLepton.push_back( acos( pfoMomentum.Dot( leptonDirection ) ) );
		m_cosAlphaChargedPFOsNotFromSLDwrtLepton.push_back( pfoMomentum.Dot( leptonDirection ) );
		m_alphaChargedPFOsNotFromSLDwrtJet.push_back( acos( pfoMomentum.Dot( jetAxis ) ) );
		m_cosAlphaChargedPFOsNotFromSLDwrtJet.push_back( pfoMomentum.Dot( jetAxis ) );
		if ( SLDStatus == 4 || SLDStatus == 5 )
		{
			m_alphaChargedPFOsNotFromSLDwrtFD.push_back( acos( pfoMomentum.Dot( recoFlightDirection ) ) );
			m_cosAlphaChargedPFOsNotFromSLDwrtFD.push_back( pfoMomentum.Dot( recoFlightDirection ) );
		}

	}
	for ( unsigned int i_par = 0 ; i_par < neutralPFOsInJetFromSLD.size() ; ++i_par )
	{
		TVector3 pfoMomentum( ( neutralPFOsInJetFromSLD[ i_par ] )->getMomentum() ); pfoMomentum.SetMag( 1.0 );
		m_alphaNeutralPFOsFromSLDwrtLepton.push_back( acos( pfoMomentum.Dot( leptonDirection ) ) );
		m_cosAlphaNeutralPFOsFromSLDwrtLepton.push_back( pfoMomentum.Dot( leptonDirection ) );
		m_alphaNeutralPFOsFromSLDwrtJet.push_back( acos( pfoMomentum.Dot( jetAxis ) ) );
		m_cosAlphaNeutralPFOsFromSLDwrtJet.push_back( pfoMomentum.Dot( jetAxis ) );
		if ( SLDStatus == 4 || SLDStatus == 5 )
		{
			m_alphaNeutralPFOsFromSLDwrtFD.push_back( acos( pfoMomentum.Dot( recoFlightDirection ) ) );
			m_cosAlphaNeutralPFOsFromSLDwrtFD.push_back( pfoMomentum.Dot( recoFlightDirection ) );
		}

	}
	for ( unsigned int i_par = 0 ; i_par < neutralPFOsInJetNotFromSLD.size() ; ++i_par )
	{
		TVector3 pfoMomentum( ( neutralPFOsInJetNotFromSLD[ i_par ] )->getMomentum() ); pfoMomentum.SetMag( 1.0 );
		m_alphaNeutralPFOsNotFromSLDwrtLepton.push_back( acos( pfoMomentum.Dot( leptonDirection ) ) );
		m_cosAlphaNeutralPFOsNotFromSLDwrtLepton.push_back( pfoMomentum.Dot( leptonDirection ) );
		m_alphaNeutralPFOsNotFromSLDwrtJet.push_back( acos( pfoMomentum.Dot( jetAxis ) ) );
		m_cosAlphaNeutralPFOsNotFromSLDwrtJet.push_back( pfoMomentum.Dot( jetAxis ) );
		if ( SLDStatus == 4 || SLDStatus == 5 )
		{
			m_alphaNeutralPFOsNotFromSLDwrtFD.push_back( acos( pfoMomentum.Dot( recoFlightDirection ) ) );
			m_cosAlphaNeutralPFOsNotFromSLDwrtFD.push_back( pfoMomentum.Dot( recoFlightDirection ) );
		}

	}

}

void SLDCorrection::fillTrueRecoFourMomentum(	TLorentzVector trueNeutralFourMomentum , TLorentzVector trueChargedFourMomentum ,
						TLorentzVector trueLeptonFourMomentum , TLorentzVector trueVisibleFourMomentum ,
						TLorentzVector trueNeutrinoFourMomentum , TLorentzVector trueHadronFourMomentum ,
						TLorentzVector recoNeutralFourMomentum , TLorentzVector recoChargedFourMomentum ,
						TLorentzVector recoLeptonFourMomentum , TLorentzVector recoVisibleFourMomentum ,
						TLorentzVector recoNeutrinoFourMomentum , TLorentzVector recoHadronFourMomentum ,
						TLorentzVector usedNeutralFourMomentum , TLorentzVector usedChargedFourMomentum ,
						TLorentzVector usedLeptonFourMomentum , TLorentzVector usedVisibleFourMomentum ,
					 	TLorentzVector cheatedPVARecoChargedFourMomentum , TLorentzVector cheatedPVARecoNeutralFourMomentum )
{
	m_trueNeutralPx.push_back( trueNeutralFourMomentum.Px() );
	m_trueNeutralPy.push_back( trueNeutralFourMomentum.Py() );
	m_trueNeutralPz.push_back( trueNeutralFourMomentum.Pz() );
	m_trueNeutralE.push_back( trueNeutralFourMomentum.E() );
	m_trueNeutralM.push_back( trueNeutralFourMomentum.M() );
	m_trueChargedPx.push_back( trueChargedFourMomentum.Px() );
	m_trueChargedPy.push_back( trueChargedFourMomentum.Py() );
	m_trueChargedPz.push_back( trueChargedFourMomentum.Pz() );
	m_trueChargedE.push_back( trueChargedFourMomentum.E() );
	m_trueChargedM.push_back( trueChargedFourMomentum.M() );
	m_trueLeptonPx.push_back( trueLeptonFourMomentum.Px() );
	m_trueLeptonPy.push_back( trueLeptonFourMomentum.Py() );
	m_trueLeptonPz.push_back( trueLeptonFourMomentum.Pz() );
	m_trueLeptonE.push_back( trueLeptonFourMomentum.E() );
	m_trueVisiblePx.push_back( trueVisibleFourMomentum.Px() );
	m_trueVisiblePy.push_back( trueVisibleFourMomentum.Py() );
	m_trueVisiblePz.push_back( trueVisibleFourMomentum.Pz() );
	m_trueVisibleE.push_back( trueVisibleFourMomentum.E() );
	m_trueVisibleM.push_back( trueVisibleFourMomentum.M() );
	m_trueNeutrinoPx.push_back( trueNeutrinoFourMomentum.Px() );
	m_trueNeutrinoPy.push_back( trueNeutrinoFourMomentum.Py() );
	m_trueNeutrinoPz.push_back( trueNeutrinoFourMomentum.Pz() );
	m_trueNeutrinoE.push_back( trueNeutrinoFourMomentum.E() );
	m_trueHadronPx.push_back( trueHadronFourMomentum.Px() );
	m_trueHadronPy.push_back( trueHadronFourMomentum.Py() );
	m_trueHadronPz.push_back( trueHadronFourMomentum.Pz() );
	m_trueHadronE.push_back( trueHadronFourMomentum.E() );
	m_trueHadronM.push_back( trueHadronFourMomentum.M() );
	m_recoNeutralPx.push_back( recoNeutralFourMomentum.Px() );
	m_recoNeutralPy.push_back( recoNeutralFourMomentum.Py() );
	m_recoNeutralPz.push_back( recoNeutralFourMomentum.Pz() );
	m_recoNeutralE.push_back( recoNeutralFourMomentum.E() );
	m_recoNeutralM.push_back( recoNeutralFourMomentum.M() );
	m_recoChargedPx.push_back( recoChargedFourMomentum.Px() );
	m_recoChargedPy.push_back( recoChargedFourMomentum.Py() );
	m_recoChargedPz.push_back( recoChargedFourMomentum.Pz() );
	m_recoChargedE.push_back( recoChargedFourMomentum.E() );
	m_recoChargedM.push_back( recoChargedFourMomentum.M() );
	m_recoLeptonPx.push_back( recoLeptonFourMomentum.Px() );
	m_recoLeptonPy.push_back( recoLeptonFourMomentum.Py() );
	m_recoLeptonPz.push_back( recoLeptonFourMomentum.Pz() );
	m_recoLeptonE.push_back( recoLeptonFourMomentum.E() );
	m_recoVisiblePx.push_back( recoVisibleFourMomentum.Px() );
	m_recoVisiblePy.push_back( recoVisibleFourMomentum.Py() );
	m_recoVisiblePz.push_back( recoVisibleFourMomentum.Pz() );
	m_recoVisibleE.push_back( recoVisibleFourMomentum.E() );
	m_recoVisibleM.push_back( recoVisibleFourMomentum.M() );
	m_recoNeutrinoPx.push_back( recoNeutrinoFourMomentum.Px() );
	m_recoNeutrinoPy.push_back( recoNeutrinoFourMomentum.Py() );
	m_recoNeutrinoPz.push_back( recoNeutrinoFourMomentum.Pz() );
	m_recoNeutrinoE.push_back( recoNeutrinoFourMomentum.E() );
	m_recoHadronPx.push_back( recoHadronFourMomentum.Px() );
	m_recoHadronPy.push_back( recoHadronFourMomentum.Py() );
	m_recoHadronPz.push_back( recoHadronFourMomentum.Pz() );
	m_recoHadronE.push_back( recoHadronFourMomentum.E() );
	m_recoHadronM.push_back( recoHadronFourMomentum.M() );
	m_usedNeutralPx.push_back( usedNeutralFourMomentum.Px() );
	m_usedNeutralPy.push_back( usedNeutralFourMomentum.Py() );
	m_usedNeutralPz.push_back( usedNeutralFourMomentum.Pz() );
	m_usedNeutralE.push_back( usedNeutralFourMomentum.E() );
	m_usedNeutralM.push_back( usedNeutralFourMomentum.M() );
	m_usedChargedPx.push_back( usedChargedFourMomentum.Px() );
	m_usedChargedPy.push_back( usedChargedFourMomentum.Py() );
	m_usedChargedPz.push_back( usedChargedFourMomentum.Pz() );
	m_usedChargedE.push_back( usedChargedFourMomentum.E() );
	m_usedChargedM.push_back( usedChargedFourMomentum.M() );
	m_usedLeptonPx.push_back( usedLeptonFourMomentum.Px() );
	m_usedLeptonPy.push_back( usedLeptonFourMomentum.Py() );
	m_usedLeptonPz.push_back( usedLeptonFourMomentum.Pz() );
	m_usedLeptonE.push_back( usedLeptonFourMomentum.E() );
	m_usedVisiblePx.push_back( usedVisibleFourMomentum.Px() );
	m_usedVisiblePy.push_back( usedVisibleFourMomentum.Py() );
	m_usedVisiblePz.push_back( usedVisibleFourMomentum.Pz() );
	m_usedVisibleE.push_back( usedVisibleFourMomentum.E() );
	m_usedVisibleM.push_back( usedVisibleFourMomentum.M() );
	m_cheatedPVARecoNeutralPx.push_back( cheatedPVARecoNeutralFourMomentum.Px() );
	m_cheatedPVARecoNeutralPy.push_back( cheatedPVARecoNeutralFourMomentum.Py() );
	m_cheatedPVARecoNeutralPz.push_back( cheatedPVARecoNeutralFourMomentum.Pz() );
	m_cheatedPVARecoNeutralE.push_back( cheatedPVARecoNeutralFourMomentum.E() );
	m_cheatedPVARecoNeutralM.push_back( cheatedPVARecoNeutralFourMomentum.M() );
	m_cheatedPVARecoChargedPx.push_back( cheatedPVARecoChargedFourMomentum.Px() );
	m_cheatedPVARecoChargedPy.push_back( cheatedPVARecoChargedFourMomentum.Py() );
	m_cheatedPVARecoChargedPz.push_back( cheatedPVARecoChargedFourMomentum.Pz() );
	m_cheatedPVARecoChargedE.push_back( cheatedPVARecoChargedFourMomentum.E() );
	m_cheatedPVARecoChargedM.push_back( cheatedPVARecoChargedFourMomentum.M() );
	m_leptonE_to_parentE.push_back( trueLeptonFourMomentum.E() / trueHadronFourMomentum.E() );
	m_otherChargedE_to_parentE.push_back( trueChargedFourMomentum.E() / trueHadronFourMomentum.E() );
	m_allChargedE_to_parentE.push_back( ( trueLeptonFourMomentum + trueChargedFourMomentum ).E() / trueHadronFourMomentum.E() );
	m_neutralE_to_parentE.push_back( trueNeutralFourMomentum.E() / trueHadronFourMomentum.E() );
	m_neutrino_to_parentE.push_back( trueNeutrinoFourMomentum.E() / trueHadronFourMomentum.E() );
}

void SLDCorrection::investigateJetEnergyContent( EVENT::ReconstructedParticle *assignedJet )
{
	double jetEnergy = assignedJet->getEnergy();
	double neutralHadronEnergy = 0.0;
	double chargedEnergy = 0.0;
	double photonEnergy = 0.0;
	double neutralsEnergy = 0.0;
	for ( unsigned int i_par = 0 ; i_par < assignedJet->getParticles().size() ; ++i_par )
	{
		if ( ( assignedJet->getParticles()[ i_par ] )->getTracks().size() != 0 )
		{
			chargedEnergy += assignedJet->getParticles()[ i_par ]->getEnergy();
		}
		else if ( ( assignedJet->getParticles()[ i_par ] )->getType() == 22 )
		{
			photonEnergy += assignedJet->getParticles()[ i_par ]->getEnergy();
			neutralsEnergy += assignedJet->getParticles()[ i_par ]->getEnergy();
		}
		else
		{
			neutralHadronEnergy += assignedJet->getParticles()[ i_par ]->getEnergy();
			neutralsEnergy += assignedJet->getParticles()[ i_par ]->getEnergy();
		}
	}
	m_jetEnergyFractionCharged.push_back( chargedEnergy / jetEnergy );
	m_jetEnergyFractionNeutralHadron.push_back( neutralHadronEnergy / jetEnergy );
	m_jetEnergyFractionPhoton.push_back( photonEnergy / jetEnergy );
	m_jetEnergyFractionNeutrals.push_back( neutralsEnergy / jetEnergy );
	m_jetEnergy.push_back( jetEnergy );
}

void SLDCorrection::plotHistograms( TLorentzVector trueFourMomentumNeutrino , TLorentzVector FourMomentumNuClose , std::vector<float> NeutrinoCovMat )
{
	double NuPxResidual = FourMomentumNuClose.Px() - trueFourMomentumNeutrino.Px(); m_NuPxResidual.push_back( NuPxResidual );
	double NuPyResidual = FourMomentumNuClose.Py() - trueFourMomentumNeutrino.Py(); m_NuPyResidual.push_back( NuPyResidual );
	double NuPzResidual = FourMomentumNuClose.Pz() - trueFourMomentumNeutrino.Pz(); m_NuPzResidual.push_back( NuPzResidual );
	double NuEResidual = FourMomentumNuClose.E() - trueFourMomentumNeutrino.E(); m_NuEResidual.push_back( NuEResidual );
	double NuPxNormalizedResidual = NuPxResidual / sqrt( NeutrinoCovMat[ 0 ] ); m_NuPxNormalizedResidual.push_back( NuPxNormalizedResidual );
	double NuPyNormalizedResidual = NuPyResidual / sqrt( NeutrinoCovMat[ 2 ] ); m_NuPyNormalizedResidual.push_back( NuPyNormalizedResidual );
	double NuPzNormalizedResidual = NuPzResidual / sqrt( NeutrinoCovMat[ 5 ] ); m_NuPzNormalizedResidual.push_back( NuPzNormalizedResidual );
	double NuENormalizedResidual = NuEResidual / sqrt( NeutrinoCovMat[ 9 ] ); m_NuENormalizedResidual.push_back( NuENormalizedResidual );
	h_NuPxResidual->Fill( NuPxResidual ); ++n_NuPxResidual;
	h_NuPyResidual->Fill( NuPyResidual ); ++n_NuPyResidual;
	h_NuPzResidual->Fill( NuPzResidual ); ++n_NuPzResidual;
	h_NuEResidual->Fill( NuEResidual ); ++n_NuEResidual;
	h_NuPxNormalizedResidual->Fill( NuPxNormalizedResidual ); ++n_NuPxNormalizedResidual;
	h_NuPyNormalizedResidual->Fill( NuPyNormalizedResidual ); ++n_NuPyNormalizedResidual;
	h_NuPzNormalizedResidual->Fill( NuPzNormalizedResidual ); ++n_NuPzNormalizedResidual;
	h_NuENormalizedResidual->Fill( NuENormalizedResidual ); ++n_NuENormalizedResidual;
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
//	float fit_range = 4.0;
//	float fit_min = -2.0;
//	float fit_max = 2.0;
//	doProperGaussianFit( histogram , fit_min , fit_max , fit_range );
//	histogram->GetFunction("gaus")->SetLineColor( color );
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
	try
	{
		streamlog_out(DEBUG4) << "	CHECK : processed events: " << m_nEvtSum << " , Number of semi-leptonic vertex: " << pLCEvent->getCollection(m_SLDVertex)->getNumberOfElements() << std::endl;
		streamlog_out(DEBUG4) << "	CHECK : processed events: " << m_nEvtSum << " , Number of semi-leptonic vertexRP: " << pLCEvent->getCollection(m_SLDVertexRP)->getNumberOfElements() << std::endl;
		streamlog_out(DEBUG4) << "	CHECK : processed events: " << m_nEvtSum << " , Number of Reconstrcuted Neutrinos: " << pLCEvent->getCollection(m_reconstructedNeutrino)->getNumberOfElements() << std::endl;
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
//		h_SLDStatus->Scale( 100.0 / n_SLDStatus );
		h_SLDStatus->GetYaxis()->SetTitle("#SLDecay [%]");
		h_SLDStatus->Write();
		h_BHadronType->Write();
		h_CHadronType->Write();
		InitializeHistogram( h_NuPxResidual , n_NuPxResidual , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_NuPyResidual , n_NuPyResidual , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_NuPzResidual , n_NuPzResidual , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_NuEResidual , n_NuEResidual , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_NuPxNormalizedResidual , n_NuPxNormalizedResidual , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_NuPyNormalizedResidual , n_NuPyNormalizedResidual , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_NuPzNormalizedResidual , n_NuPzNormalizedResidual , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_NuENormalizedResidual , n_NuENormalizedResidual , 4 , 1 , 1.0 , 1 );
		h_SLDecayFlavour->Write();
		h_SLDecayModeB->Write();
		h_SLDecayModeC->Write();
		h_SLDecayOrder->Write();
		h_NuPxResidual->Write();
		h_NuPyResidual->Write();
		h_NuPzResidual->Write();
		h_NuEResidual->Write();
		h_NuPxNormalizedResidual->Write();
		h_NuPyNormalizedResidual->Write();
		h_NuPzNormalizedResidual->Write();
		h_NuENormalizedResidual->Write();
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
