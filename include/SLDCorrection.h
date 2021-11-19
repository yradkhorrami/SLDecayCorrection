#ifndef SLDCorrection_h
#define SLDCorrection_h 1
#include <marlin/Processor.h>
#include <marlin/Global.h>
#include "EVENT/LCStrVec.h"
#include <EVENT/MCParticle.h>
#include <EVENT/Track.h>
#include "IMPL/LCCollectionVec.h"
#include "UTIL/LCRelationNavigator.h"
#include <IMPL/ReconstructedParticleImpl.h>
#include "lcio.h"
#include <string>
#include <vector>
#include <math.h>
#include <set>
#include <vector>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMatrixD.h"
#include "DD4hep/Detector.h"
class TFile;
class TDirectory;
class TH1F;
class TH1I;
class TH2I;
class TH2F;
class TTree;
class TF1;

using namespace lcio ;
using namespace marlin ;

class SLDCorrection : public Processor
{
public:
	virtual Processor *newProcessor()
	{
		return new SLDCorrection;
	}
	SLDCorrection();
	virtual ~SLDCorrection() = default;
	SLDCorrection( const SLDCorrection& ) = delete;
	SLDCorrection &operator = ( const SLDCorrection& ) = delete;
	virtual void init();
	virtual void Clear();
	virtual void processRunHeader();
	virtual void processEvent( EVENT::LCEvent *pLCEvent );

	bool hasPrimarySLDecay( MCParticle *parentHadron );
	bool hasDownStreamSLDecay( MCParticle *parentHadron );
	bool hasUpStreamSLDecay( MCParticle *parentHadron );
	bool checkBHadronSLDecay( MCParticle *SLDLepton );
	bool checkCHadronSLDecay( MCParticle *SLDLepton );
	bool checkTauLeptonSLDecay( MCParticle *SLDLepton );
	virtual void doSLDCorrection( EVENT::LCEvent *pLCEvent , MCParticle *SLDLepton );
	void showTrueParameters( MCParticle *SLDLepton );
	TLorentzVector getNeutrinoFourMomentum( TVector3 flightDirection , TLorentzVector visibleFourMomentum , double parentHadronMass , float solutionSign );
	MCParticle* getTrueNeutrino( MCParticle *SLDLepton , TLorentzVector& InisibleFourMomentum );
	void fillTrueRecoFourMomentum( TLorentzVector trueNeutralFourMomentum , TLorentzVector trueChargedFourMomentum , TLorentzVector trueLeptonFourMomentum , TLorentzVector trueVisibleFourMomentum , TLorentzVector trueNeutrinoFourMomentum , TLorentzVector trueHadronFourMomentum ,  TLorentzVector recoNeutralFourMomentum , TLorentzVector recoChargedFourMomentum , TLorentzVector recoLeptonFourMomentum , TLorentzVector recoVisibleFourMomentum , TLorentzVector recoNeutrinoFourMomentum , TLorentzVector recoHadronFourMomentum );
	virtual void plotHistograms( TLorentzVector trueFourMomentumNeutrino , TLorentzVector FourMomentumNuClose , std::vector<float> NeutrinoCovMat );
	virtual void InitializeHistogram( TH1F *histogram , int scale , int color , int lineWidth , int markerSize , int markerStyle );
	virtual void doProperGaussianFit( TH1F *histogram , float fitMin , float fitMax , float fitRange );
	int getVertexInJetsDistribution( Vertex* testVertex , std::vector<EVENT::ReconstructedParticle*> jetVector );

	virtual void check( EVENT::LCEvent *pLCEvent );
	virtual void end();
	dd4hep::Detector& _theDetector = dd4hep::Detector::getInstance();

private:

	typedef std::vector<int>				IntVector;
	typedef std::vector<double>				DoubleVector;
	typedef std::vector<float>				FloatVector;
	typedef std::vector<EVENT::MCParticle*>			mcpVector;
	typedef std::vector<EVENT::ReconstructedParticle*>	pfoVector;
	typedef std::vector<EVENT::Vertex*>			vtxVector;


	std::string				m_mcParticleCollection{};
	std::string				m_inputPfoCollection{};
	std::string				m_inputJetCollection{};
	std::string				m_jetFlavour{};
	std::string				m_inputPrimaryVertex{};
	std::string				m_inputBuildUpVertex{};
	std::string				m_RecoMCTruthLinkCollection{};
	std::string				m_MCTruthRecoLinkCollection{};
	std::string				m_TrackMCTruthLinkCollection{};
	std::string				m_MCTruthTrackLinkCollection{};
	std::string				m_ClusterMCTruthLinkCollection{};
	std::string				m_MCTruthClusterLinkCollection{};
	std::string				m_SLDNuCollection{};
	std::string				m_rootFile{};

	bool					m_includbJets = true;
	bool					m_includcJets = true;
	bool					m_includgJets = true;
	bool					m_includOthers = true;
	bool					m_includeBSLD = true;
	bool					m_includeCSLD = true;
	bool					m_includeTSLD = true;
	bool					m_cheatSLDLeptons = true;
	bool					m_cheatFlightDirection = true;
	int					m_vertexingScenario = 1;
	bool					m_useJetAxisAsFlightDirection = true;
	bool					m_considerParentCharge = true;
	bool					m_cheatVertices = true;
	bool					m_cheatLepton4momentum = true;
	bool					m_cheatCharged4momentum = true;
	bool					m_cheatNeutral4momentum = true;
	bool					m_cheatPIDcharged = true;
	bool					m_cheatPIDneutral = true;
	int					m_nIterFlightDirCorrection = 0;
	int					m_recoFourMomentumOfVisibles = 0;
	bool					m_displayEvent = true;
	bool					m_fillRootTree = true;

	int					m_nRun;
	int					m_nEvt;
	int					m_nRunSum;
	int					m_nEvtSum;
	double					m_Bfield;
	double					c;
	double					mm2m;
	double					eV2GeV;
	double					eB;
	bool					foundFlightDirection;
	int					m_nTauSLDecay;
	int					m_nTauNeutrino;
	int					m_nNeutrino;
	int					m_nChargedPFOwoTrack;
	IntVector				m_jetFlavourPDG{};
	IntVector				m_nSLD_chargedMCPwoTrack{};
	IntVector				m_GenStatParentHadron{};
	IntVector				m_ChargeParentHadron{};
	IntVector				m_nTrueNeutralDecayProducts{};
	IntVector				m_nTrueAloneChargedDecayProducts{};
	IntVector				m_nTrueVertices{};
	IntVector				m_nJetsVerticesDistributedIn{};
	IntVector				m_nRecoVerticesInJet{};
	IntVector				m_nAloneChargedPFOs{};
	IntVector				m_nAloneChargedPFOsFromSLD{};
	DoubleVector				m_distLeptonAlonePFOsNotFromSLD{};
	DoubleVector				m_distLeptonAlonePFOsFromSLD{};
	IntVector				m_foundRecoLepton{};
	IntVector				m_foundBuildUpVertex{};
	IntVector				m_foundRecoLeptonInBuildUpVertex{};
	IntVector				m_foundRecoLeptonInPrimaryVertex{};
	DoubleVector				m_lostChargedMCP_CosTheta{};
	DoubleVector				m_lostChargedMCP_Energy{};
	DoubleVector				m_lostChargedMCP_Pt{};
	int					n_SLDStatus;
	int					n_NuPxResidual;
	int					n_NuPyResidual;
	int					n_NuPzResidual;
	int					n_NuEResidual;
	int					n_NuPxNormalizedResidual;
	int					n_NuPyNormalizedResidual;
	int					n_NuPzNormalizedResidual;
	int					n_NuENormalizedResidual;
	int					n_secondaryVertex;
	DoubleVector				m_SLDecayXi{};
	DoubleVector				m_SLDecayYi{};
	DoubleVector				m_SLDecayZi{};
	DoubleVector				m_SLDecayRi{};
	DoubleVector				m_SLDecayXf{};
	DoubleVector				m_SLDecayYf{};
	DoubleVector				m_SLDecayZf{};
	DoubleVector				m_SLDecayRf{};
	DoubleVector				m_trueNuPx{};
	DoubleVector				m_trueNuPy{};
	DoubleVector				m_trueNuPz{};
	DoubleVector				m_trueNuE{};
	DoubleVector				m_recoNuCloseInitialPx{};
	DoubleVector				m_recoNuCloseInitialPy{};
	DoubleVector				m_recoNuCloseInitialPz{};
	DoubleVector				m_recoNuCloseInitialE{};
	DoubleVector				m_recoNuClosePx{};
	DoubleVector				m_recoNuClosePy{};
	DoubleVector				m_recoNuClosePz{};
	DoubleVector				m_recoNuCloseE{};
	DoubleVector				m_recoNuPosPx{};
	DoubleVector				m_recoNuPosPy{};
	DoubleVector				m_recoNuPosPz{};
	DoubleVector				m_recoNuPosE{};
	DoubleVector				m_recoNuNegPx{};
	DoubleVector				m_recoNuNegPy{};
	DoubleVector				m_recoNuNegPz{};
	DoubleVector				m_recoNuNegE{};
	DoubleVector				m_NuPxResidual{};
	DoubleVector				m_NuPyResidual{};
	DoubleVector				m_NuPzResidual{};
	DoubleVector				m_NuEResidual{};
	DoubleVector				m_NuPxNormalizedResidual{};
	DoubleVector				m_NuPyNormalizedResidual{};
	DoubleVector				m_NuPzNormalizedResidual{};
	DoubleVector				m_NuENormalizedResidual{};
	IntVector				m_solutionSign{};
	DoubleVector				m_E_vis{};
	DoubleVector				m_E_vis_prime{};
	DoubleVector				m_P_vis_par{};
	DoubleVector				m_P_vis_par_prime{};
	DoubleVector				m_P_vis_nor{};
	DoubleVector				m_P_vis_nor_prime{};
	IntVector				m_flightDirectionStatus{};
	DoubleVector				m_FlightDirectionErrorSinAlpha{};
	DoubleVector				m_FlightDirectionErrorCosAlpha{};
	DoubleVector				m_FlightDirectionErrorAlpha{};
	DoubleVector				m_distRecoLeptonToDownStreamVertex{};
	DoubleVector				m_dsVertexResidualX{};
	DoubleVector				m_dsVertexResidualY{};
	DoubleVector				m_dsVertexResidualZ{};
	DoubleVector				m_SecVertexResidualX{};
	DoubleVector				m_SecVertexResidualY{};
	DoubleVector				m_SecVertexResidualZ{};
	DoubleVector				m_parentHadronMass{};
	DoubleVector				m_parentHadronFlightDistance{};
	DoubleVector				m_daughterHadronMass{};
	DoubleVector				m_widestConeAlphaNeutrals{};
	DoubleVector				m_widestConeAlphaVertices{};
	DoubleVector				m_widestConeAlphaCharged{};
	DoubleVector				m_widestConeCosAlphaNeutrals{};
	DoubleVector				m_widestConeCosAlphaVertices{};
	DoubleVector				m_widestConeCosAlphaCharged{};
	DoubleVector				m_widestConeAlphaAlonePFOsFromSLDwrtLepton{};
	DoubleVector				m_widestConeCosAlphaAlonePFOsFromSLDwrtLepton{};
	DoubleVector				m_widestConeAlphaAlonePFOsFromSLDwrtFD{};
	DoubleVector				m_widestConeCosAlphaAlonePFOsFromSLDwrtFD{};
	DoubleVector				m_widestConeAlphaAlonePFOsFromSLDwrtJet{};
	DoubleVector				m_widestConeCosAlphaAlonePFOsFromSLDwrtJet{};
	DoubleVector				m_widestConeAlphaAlonePFOsNotFromSLDwrtLepton{};
	DoubleVector				m_widestConeCosAlphaAlonePFOsNotFromSLDwrtLepton{};
	DoubleVector				m_widestConeAlphaAlonePFOsNotFromSLDwrtFD{};
	DoubleVector				m_widestConeCosAlphaAlonePFOsNotFromSLDwrtFD{};
	DoubleVector				m_widestConeAlphaAlonePFOsNotFromSLDwrtJet{};
	DoubleVector				m_widestConeCosAlphaAlonePFOsNotFromSLDwrtJet{};

	DoubleVector				m_trueNeutralPx{};
	DoubleVector				m_trueNeutralPy{};
	DoubleVector				m_trueNeutralPz{};
	DoubleVector				m_trueNeutralE{};
	DoubleVector				m_trueChargedPx{};
	DoubleVector				m_trueChargedPy{};
	DoubleVector				m_trueChargedPz{};
	DoubleVector				m_trueChargedE{};
	DoubleVector				m_trueLeptonPx{};
	DoubleVector				m_trueLeptonPy{};
	DoubleVector				m_trueLeptonPz{};
	DoubleVector				m_trueLeptonE{};
	DoubleVector				m_trueVisiblePx{};
	DoubleVector				m_trueVisiblePy{};
	DoubleVector				m_trueVisiblePz{};
	DoubleVector				m_trueVisibleE{};
	DoubleVector				m_trueNeutrinoPx{};
	DoubleVector				m_trueNeutrinoPy{};
	DoubleVector				m_trueNeutrinoPz{};
	DoubleVector				m_trueNeutrinoE{};
	DoubleVector				m_trueHadronPx{};
	DoubleVector				m_trueHadronPy{};
	DoubleVector				m_trueHadronPz{};
	DoubleVector				m_trueHadronE{};
	DoubleVector				m_recoNeutralPx{};
	DoubleVector				m_recoNeutralPy{};
	DoubleVector				m_recoNeutralPz{};
	DoubleVector				m_recoNeutralE{};
	DoubleVector				m_recoChargedPx{};
	DoubleVector				m_recoChargedPy{};
	DoubleVector				m_recoChargedPz{};
	DoubleVector				m_recoChargedE{};
	DoubleVector				m_recoLeptonPx{};
	DoubleVector				m_recoLeptonPy{};
	DoubleVector				m_recoLeptonPz{};
	DoubleVector				m_recoLeptonE{};
	DoubleVector				m_recoVisiblePx{};
	DoubleVector				m_recoVisiblePy{};
	DoubleVector				m_recoVisiblePz{};
	DoubleVector				m_recoVisibleE{};
	DoubleVector				m_recoNeutrinoPx{};
	DoubleVector				m_recoNeutrinoPy{};
	DoubleVector				m_recoNeutrinoPz{};
	DoubleVector				m_recoNeutrinoE{};
	DoubleVector				m_recoHadronPx{};
	DoubleVector				m_recoHadronPy{};
	DoubleVector				m_recoHadronPz{};
	DoubleVector				m_recoHadronE{};

	IntVector				m_SLDStatus{};
	FloatVector				m_weightPFOtoMCP_Lepton{};
	FloatVector				m_weightMCPtoPFO_Lepton{};
	FloatVector				m_weightPFOtoMCP_Neutral{};
	FloatVector				m_weightMCPtoPFO_Neutral{};
	FloatVector				m_weightPFOtoMCP_Charged{};
	FloatVector				m_weightMCPtoPFO_Charged{};
	TH1F					*h_SLDStatus{};
	TH1F					*h_NuPxResidual{};
	TH1F					*h_NuPyResidual{};
	TH1F					*h_NuPzResidual{};
	TH1F					*h_NuEResidual{};
	TH1F					*h_NuPxNormalizedResidual{};
	TH1F					*h_NuPyNormalizedResidual{};
	TH1F					*h_NuPzNormalizedResidual{};
	TH1F					*h_NuENormalizedResidual{};
	TH2F					*h_recoNuPx_mcNuPx{};
	TH2F					*h_recoNuPy_mcNuPy{};
	TH2F					*h_recoNuPz_mcNuPz{};
	TH2F					*h_recoNuE_mcNuE{};
	TH2F					*h_parentPx_daughtersPx{};
	TH2F					*h_parentPy_daughtersPy{};
	TH2F					*h_parentPz_daughtersPz{};
	TH2F					*h_parentE_daughtersE{};
	TH1I					*h_recoPFOLinkedToElectron_Type{};
	TH1I					*h_recoPFOLinkedToMuon_Type{};
	TH1I					*h_SLDecayOrder{};
	TH2I					*h_foundVertex{};
	TH1F					*h_secondaryVertex{};
	TH1I					*h_parentHadronCharge{};
	TH1I					*h_MCPTracks{};
	TH1I					*h_MCPTracks_Eweighted{};
	TH1I					*h_MCPTracks_Ptweighted{};
	TH1F					*h_FlightDirectionError{};
	TH1F					*h_distRecoLeptonToDownStreamVertex{};
	TFile					*m_pTFile{};
	TTree					*m_pTTree1{};
	TTree					*m_pTTree2{};


};
#endif
