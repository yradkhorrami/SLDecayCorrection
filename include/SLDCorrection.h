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
	TLorentzVector getNeutrinoFourMomentum( TVector3 flightDirection , TLorentzVector fourMomentumLepton , TLorentzVector visibleFourMomentumCharged , TLorentzVector visibleFourMomentumNeutral , double parentHadronMass , int solutionSign );
	TLorentzVector getTrueNeutrinoFourMomentum( MCParticle *SLDLepton );
	virtual void plotHistograms( TLorentzVector trueFourMomentumNeutrino , TLorentzVector FourMomentumNuClose , std::vector<float> NeutrinoCovMat );
	virtual void InitializeHistogram( TH1F *histogram , int scale , int color , int lineWidth , int markerSize , int markerStyle );
	virtual void doProperGaussianFit( TH1F *histogram , float fitMin , float fitMax , float fitRange );

	virtual void check( EVENT::LCEvent *pLCEvent );
	virtual void end();
	dd4hep::Detector& _theDetector = dd4hep::Detector::getInstance();

private:

	typedef std::vector<int>		IntVector;
	typedef std::vector<double>		DoubleVector;
	typedef std::vector<float>		FloatVector;

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
	IntVector				m_foundRecoLepton{};
	IntVector				m_foundBuildUpVertex{};
	IntVector				m_foundRecoLeptonInBuildUpVertex{};
	IntVector				m_foundRecoLeptonInPrimaryVertex{};
	DoubleVector				m_lostChargedMCP_CosTheta{};
	DoubleVector				m_lostChargedMCP_Energy{};
	DoubleVector				m_lostChargedMCP_Pt{};
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
	TTree					*m_pTTree{};


};
#endif
