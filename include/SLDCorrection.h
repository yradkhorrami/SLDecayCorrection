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
#include "IMPL/ParticleIDImpl.h"
#include <IMPL/VertexImpl.h>
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

	typedef std::vector<int>				IntVector;
	typedef std::vector<double>				DoubleVector;
	typedef std::vector<float>				FloatVector;
	typedef std::vector<EVENT::MCParticle*>			mcpVector;
	typedef std::vector<EVENT::ReconstructedParticle*>	pfoVector;
	typedef std::vector<std::vector<EVENT::ReconstructedParticle*>>	pfoVectorVector;
	typedef std::vector<EVENT::Vertex*>			vtxVector;
	bool hasPrimarySLDecay( MCParticle *parentHadron );
	bool hasDownStreamSLDecay( MCParticle *parentHadron );
	bool hasUpStreamSLDecay( MCParticle *parentHadron );
	bool checkBHadronSLDecay( MCParticle *SLDLepton );
	bool checkCHadronSLDecay( MCParticle *SLDLepton );
	bool checkTauLeptonSLDecay( MCParticle *SLDLepton );
	virtual void doSLDCorrection( EVENT::LCEvent *pLCEvent , MCParticle *SLDLepton , vtxVector& semiLeptonicVertices , pfoVector& semiLeptonicVertexRecoParticles , pfoVector& jetsOfSemiLeptonicDecays , pfoVectorVector& neutrinos , IntVector &sldStatus , IntVector &pvaStatus , IntVector &solutionSigns , mcpVector &trueNeutrinos );
	void showTrueParameters( MCParticle *SLDLepton );
	TLorentzVector getNeutrinoFourMomentum( TVector3 flightDirection , TLorentzVector visibleFourMomentum , double parentHadronMass , float solutionSign );
	TLorentzVector getNeutrinoFourMomentumModified( TVector3 &flightDirection , TLorentzVector visibleFourMomentum , double parentHadronMass , float solutionSign );
	TLorentzVector getNeutrinoFourMomentumStandardMethod( TVector3 flightDirection , TLorentzVector visibleFourMomentum , double parentHadronMass , float solutionSign );
	MCParticle* getTrueNeutrino( MCParticle *SLDLepton );
	void fillTrueRecoFourMomentum(	TLorentzVector trueVisibleFourMomentumAtSLDVertex , TLorentzVector truePVATrueFourMomentum , TLorentzVector truePVARecoFourMomentum , TLorentzVector recoPVARecoFourMomentum , TLorentzVector visibleFourMomentum ,
								TLorentzVector trueLeptonFourMomentum , TLorentzVector recoLeptonFourMomentum , TLorentzVector leptonFourMomentum ,
								TLorentzVector truePVATrueChargedFourMomentum , TLorentzVector truePVARecoChargedFourMomentum , TLorentzVector recoSLDVertexChargedFourMomentum , TLorentzVector recoPVARecoChargedFourMomentum , TLorentzVector chargedFourMomentum ,
								TLorentzVector truePVATrueNeutralFourMomentum , TLorentzVector truePVARecoNeutralFourMomentum , TLorentzVector recoPVARecoNeutralFourMomentum , TLorentzVector neutralFourMomentum ,
								TLorentzVector trueNeutrinoFourMomentum , TLorentzVector recoNeutrinoFourMomentumClose , TLorentzVector trueHadronFourMomentum , TLorentzVector recoHadronFourMomentum ,
								TVector3 trueFlightDirection , TVector3 recoFlightDirection , TVector3 flightDirection );

//	virtual void getCovMatPVA( std::vector<EVENT::ReconstructedParticle*> decayProducts , std::vector<EVENT::ReconstructedParticle*> associatedParticles , int SLDStatus , std::vector< float > &CovMatrixPVA );
//	virtual void getCovMatDetFlightDirection( std::vector<EVENT::ReconstructedParticle*> decayProducts , TVector3 flightDirection , std::vector< float > CovMatrixFlightDirection , std::vector< float > &CovMatrixDetector , EVENT::ReconstructedParticle* linkedRecoLepton , std::vector< float > &CovMatrixDetPar , std::vector< float > &CovMatrixDetNor );
	virtual void getCovMatFlightDirection( TVector3 flightDirection , float sigmaAlpha , std::vector< float > &CovMatrixFlightDirection );
//	virtual void getCovMatrixDetPar( TVector3 flightDirection , TLorentzVector visibleFourMomentum , std::vector< float > CovMatrixFlightDirection , std::vector< float > initialCovMatrixDetector , std::vector< float > &CovMatrixDetPar );
//	virtual void getCovMatrixDetNor( TVector3 flightDirection , TLorentzVector visibleFourMomentum , std::vector< float > CovMatrixFlightDirection , std::vector< float > initialCovMatrixDetector , std::vector< float > &CovMatrixDetNor );
//	virtual void getNeutrinoCovMat( TLorentzVector recoNeutrinoFourMomentum , TLorentzVector visibleFourMomentum , TVector3 flightDirection , double parentHadronMass , std::vector< float > CovMatrixPVA , std::vector< float > CovMatrixDetector , std::vector< float > CovMatrixDetPar , std::vector< float > CovMatrixDetNor , std::vector< float > &NeutrinoCovMatrix );

//	virtual void getNeutrinoCovMatET( std::vector<EVENT::ReconstructedParticle*> decayProducts , std::vector<EVENT::ReconstructedParticle*> associatedParticles , EVENT::ReconstructedParticle* linkedRecoLepton , TVector3 flightDirection , TLorentzVector recoNeutrinoFourMomentum , std::vector< float > &NeutrinoCovMatrix , int SLDStatus );
//	virtual void getNeutrinoCovarianceMatrix( std::vector< float > &NeutrinoCovMatrix , int SLDStatus , std::vector<EVENT::ReconstructedParticle*> decayProducts , std::vector<EVENT::ReconstructedParticle*> associatedParticles , EVENT::ReconstructedParticle* linkedRecoLepton , TVector3 flightDirection , TLorentzVector recoNeutrinoFourMomentum , double parentHadronMass , float solutionSign );
	virtual void getNeutrinoCovarianceMatrix( std::vector< float > &NeutrinoCovMatrixPos , std::vector< float > &NeutrinoCovMatrixNeg , int SLDStatus , EVENT::MCParticle *SLDLepton , EVENT::ReconstructedParticle *linkedRecoLepton , mcpVector trueChargedDecayProducts , mcpVector trueNeutralDecayProducts , std::vector<EVENT::ReconstructedParticle*> truePVAChargedDecayProducts , std::vector<EVENT::ReconstructedParticle*> truePVANeutralDecayProducts , std::vector<EVENT::ReconstructedParticle*> recoSLDVertexDecayProducts , std::vector<EVENT::ReconstructedParticle*> recoPVAVertexDecayProducts , std::vector<EVENT::ReconstructedParticle*> recoPVAChargedDecayProducts , std::vector<EVENT::ReconstructedParticle*> recoPVANeutralDecayProducts , TVector3 flightDirection , TLorentzVector recoNeutrinoFourMomentumPos , TLorentzVector recoNeutrinoFourMomentumNeg , double parentHadronMass , std::vector< float > &CovMatDetector , std::vector< float > &CovMatFlightDirection );
	virtual void addNeutrinoCovarianceMatrix( TLorentzVector neutrinoFourMomentum , std::vector< float > &NuCovMat );
	virtual void evaluateInputCovMat( TLorentzVector trueVisibleFourMomentum , TVector3 trueFlightDirection , TLorentzVector trueNeutrinoFourMomentum , TLorentzVector visibleFourMomentum , TVector3 flightDirection , TLorentzVector recoNeutrinoFourMomentum , std::vector< float > CovMatDetector , std::vector< float > CovMatFlightDirection , std::vector< float > CovMatNeutrino );
	virtual void getNuCovMatP4( TLorentzVector visibleFourMomentum , TVector3 flightDirection , double parentHadronMass , std::vector< float > CovMatrixDetectorPVA , std::vector< float > CovMatrixFlightDirection , TLorentzVector recoNeutrinoFourMomentum , std::vector< float > &NeutrinoCovMatrix , float solutionSign );
	virtual void getNuCovMatP4( TLorentzVector visibleFourMomentum , TVector3 flightDirection , double parentHadronMass , std::vector< float > CovMatrixDetector , TLorentzVector recoNeutrinoFourMomentum , std::vector< float > &NeutrinoCovMatrix , float solutionSign );
//	virtual void getNuCovMatP4( std::vector< float > &NeutrinoCovMatrix , float NeutrinoEnergyResolution , TLorentzVector recoNeutrinoFourMomentum , TLorentzVector visibleFourMomentum , TVector3 flightDirection , std::vector< float > CovMatrixDetectorPVA , std::vector< float > CovMatrixFlightDirection );
	virtual void getCovMatPVA( TLorentzVector recoNeutrinoFourMomentum , int PVAStatus , std::vector< float > &CovMatrixPVA );
	virtual void getPVACovMat( std::vector<EVENT::ReconstructedParticle*> recoSLDVertexDecayProducts , std::vector<EVENT::ReconstructedParticle*> recoPVAVertexDecayProducts , std::vector<EVENT::ReconstructedParticle*> recoPVAChargedDecayProducts , std::vector<EVENT::ReconstructedParticle*> recoPVANeutralDecayProducts , int SLDStatus , std::vector< float > &CovMatrixChargedPVA , std::vector< float > &CovMatrixNeutralPVA );
	virtual void getCovarianceMatrixPVA( TLorentzVector PVAFourMomentum , float sigmaAlpha , float sigmaEnergy , std::vector< float > &CovMatrixPVA );

	virtual void plotHistograms( TLorentzVector trueFourMomentumNeutrino , TLorentzVector FourMomentumNuClose , std::vector<float> NeutrinoCovMat );
	virtual void InitializeHistogram( TH1F *histogram , int scale , int color , int lineWidth , int markerSize , int markerStyle );
	virtual void doProperGaussianFit( TH1F *histogram , float fitMin , float fitMax , float fitRange );
	int getVertexInJetsDistribution( Vertex* testVertex , std::vector<EVENT::ReconstructedParticle*> jetVector );
	void investigateJetEnergyContent( EVENT::ReconstructedParticle *assignedJet );
	void checkSLDInput( MCParticle *SLDHadron );

	virtual void check( EVENT::LCEvent *pLCEvent );
	virtual void end();
	dd4hep::Detector& _theDetector = dd4hep::Detector::getInstance();

private:



	std::string				m_mcParticleCollection{};
	std::string				m_inputPfoCollection{};
	std::string				m_inputJetCollection{};
	std::string				m_inputPrimaryVertex{};
	std::string				m_inputBuildUpVertex{};
	std::string				m_RecoMCTruthLinkCollection{};
	std::string				m_MCTruthRecoLinkCollection{};
	std::string				m_ClusterMCTruthLinkCollection{};
	std::string				m_MCTruthClusterLinkCollection{};
	std::string				m_SLDVertex{};
	std::string				m_SLDVertexRP{};
	std::string				m_reconstructedNeutrino{};
	std::string				m_JetSLDLinkName{};
	std::string				m_SLDJetLinkName{};
	std::string				m_mcNurecoNuLinkName{};
	std::string				m_recoNumcNuLinkName{};
	std::string				m_SLDNuLinkName{};
	std::string				m_NuSLDLinkName{};
	std::string				m_rootFile{};

	bool					m_includeBSLD = true;
	bool					m_includeCSLD = true;
	bool					m_includeTSLD = true;
	bool					m_cheatSLDLeptons = true;
	bool					m_cheatFlightDirection = true;
	int						m_vertexingScenario = 1;
	bool					m_cheatLepton4momentum = true;
	bool					m_cheatCharged4momentum = true;
	bool					m_cheatNeutral4momentum = true;
	bool					m_cheatPVAcharged = true;
	float					m_chargedCosAcceptanceAngleSLD4 = 0.0;
	float					m_chargedCosAcceptanceAngleSLD5 = 0.0;
	bool					m_cheatPVAneutral = true;
	bool					m_cheatSolutionSign = false;
	float					m_neutralCosAcceptanceAngle = 0.0;
	float					m_BSLDChargedSLD4InvMassCut = 0.0;
	float					m_BSLDChargedSLD5InvMassCut = 0.0;
	float					m_BSLDNeutralSLD4InvMassCut = 0.0;
	float					m_BSLDNeutralSLD5InvMassCut = 0.0;
	float					m_CSLDChargedSLD4InvMassCut = 0.0;
	float					m_CSLDChargedSLD5InvMassCut = 0.0;
	float					m_CSLDNeutralSLD4InvMassCut = 0.0;
	float					m_CSLDNeutralSLD5InvMassCut = 0.0;
	int						m_nIterFlightDirCorrection = 0;
	float					m_BSLD4SigmaAlpha = 0.004;
	float					m_BSLD5SigmaAlpha = 0.010;
	float					m_BSLD4SigmaECPVA = 0;
	float					m_BSLD4SigmaENPVA = 0;
	float					m_BSLD4SigmaAlphaCPVA = 0.0;
	float					m_BSLD4SigmaAlphaNPVA = 0.0;
	float					m_BSLD5SigmaECPVA = 0;
	float					m_BSLD5SigmaENPVA = 0;
	float					m_BSLD5SigmaAlphaCPVA = 0.0;
	float					m_BSLD5SigmaAlphaNPVA = 0.0;
	float					m_sigmaAlphaNu = 0.1;
	float					m_sigmaENu = 4.0;
	bool					m_displayEvent = true;
	bool					m_fillRootTree = true;
	bool					m_traceEvent = false;
	int						m_BSLDMode = 0;
	int						m_CSLDMode = 0;
	int						m_TSLDMode = 0;
	int						m_SLDMode = 0;

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
	IntVector				m_SLDFlavour{}; //4: SLDecayOfCHadron, 5: SLDecayOfBHadron, 15: SLDecayOfTauLepton
	IntVector				m_SLDType{}; //0: SLDecay with DownStream/UpStream semi-leptonic decay(s), 1: SLDecay without DownStream/UpStream semi-leptonic decay(s)
	IntVector				m_SLDLeptonID{}; //+11/-11: SLDecay to electron/positron, +13/-13: SLDecay to muon/anti-muon, +15/-15: SLDecay to tau/anti-tau
	FloatVector				m_leptonE_to_parentE{};
	FloatVector				m_otherChargedE_to_parentE{};
	FloatVector				m_allChargedE_to_parentE{};
	FloatVector				m_neutralE_to_parentE{};
	FloatVector				m_neutrino_to_parentE{};
	unsigned int			m_nSLDecayOfBHadron;
	unsigned int			m_nSLDecayOfCHadron;
	unsigned int			m_nSLDecayOfTauLepton;
	unsigned int			m_nSLDecayTotal;
	int						m_nSLDecayToElectron;
	int						m_nSLDecayToMuon;
	int						m_nSLDecayToTau;
	int						m_nTauNeutrino;
	int						m_nNeutrino;
	int						m_nChargedPFOwoTrack;
	IntVector				BHadPDGs{};
	IntVector				CHadPDGs{};
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
	DoubleVector				m_jetEnergy{};
	DoubleVector				m_jetEnergyFractionCharged{};
	DoubleVector				m_jetEnergyFractionNeutralHadron{};
	DoubleVector				m_jetEnergyFractionPhoton{};
	DoubleVector				m_jetEnergyFractionNeutrals{};
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
	DoubleVector				m_true_E_vis{};
	DoubleVector				m_true_E_vis_prime{};
	DoubleVector				m_true_P_vis_par{};
	DoubleVector				m_true_P_vis_par_prime{};
	DoubleVector				m_true_P_vis_nor{};
	DoubleVector				m_E_vis{};
	DoubleVector				m_E_vis_prime{};
	DoubleVector				m_P_vis_par{};
	DoubleVector				m_P_vis_par_prime{};
	DoubleVector				m_P_vis_par_prime_squared{};
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
	IntVector					m_parentHadronPDG{};
	DoubleVector				m_trueParentHadronFlightDistance{};
	DoubleVector				m_recoParentHadronFlightDistance{};
	DoubleVector				m_daughterHadronMass{};
	IntVector					m_daughterHadronPDG{};
	DoubleVector				m_daughterHadronFlightDistance{};
	FloatVector					m_visibleChargedInvMassCut{};
	FloatVector					m_visibleNeutralInvMassCut{};

	DoubleVector				m_trueVisibleFourMomentumAtSLDVertex_Px{};
	DoubleVector				m_trueVisibleFourMomentumAtSLDVertex_Py{};
	DoubleVector				m_trueVisibleFourMomentumAtSLDVertex_Pz{};
	DoubleVector				m_trueVisibleFourMomentumAtSLDVertex_E{};
	DoubleVector				m_trueVisibleFourMomentumAtSLDVertex_M{};
	DoubleVector				m_truePVATrueFourMomentum_Px{};
	DoubleVector				m_truePVATrueFourMomentum_Py{};
	DoubleVector				m_truePVATrueFourMomentum_Pz{};
	DoubleVector				m_truePVATrueFourMomentum_E{};
	DoubleVector				m_truePVATrueFourMomentum_M{};
	DoubleVector				m_truePVARecoFourMomentum_Px{};
	DoubleVector				m_truePVARecoFourMomentum_Py{};
	DoubleVector				m_truePVARecoFourMomentum_Pz{};
	DoubleVector				m_truePVARecoFourMomentum_E{};
	DoubleVector				m_truePVARecoFourMomentum_M{};
	DoubleVector				m_recoPVARecoFourMomentum_Px{};
	DoubleVector				m_recoPVARecoFourMomentum_Py{};
	DoubleVector				m_recoPVARecoFourMomentum_Pz{};
	DoubleVector				m_recoPVARecoFourMomentum_E{};
	DoubleVector				m_recoPVARecoFourMomentum_M{};
	DoubleVector				m_usedVisibleFourMomentum_Px{};
	DoubleVector				m_usedVisibleFourMomentum_Py{};
	DoubleVector				m_usedVisibleFourMomentum_Pz{};
	DoubleVector				m_usedVisibleFourMomentum_E{};
	DoubleVector				m_usedVisibleFourMomentum_M{};
	DoubleVector				m_trueLeptonFourMomentum_Px{};
	DoubleVector				m_trueLeptonFourMomentum_Py{};
	DoubleVector				m_trueLeptonFourMomentum_Pz{};
	DoubleVector				m_trueLeptonFourMomentum_E{};
	DoubleVector				m_trueLeptonFourMomentum_M{};
	DoubleVector				m_recoLeptonFourMomentum_Px{};
	DoubleVector				m_recoLeptonFourMomentum_Py{};
	DoubleVector				m_recoLeptonFourMomentum_Pz{};
	DoubleVector				m_recoLeptonFourMomentum_E{};
	DoubleVector				m_recoLeptonFourMomentum_M{};
	DoubleVector				m_usedLeptonFourMomentum_Px{};
	DoubleVector				m_usedLeptonFourMomentum_Py{};
	DoubleVector				m_usedLeptonFourMomentum_Pz{};
	DoubleVector				m_usedLeptonFourMomentum_E{};
	DoubleVector				m_usedLeptonFourMomentum_M{};
	DoubleVector				m_truePVATrueChargedFourMomentum_Px{};
	DoubleVector				m_truePVATrueChargedFourMomentum_Py{};
	DoubleVector				m_truePVATrueChargedFourMomentum_Pz{};
	DoubleVector				m_truePVATrueChargedFourMomentum_E{};
	DoubleVector				m_truePVATrueChargedFourMomentum_M{};
	DoubleVector				m_truePVARecoChargedFourMomentum_Px{};
	DoubleVector				m_truePVARecoChargedFourMomentum_Py{};
	DoubleVector				m_truePVARecoChargedFourMomentum_Pz{};
	DoubleVector				m_truePVARecoChargedFourMomentum_E{};
	DoubleVector				m_truePVARecoChargedFourMomentum_M{};
	DoubleVector				m_recoSLDVertexChargedFourMomentum_Px{};
	DoubleVector				m_recoSLDVertexChargedFourMomentum_Py{};
	DoubleVector				m_recoSLDVertexChargedFourMomentum_Pz{};
	DoubleVector				m_recoSLDVertexChargedFourMomentum_E{};
	DoubleVector				m_recoSLDVertexChargedFourMomentum_M{};
	DoubleVector				m_recoPVARecoChargedFourMomentum_Px{};
	DoubleVector				m_recoPVARecoChargedFourMomentum_Py{};
	DoubleVector				m_recoPVARecoChargedFourMomentum_Pz{};
	DoubleVector				m_recoPVARecoChargedFourMomentum_E{};
	DoubleVector				m_recoPVARecoChargedFourMomentum_M{};
	DoubleVector				m_usedChargedFourMomentum_Px{};
	DoubleVector				m_usedChargedFourMomentum_Py{};
	DoubleVector				m_usedChargedFourMomentum_Pz{};
	DoubleVector				m_usedChargedFourMomentum_E{};
	DoubleVector				m_usedChargedFourMomentum_M{};
	DoubleVector				m_truePVATrueNeutralFourMomentum_Px{};
	DoubleVector				m_truePVATrueNeutralFourMomentum_Py{};
	DoubleVector				m_truePVATrueNeutralFourMomentum_Pz{};
	DoubleVector				m_truePVATrueNeutralFourMomentum_E{};
	DoubleVector				m_truePVATrueNeutralFourMomentum_M{};
	DoubleVector				m_truePVARecoNeutralFourMomentum_Px{};
	DoubleVector				m_truePVARecoNeutralFourMomentum_Py{};
	DoubleVector				m_truePVARecoNeutralFourMomentum_Pz{};
	DoubleVector				m_truePVARecoNeutralFourMomentum_E{};
	DoubleVector				m_truePVARecoNeutralFourMomentum_M{};
	DoubleVector				m_recoPVARecoNeutralFourMomentum_Px{};
	DoubleVector				m_recoPVARecoNeutralFourMomentum_Py{};
	DoubleVector				m_recoPVARecoNeutralFourMomentum_Pz{};
	DoubleVector				m_recoPVARecoNeutralFourMomentum_E{};
	DoubleVector				m_recoPVARecoNeutralFourMomentum_M{};
	DoubleVector				m_usedNeutralFourMomentum_Px{};
	DoubleVector				m_usedNeutralFourMomentum_Py{};
	DoubleVector				m_usedNeutralFourMomentum_Pz{};
	DoubleVector				m_usedNeutralFourMomentum_E{};
	DoubleVector				m_usedNeutralFourMomentum_M{};
	DoubleVector				m_recoPVARecoFourMomentum_minus_truePVARecoNeutralFourMomentum_Px{};
	DoubleVector				m_recoPVARecoFourMomentum_minus_truePVARecoNeutralFourMomentum_Py{};
	DoubleVector				m_recoPVARecoFourMomentum_minus_truePVARecoNeutralFourMomentum_Pz{};
	DoubleVector				m_recoPVARecoFourMomentum_minus_truePVARecoNeutralFourMomentum_E{};
	DoubleVector				m_recoPVARecoFourMomentum_minus_truePVARecoNeutralFourMomentum_M{};
	DoubleVector				m_recoPVARecoFourMomentum_minus_truePVARecoChargedFourMomentum_Px{};
	DoubleVector				m_recoPVARecoFourMomentum_minus_truePVARecoChargedFourMomentum_Py{};
	DoubleVector				m_recoPVARecoFourMomentum_minus_truePVARecoChargedFourMomentum_Pz{};
	DoubleVector				m_recoPVARecoFourMomentum_minus_truePVARecoChargedFourMomentum_E{};
	DoubleVector				m_recoPVARecoFourMomentum_minus_truePVARecoChargedFourMomentum_M{};
	DoubleVector				m_trueNeutrinoFourMomentum_Px{};
	DoubleVector				m_trueNeutrinoFourMomentum_Py{};
	DoubleVector				m_trueNeutrinoFourMomentum_Pz{};
	DoubleVector				m_trueNeutrinoFourMomentum_E{};
	DoubleVector				m_trueNeutrinoFourMomentum_M{};
	DoubleVector				m_recoNeutrinoFourMomentumClose_Px{};
	DoubleVector				m_recoNeutrinoFourMomentumClose_Py{};
	DoubleVector				m_recoNeutrinoFourMomentumClose_Pz{};
	DoubleVector				m_recoNeutrinoFourMomentumClose_E{};
	DoubleVector				m_recoNeutrinoFourMomentumClose_M{};
	DoubleVector				m_trueHadronFourMomentum_Px{};
	DoubleVector				m_trueHadronFourMomentum_Py{};
	DoubleVector				m_trueHadronFourMomentum_Pz{};
	DoubleVector				m_trueHadronFourMomentum_E{};
	DoubleVector				m_trueHadronFourMomentum_M{};
	DoubleVector				m_recoHadronFourMomentum_Px{};
	DoubleVector				m_recoHadronFourMomentum_Py{};
	DoubleVector				m_recoHadronFourMomentum_Pz{};
	DoubleVector				m_recoHadronFourMomentum_E{};
	DoubleVector				m_recoHadronFourMomentum_M{};
	DoubleVector				m_trueFlightDirection_X{};
	DoubleVector				m_trueFlightDirection_Y{};
	DoubleVector				m_trueFlightDirection_Z{};
	DoubleVector				m_recoFlightDirection_X{};
	DoubleVector				m_recoFlightDirection_Y{};
	DoubleVector				m_recoFlightDirection_Z{};
	DoubleVector				m_usedFlightDirection_X{};
	DoubleVector				m_usedFlightDirection_Y{};
	DoubleVector				m_usedFlightDirection_Z{};

	DoubleVector				m_NeutralPVAAlpha{};
	DoubleVector				m_NeutralPVASinAlpha{};
	DoubleVector				m_NeutralPVACosAlpha{};
	DoubleVector				m_ChargedPVAAlpha{};
	DoubleVector				m_ChargedPVASinAlpha{};
	DoubleVector				m_ChargedPVACosAlpha{};
	DoubleVector				m_PVAAlpha{};
	DoubleVector				m_PVASinAlpha{};
	DoubleVector				m_PVACosAlpha{};
	DoubleVector				m_Alpha_RecoPVARecoAll_minus_TruePVARecoCharged_vs_recoPVArecoNeutral{};
	DoubleVector				m_SinAlpha_RecoPVARecoAll_minus_TruePVARecoCharged_vs_recoPVArecoNeutral{};
	DoubleVector				m_CosAlpha_RecoPVARecoAll_minus_TruePVARecoCharged_vs_recoPVArecoNeutral{};
	DoubleVector				m_Alpha_RecoPVARecoAll_minus_TruePVARecoNeutral_vs_recoPVArecoCharged{};
	DoubleVector				m_SinAlpha_RecoPVARecoAll_minus_TruePVARecoNeutral_vs_recoPVArecoCharged{};
	DoubleVector				m_CosAlpha_RecoPVARecoAll_minus_TruePVARecoNeutral_vs_recoPVArecoCharged{};

	DoubleVector				m_visiblePxBeforePVA{};
	DoubleVector				m_visiblePyBeforePVA{};
	DoubleVector				m_visiblePzBeforePVA{};
	DoubleVector				m_visibleEnergyBeforePVA{};
	DoubleVector				m_visibleMassBeforePVA{};
	DoubleVector				m_visiblePxAfterVertexPVA{};
	DoubleVector				m_visiblePyAfterVertexPVA{};
	DoubleVector				m_visiblePzAfterVertexPVA{};
	DoubleVector				m_visibleEnergyAfterVertexPVA{};
	DoubleVector				m_visibleMassAfterVertexPVA{};
	DoubleVector				m_visiblePxAfterChargedPVA{};
	DoubleVector				m_visiblePyAfterChargedPVA{};
	DoubleVector				m_visiblePzAfterChargedPVA{};
	DoubleVector				m_visibleEnergyAfterChargedPVA{};
	DoubleVector				m_visibleMassAfterChargedPVA{};
	DoubleVector				m_visiblePxAfterNeutralPVA{};
	DoubleVector				m_visiblePyAfterNeutralPVA{};
	DoubleVector				m_visiblePzAfterNeutralPVA{};
	DoubleVector				m_visibleEnergyAfterNeutralPVA{};
	DoubleVector				m_visibleMassAfterNeutralPVA{};

	DoubleVector				m_trueVisibleFourMomentumPx{};
	DoubleVector				m_trueVisibleFourMomentumPy{};
	DoubleVector				m_trueVisibleFourMomentumPz{};
	DoubleVector				m_trueVisibleFourMomentumE{};
	DoubleVector				m_recoVisibleFourMomentumPx{};
	DoubleVector				m_recoVisibleFourMomentumPy{};
	DoubleVector				m_recoVisibleFourMomentumPz{};
	DoubleVector				m_recoVisibleFourMomentumE{};
	DoubleVector				m_residualVisibleFourMomentumPx{};
	DoubleVector				m_residualVisibleFourMomentumPy{};
	DoubleVector				m_residualVisibleFourMomentumPz{};
	DoubleVector				m_residualVisibleFourMomentumE{};
	DoubleVector				m_sigmaPxPx_Det{};
	DoubleVector				m_sigmaPxPy_Det{};
	DoubleVector				m_sigmaPyPy_Det{};
	DoubleVector				m_sigmaPxPz_Det{};
	DoubleVector				m_sigmaPyPz_Det{};
	DoubleVector				m_sigmaPzPz_Det{};
	DoubleVector				m_sigmaPxE_Det{};
	DoubleVector				m_sigmaPyE_Det{};
	DoubleVector				m_sigmaPzE_Det{};
	DoubleVector				m_sigmaEE_Det{};
	DoubleVector				m_normalizedResidualVisibleFourMomentumPx{};
	DoubleVector				m_normalizedResidualVisibleFourMomentumPy{};
	DoubleVector				m_normalizedResidualVisibleFourMomentumPz{};
	DoubleVector				m_normalizedResidualVisibleFourMomentumE{};
	DoubleVector				m_trueFlightDirectionUx{};
	DoubleVector				m_trueFlightDirectionUy{};
	DoubleVector				m_trueFlightDirectionUz{};
	DoubleVector				m_recoFlightDirectionUx{};
	DoubleVector				m_recoFlightDirectionUy{};
	DoubleVector				m_recoFlightDirectionUz{};
	DoubleVector				m_residualFlightDirectionUx{};
	DoubleVector				m_residualFlightDirectionUy{};
	DoubleVector				m_residualFlightDirectionUz{};
	DoubleVector				m_sigmaUxUx{};
	DoubleVector				m_sigmaUxUy{};
	DoubleVector				m_sigmaUyUy{};
	DoubleVector				m_sigmaUxUz{};
	DoubleVector				m_sigmaUyUz{};
	DoubleVector				m_sigmaUzUz{};
	DoubleVector				m_normalizedResidualFlightDirectionUx{};
	DoubleVector				m_normalizedResidualFlightDirectionUy{};
	DoubleVector				m_normalizedResidualFlightDirectionUz{};
	DoubleVector				m_trueNeutrinoFourMomentumPx{};
	DoubleVector				m_trueNeutrinoFourMomentumPy{};
	DoubleVector				m_trueNeutrinoFourMomentumPz{};
	DoubleVector				m_trueNeutrinoFourMomentumE{};
	DoubleVector				m_recoNeutrinoFourMomentumClosePx{};
	DoubleVector				m_recoNeutrinoFourMomentumClosePy{};
	DoubleVector				m_recoNeutrinoFourMomentumClosePz{};
	DoubleVector				m_recoNeutrinoFourMomentumCloseE{};
	DoubleVector				m_sigmaNeutrinoPxPx{};
	DoubleVector				m_sigmaNeutrinoPxPy{};
	DoubleVector				m_sigmaNeutrinoPyPy{};
	DoubleVector				m_sigmaNeutrinoPxPz{};
	DoubleVector				m_sigmaNeutrinoPyPz{};
	DoubleVector				m_sigmaNeutrinoPzPz{};
	DoubleVector				m_sigmaNeutrinoPxE{};
	DoubleVector				m_sigmaNeutrinoPyE{};
	DoubleVector				m_sigmaNeutrinoPzE{};
	DoubleVector				m_sigmaNeutrinoEE{};
	DoubleVector				m_residualNeutrinoFourMomentumPx{};
	DoubleVector				m_residualNeutrinoFourMomentumPy{};
	DoubleVector				m_residualNeutrinoFourMomentumPz{};
	DoubleVector				m_residualNeutrinoFourMomentumE{};
	DoubleVector				m_normalizedResidualNeutrinoFourMomentumPx{};
	DoubleVector				m_normalizedResidualNeutrinoFourMomentumPy{};
	DoubleVector				m_normalizedResidualNeutrinoFourMomentumPz{};
	DoubleVector				m_normalizedResidualNeutrinoFourMomentumE{};
	DoubleVector				m_recoNeutrinoDirectionError{};
	DoubleVector				m_PCAatLeptonX{};
	DoubleVector				m_PCAatLeptonY{};
	DoubleVector				m_PCAatLeptonZ{};
	DoubleVector				m_PCAatOtherParticleX{};
	DoubleVector				m_PCAatOtherParticleY{};
	DoubleVector				m_PCAatOtherParticleZ{};
	DoubleVector				m_JetAxisX{};
	DoubleVector				m_JetAxisY{};
	DoubleVector				m_JetAxisZ{};
	DoubleVector				m_PrimaryVertexX{};
	DoubleVector				m_PrimaryVertexY{};
	DoubleVector				m_PrimaryVertexZ{};
	DoubleVector				m_AngleTrueFlightDirectionJet{};
	DoubleVector				m_AngleRecoFlightDirectionJet{};
	DoubleVector				m_AngleLeptonJet{};
	DoubleVector				m_AngleDSVertexJet{};
	DoubleVector				m_AngleRecoNeutrinoJet{};
	DoubleVector				m_AngleTrueNeutrinoJet{};
	DoubleVector				m_LeptonDistanceFromPV{};
	DoubleVector				m_DSVDistanceFromPV{};
	DoubleVector				m_Lepton3DImpactParameter{};
	DoubleVector				m_OtherParticle3DImpactParameter{};
	IntVector				m_SLDStatus{};
	IntVector				m_PVAStatus{};
	IntVector				m_trueSolutionSign{};
	IntVector				m_recoSolutionSign{};
	IntVector				m_nChargedParticlesInPrimVertex{};
	IntVector				m_nChargedParticlesNotInPrimVertex{};
	FloatVector				m_weightPFOtoMCP_Lepton{};
	FloatVector				m_weightMCPtoPFO_Lepton{};
	FloatVector				m_weightPFOtoMCP_Neutral{};
	FloatVector				m_weightMCPtoPFO_Neutral{};
	FloatVector				m_weightPFOtoMCP_Charged{};
	FloatVector				m_weightMCPtoPFO_Charged{};
	TH1I					*h_SLDStatus{};
	TH1F					*h_BHadronType{};
	TH1F					*h_CHadronType{};
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
	TH1I					*h_SLDecayFlavour{};
	TH1I					*h_SLDecayModeB{};
	TH1I					*h_SLDecayModeC{};
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
	TTree					*m_pTTree3{};


};
#endif
