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
	typedef std::vector<EVENT::Vertex*>			vtxVector;
	bool hasPrimarySLDecay( MCParticle *parentHadron );
	bool hasDownStreamSLDecay( MCParticle *parentHadron );
	bool hasUpStreamSLDecay( MCParticle *parentHadron );
	bool checkBHadronSLDecay( MCParticle *SLDLepton );
	bool checkCHadronSLDecay( MCParticle *SLDLepton );
	bool checkTauLeptonSLDecay( MCParticle *SLDLepton );
	virtual void doSLDCorrection( EVENT::LCEvent *pLCEvent , MCParticle *SLDLepton , vtxVector& semiLeptonicVertices , pfoVector& semiLeptonicVertexRecoParticles , pfoVector& jetsOfSemiLeptonicDecays , pfoVector& neutrinos );
	void showTrueParameters( MCParticle *SLDLepton );
	TLorentzVector getNeutrinoFourMomentum( TVector3 flightDirection , TLorentzVector visibleFourMomentum , double parentHadronMass , float solutionSign );
	TLorentzVector getNeutrinoFourMomentumModified( TVector3 flightDirection , TLorentzVector visibleFourMomentum , double parentHadronMass , float solutionSign );
	TLorentzVector getNeutrinoFourMomentumStandardMethod( TVector3 flightDirection , TLorentzVector visibleFourMomentum , double parentHadronMass , float solutionSign );
	MCParticle* getTrueNeutrino( MCParticle *SLDLepton );
	void fillTrueRecoFourMomentum( TLorentzVector trueNeutralFourMomentum , TLorentzVector trueChargedFourMomentum , TLorentzVector trueLeptonFourMomentum , TLorentzVector trueVisibleFourMomentum , TLorentzVector trueNeutrinoFourMomentum , TLorentzVector trueHadronFourMomentum ,  TLorentzVector recoNeutralFourMomentum , TLorentzVector recoChargedFourMomentum , TLorentzVector recoLeptonFourMomentum , TLorentzVector recoVisibleFourMomentum , TLorentzVector recoNeutrinoFourMomentum , TLorentzVector recoHadronFourMomentum , TLorentzVector usedNeutralFourMomentum , TLorentzVector usedChargedFourMomentum , TLorentzVector usedLeptonFourMomentum , TLorentzVector usedVisibleFourMomentum , TLorentzVector cheatedPVARecoChargedFourMomentum , TLorentzVector cheatedPVARecoNeutralFourMomentum );
	virtual void getCovMatPVA( std::vector<EVENT::ReconstructedParticle*> decayProducts , std::vector<EVENT::ReconstructedParticle*> associatedParticles , int SLDStatus , std::vector< float > &CovMatrixPVA );
	virtual void getCovMatDetFlightDirection( std::vector<EVENT::ReconstructedParticle*> decayProducts , TVector3 flightDirection , std::vector< float > CovMatrixFlightDirection , std::vector< float > &CovMatrixDetector , EVENT::ReconstructedParticle* linkedRecoLepton , std::vector< float > &CovMatrixDetPar , std::vector< float > &CovMatrixDetNor );
	virtual void getCovMatFlightDirection( TVector3 flightDirection , float sigmaTheta , float sigmaPhi , std::vector< float > &CovMatrixFlightDirection );
	virtual void getCovMatrixDetPar( TVector3 flightDirection , TLorentzVector visibleFourMomentum , std::vector< float > CovMatrixFlightDirection , std::vector< float > initialCovMatrixDetector , std::vector< float > &CovMatrixDetPar );
	virtual void getCovMatrixDetNor( TVector3 flightDirection , TLorentzVector visibleFourMomentum , std::vector< float > CovMatrixFlightDirection , std::vector< float > initialCovMatrixDetector , std::vector< float > &CovMatrixDetNor );
	virtual void getNeutrinoCovMat( TLorentzVector recoNeutrinoFourMomentum , TLorentzVector visibleFourMomentum , TVector3 flightDirection , double parentHadronMass , std::vector< float > CovMatrixPVA , std::vector< float > CovMatrixDetector , std::vector< float > CovMatrixDetPar , std::vector< float > CovMatrixDetNor , std::vector< float > &NeutrinoCovMatrix );
	virtual void plotHistograms( TLorentzVector trueFourMomentumNeutrino , TLorentzVector FourMomentumNuClose , std::vector<float> NeutrinoCovMat );
	virtual void InitializeHistogram( TH1F *histogram , int scale , int color , int lineWidth , int markerSize , int markerStyle );
	virtual void doProperGaussianFit( TH1F *histogram , float fitMin , float fitMax , float fitRange );
	int getVertexInJetsDistribution( Vertex* testVertex , std::vector<EVENT::ReconstructedParticle*> jetVector );
	void evaluatePFOsAngle( std::vector<EVENT::ReconstructedParticle*> aloneChargedPFOsInJetFromSLD , std::vector<EVENT::ReconstructedParticle*> aloneChargedPFOsInJetNotFromSLD , std::vector<EVENT::ReconstructedParticle*> chargedPFOsInJetFromSLD , std::vector<EVENT::ReconstructedParticle*> chargedPFOsInJetNotFromSLD , std::vector<EVENT::ReconstructedParticle*> neutralPFOsInJetFromSLD , std::vector<EVENT::ReconstructedParticle*> neutralPFOsInJetNotFromSLD , TVector3 leptonDirection , TVector3 jetAxis , TVector3 recoFlightDirection , int SLDStatus );
	void investigateJetEnergyContent( EVENT::ReconstructedParticle *assignedJet );

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
	int					m_vertexingScenario = 1;
	bool					m_cheatLepton4momentum = true;
	bool					m_cheatCharged4momentum = true;
	bool					m_cheatNeutral4momentum = true;
	bool					m_cheatPVAcharged = true;
	float					m_chargedCosAcceptanceAngleSLD4 = 0.0;
	float					m_chargedCosAcceptanceAngleSLD5 = 0.0;
	bool					m_cheatPVAneutral = true;
	float					m_neutralCosAcceptanceAngle = 0.0;
	float					m_BSLDChargedSLD4InvMassCut = 0.0;
	float					m_BSLDChargedSLD5InvMassCut = 0.0;
	float					m_BSLDNeutralSLD4InvMassCut = 0.0;
	float					m_BSLDNeutralSLD5InvMassCut = 0.0;
	float					m_CSLDChargedSLD4InvMassCut = 0.0;
	float					m_CSLDChargedSLD5InvMassCut = 0.0;
	float					m_CSLDNeutralSLD4InvMassCut = 0.0;
	float					m_CSLDNeutralSLD5InvMassCut = 0.0;
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
	DoubleVector				m_P_vis_nor{};
	DoubleVector				m_P_vis_nor_prime{};
	IntVector				m_flightDirectionStatus{};
	DoubleVector				m_FlightDirectionErrorSinAlpha{};
	DoubleVector				m_FlightDirectionErrorCosAlpha{};
	DoubleVector				m_FlightDirectionErrorAlpha{};
	DoubleVector				m_FlightDirectionErrorDeltaTheta{};
	DoubleVector				m_FlightDirectionErrorDeltaPhi{};
	DoubleVector				m_distRecoLeptonToDownStreamVertex{};
	DoubleVector				m_dsVertexResidualX{};
	DoubleVector				m_dsVertexResidualY{};
	DoubleVector				m_dsVertexResidualZ{};
	DoubleVector				m_SecVertexResidualX{};
	DoubleVector				m_SecVertexResidualY{};
	DoubleVector				m_SecVertexResidualZ{};
	DoubleVector				m_parentHadronMass{};
	IntVector				m_parentHadronPDG{};
	DoubleVector				m_trueParentHadronFlightDistance{};
	DoubleVector				m_recoParentHadronFlightDistance{};
	DoubleVector				m_daughterHadronMass{};
	IntVector				m_daughterHadronPDG{};
	DoubleVector				m_daughterHadronFlightDistance{};
	FloatVector				m_visibleChargedInvMassCut{};
	FloatVector				m_visibleNeutralInvMassCut{};
	DoubleVector				m_alphaTrueNeutrals{};
	DoubleVector				m_alphaTrueCharged{};
	DoubleVector				m_cosAlphaNeutrals{};
	DoubleVector				m_cosAlphaCharged{};
	DoubleVector				m_alphaChargedPFOsFromSLDwrtLepton{};
	DoubleVector				m_cosAlphaChargedPFOsFromSLDwrtLepton{};
	DoubleVector				m_alphaChargedPFOsFromSLDwrtFD{};
	DoubleVector				m_cosAlphaChargedPFOsFromSLDwrtFD{};
	DoubleVector				m_alphaChargedPFOsFromSLDwrtJet{};
	DoubleVector				m_cosAlphaChargedPFOsFromSLDwrtJet{};
	DoubleVector				m_alphaChargedPFOsNotFromSLDwrtLepton{};
	DoubleVector				m_cosAlphaChargedPFOsNotFromSLDwrtLepton{};
	DoubleVector				m_alphaChargedPFOsNotFromSLDwrtFD{};
	DoubleVector				m_cosAlphaChargedPFOsNotFromSLDwrtFD{};
	DoubleVector				m_alphaChargedPFOsNotFromSLDwrtJet{};
	DoubleVector				m_cosAlphaChargedPFOsNotFromSLDwrtJet{};

	DoubleVector				m_alphaAloneChargedPFOsFromSLDwrtLepton{};
	DoubleVector				m_cosAlphaAloneChargedPFOsFromSLDwrtLepton{};
	DoubleVector				m_alphaAloneChargedPFOsFromSLDwrtFD{};
	DoubleVector				m_cosAlphaAloneChargedPFOsFromSLDwrtFD{};
	DoubleVector				m_alphaAloneChargedPFOsFromSLDwrtJet{};
	DoubleVector				m_cosAlphaAloneChargedPFOsFromSLDwrtJet{};
	DoubleVector				m_alphaAloneChargedPFOsNotFromSLDwrtLepton{};
	DoubleVector				m_cosAlphaAloneChargedPFOsNotFromSLDwrtLepton{};
	DoubleVector				m_alphaAloneChargedPFOsNotFromSLDwrtFD{};
	DoubleVector				m_cosAlphaAloneChargedPFOsNotFromSLDwrtFD{};
	DoubleVector				m_alphaAloneChargedPFOsNotFromSLDwrtJet{};
	DoubleVector				m_cosAlphaAloneChargedPFOsNotFromSLDwrtJet{};
	DoubleVector				m_alphaNeutralPFOsFromSLDwrtLepton{};
	DoubleVector				m_cosAlphaNeutralPFOsFromSLDwrtLepton{};
	DoubleVector				m_alphaNeutralPFOsFromSLDwrtFD{};
	DoubleVector				m_cosAlphaNeutralPFOsFromSLDwrtFD{};
	DoubleVector				m_alphaNeutralPFOsFromSLDwrtJet{};
	DoubleVector				m_cosAlphaNeutralPFOsFromSLDwrtJet{};
	DoubleVector				m_alphaNeutralPFOsNotFromSLDwrtLepton{};
	DoubleVector				m_cosAlphaNeutralPFOsNotFromSLDwrtLepton{};
	DoubleVector				m_alphaNeutralPFOsNotFromSLDwrtFD{};
	DoubleVector				m_cosAlphaNeutralPFOsNotFromSLDwrtFD{};
	DoubleVector				m_alphaNeutralPFOsNotFromSLDwrtJet{};
	DoubleVector				m_cosAlphaNeutralPFOsNotFromSLDwrtJet{};

	DoubleVector				m_trueNeutralPx{};
	DoubleVector				m_trueNeutralPy{};
	DoubleVector				m_trueNeutralPz{};
	DoubleVector				m_trueNeutralE{};
	DoubleVector				m_trueNeutralM{};
	DoubleVector				m_trueChargedPx{};
	DoubleVector				m_trueChargedPy{};
	DoubleVector				m_trueChargedPz{};
	DoubleVector				m_trueChargedE{};
	DoubleVector				m_trueChargedM{};
	DoubleVector				m_trueLeptonPx{};
	DoubleVector				m_trueLeptonPy{};
	DoubleVector				m_trueLeptonPz{};
	DoubleVector				m_trueLeptonE{};
	DoubleVector				m_trueVisiblePx{};
	DoubleVector				m_trueVisiblePy{};
	DoubleVector				m_trueVisiblePz{};
	DoubleVector				m_trueVisibleE{};
	DoubleVector				m_trueVisibleM{};
	DoubleVector				m_trueNeutrinoPx{};
	DoubleVector				m_trueNeutrinoPy{};
	DoubleVector				m_trueNeutrinoPz{};
	DoubleVector				m_trueNeutrinoE{};
	DoubleVector				m_trueHadronPx{};
	DoubleVector				m_trueHadronPy{};
	DoubleVector				m_trueHadronPz{};
	DoubleVector				m_trueHadronE{};
	DoubleVector				m_trueHadronM{};
	DoubleVector				m_recoNeutralPx{};
	DoubleVector				m_recoNeutralPy{};
	DoubleVector				m_recoNeutralPz{};
	DoubleVector				m_recoNeutralE{};
	DoubleVector				m_recoNeutralM{};
	DoubleVector				m_recoChargedPx{};
	DoubleVector				m_recoChargedPy{};
	DoubleVector				m_recoChargedPz{};
	DoubleVector				m_recoChargedE{};
	DoubleVector				m_recoChargedM{};
	DoubleVector				m_recoLeptonPx{};
	DoubleVector				m_recoLeptonPy{};
	DoubleVector				m_recoLeptonPz{};
	DoubleVector				m_recoLeptonE{};
	DoubleVector				m_recoVisiblePx{};
	DoubleVector				m_recoVisiblePy{};
	DoubleVector				m_recoVisiblePz{};
	DoubleVector				m_recoVisibleE{};
	DoubleVector				m_recoVisibleM{};
	DoubleVector				m_recoNeutrinoPx{};
	DoubleVector				m_recoNeutrinoPy{};
	DoubleVector				m_recoNeutrinoPz{};
	DoubleVector				m_recoNeutrinoE{};
	DoubleVector				m_recoHadronPx{};
	DoubleVector				m_recoHadronPy{};
	DoubleVector				m_recoHadronPz{};
	DoubleVector				m_recoHadronE{};
	DoubleVector				m_recoHadronM{};
	DoubleVector				m_usedNeutralPx{};
	DoubleVector				m_usedNeutralPy{};
	DoubleVector				m_usedNeutralPz{};
	DoubleVector				m_usedNeutralE{};
	DoubleVector				m_usedNeutralM{};
	DoubleVector				m_usedChargedPx{};
	DoubleVector				m_usedChargedPy{};
	DoubleVector				m_usedChargedPz{};
	DoubleVector				m_usedChargedE{};
	DoubleVector				m_usedChargedM{};
	DoubleVector				m_usedLeptonPx{};
	DoubleVector				m_usedLeptonPy{};
	DoubleVector				m_usedLeptonPz{};
	DoubleVector				m_usedLeptonE{};
	DoubleVector				m_usedVisiblePx{};
	DoubleVector				m_usedVisiblePy{};
	DoubleVector				m_usedVisiblePz{};
	DoubleVector				m_usedVisibleE{};
	DoubleVector				m_usedVisibleM{};
	DoubleVector				m_cheatedPVARecoNeutralPx{};
	DoubleVector				m_cheatedPVARecoNeutralPy{};
	DoubleVector				m_cheatedPVARecoNeutralPz{};
	DoubleVector				m_cheatedPVARecoNeutralE{};
	DoubleVector				m_cheatedPVARecoNeutralM{};
	DoubleVector				m_cheatedPVARecoChargedPx{};
	DoubleVector				m_cheatedPVARecoChargedPy{};
	DoubleVector				m_cheatedPVARecoChargedPz{};
	DoubleVector				m_cheatedPVARecoChargedE{};
	DoubleVector				m_cheatedPVARecoChargedM{};
	DoubleVector				m_neutralEnergy{};
	DoubleVector				m_neutralEnergyFromVertexing{};
	DoubleVector				m_neutralEnergyFromPVA{};
	DoubleVector				m_neutralMomentum{};
	DoubleVector				m_neutralMomentumFromVertexing{};
	DoubleVector				m_neutralMomentumFromPVA{};
	DoubleVector				m_chargedEnergy{};
	DoubleVector				m_chargedEnergyFromVertexing{};
	DoubleVector				m_chargedEnergyFromPVA{};
	DoubleVector				m_chargedMomentum{};
	DoubleVector				m_chargedMomentumFromVertexing{};
	DoubleVector				m_chargedMomentumFromPVA{};
	DoubleVector				m_expectedNeutralEnergy{};
	DoubleVector				m_expectedNeutralEnergyFromVertexing{};
	DoubleVector				m_expectedNeutralEnergyFromPVA{};
	DoubleVector				m_expectedNeutralMomentum{};
	DoubleVector				m_expectedNeutralMomentumFromVertexing{};
	DoubleVector				m_expectedNeutralMomentumFromPVA{};
	DoubleVector				m_expectedChargedEnergy{};
	DoubleVector				m_expectedChargedEnergyFromVertexing{};
	DoubleVector				m_expectedChargedEnergyFromPVA{};
	DoubleVector				m_expectedChargedMomentum{};
	DoubleVector				m_expectedChargedMomentumFromVertexing{};
	DoubleVector				m_expectedChargedMomentumFromPVA{};


	IntVector				m_SLDStatus{};
	FloatVector				m_weightPFOtoMCP_Lepton{};
	FloatVector				m_weightMCPtoPFO_Lepton{};
	FloatVector				m_weightPFOtoMCP_Neutral{};
	FloatVector				m_weightMCPtoPFO_Neutral{};
	FloatVector				m_weightPFOtoMCP_Charged{};
	FloatVector				m_weightMCPtoPFO_Charged{};
	TH1F					*h_SLDStatus{};
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
