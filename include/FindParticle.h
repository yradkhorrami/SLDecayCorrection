#ifndef FindParticle_h_1
#define FindParticle_h_1


#include "lcio.h"
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include "UTIL/LCRelationNavigator.h"
#include "DDMarlinCED.h"
#include "TVector3.h"

typedef EVENT::ReconstructedParticle* pfo;

typedef std::vector<EVENT::ReconstructedParticle*> pfoVector;

typedef std::vector<EVENT::MCParticle*> mcpVector;

void getTrueDecayProducts( EVENT::MCParticle *parentMCP , EVENT::MCParticle *SLDLepton , EVENT::MCParticle *trueNeutrino , mcpVector &trueNeutralDecayProducts , mcpVector &trueChargedDecayProducts );

void getTruePVADecayProducts( EVENT::MCParticle *parentMCP , EVENT::MCParticle *SLDLepton , EVENT::MCParticle *trueNeutrino , pfo &linkedRecoLepton , float &weightRecoLeptoMCLep , float &weightMCLeptoRecoLep , pfoVector &truePVANeutralDecayProducts , pfoVector &truePVAChargedDecayProducts , LCRelationNavigator RecoMCParticleNav , LCRelationNavigator MCParticleRecoNav );

void getLinkedRecoLepton( EVENT::MCParticle *mcParticle , pfo &linkedRecoLepton , LCRelationNavigator RecoMCParticleNav , LCRelationNavigator MCParticleRecoNav , float &weightRecoLeptoMCLep , float &weightMCLeptoRecoLep );

void getDecayProducts( mcpVector trueDecayProducts , mcpVector &trueNeutralDecayProducts , mcpVector &trueChargedDecayProducts , pfoVector &cheatedPVANeutralDecayProducts , pfoVector &cheatedPVAChargedDecayProducts , LCRelationNavigator RecoMCParticleNav , LCRelationNavigator MCParticleRecoNav );

void investigateDecayChain( EVENT::MCParticle *parentMCP , mcpVector &trueNeutralDecayProducts , mcpVector &trueChargedDecayProducts , pfoVector &cheatedPVANeutralDecayProducts , pfoVector &cheatedPVAChargedDecayProducts , LCRelationNavigator RecoMCParticleNav , LCRelationNavigator MCParticleRecoNav );

void getDecayProducts( EVENT::MCParticle *parentMCP , EVENT::MCParticle *SLDLepton , EVENT::MCParticle *trueNeutrino , mcpVector &trueNeutralDecayProducts , mcpVector &trueChargedDecayProducts , pfo &linkedRecoLepton , pfoVector &cheatedPVANeutralDecayProducts , pfoVector &cheatedPVAChargedDecayProducts , LCRelationNavigator RecoMCParticleNav , LCRelationNavigator MCParticleRecoNav , float &weightRecoLep2MCLep , float &weightMCLep2RecoLep , int generatorStatus );

EVENT::ReconstructedParticle* getLinkedRecoLepton( EVENT::MCParticle *mcParticle , LCRelationNavigator RecoMCParticleNav , LCRelationNavigator MCParticleRecoNav , float &weightPFOtoMCP , float &weightMCPtoPFO , pfoVector &cheatedPVANeutralDecayProducts );

int getNeutralMCPs( EVENT::MCParticle *parentMCP , mcpVector &neutralMCPs );

int getChargedMCPs( EVENT::MCParticle *SLDLepton , EVENT::MCParticle *parentMCP , mcpVector &chargedMCPs );

int getTrueVertices( EVENT::MCParticle *parentMCP , mcpVector &MCPsWithVertex );

EVENT::MCParticle* getLinkedMCP( EVENT::ReconstructedParticle *recoParticle , LCRelationNavigator RecoMCParticleNav , LCRelationNavigator MCParticleRecoNav , bool getChargedMCP , bool getNeutralMCP , float &weightPFOtoMCP , float &weightMCPtoPFO );

EVENT::ReconstructedParticle* getLinkedPFO( EVENT::MCParticle *mcParticle , LCRelationNavigator RecoMCParticleNav , LCRelationNavigator MCParticleRecoNav , bool getChargedPFO , bool getNeutralPFO , float &weightPFOtoMCP , float &weightMCPtoPFO );

EVENT::ReconstructedParticle* getLinkedPFO( EVENT::MCParticle *mcParticle , LCRelationNavigator RecoMCParticleNav , LCRelationNavigator MCParticleRecoNav , bool getChargedPFO , bool getNeutralPFO , float &weightPFOtoMCP , float &weightMCPtoPFO , bool &foundPFO );

float getWidestCosAlphaOfDecayProducts( std::vector<EVENT::MCParticle*> mcpVec , TVector3 direction );

float getWidestCosAlphaOfVertices( std::vector<EVENT::MCParticle*> mcpVec , TVector3 direction , std::vector<float> truePrimaryVertex );

float getWidestCosAlphaOfChargedPFOs( std::vector<EVENT::ReconstructedParticle*> pfoVec , TVector3 direction );

void isMCParticleFromSLD( EVENT::MCParticle* parentHadron , EVENT::MCParticle* testMCParticle , bool &MCPisFromSLD );

#endif
