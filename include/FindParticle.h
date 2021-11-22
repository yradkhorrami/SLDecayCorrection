#ifndef FindParticle_h_1
#define FindParticle_h_1


#include "lcio.h"
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include "UTIL/LCRelationNavigator.h"
#include "DDMarlinCED.h"
#include "TVector3.h"

typedef std::vector<EVENT::ReconstructedParticle*> pfoVector;

typedef std::vector<EVENT::MCParticle*> mcpVector;

int getNeutralMCPs( EVENT::MCParticle *parentMCP , mcpVector &neutralMCPs );

int getChargedMCPs( EVENT::MCParticle *SLDLepton , EVENT::MCParticle *parentMCP , mcpVector &chargedMCPs );

int getTrueVertices( EVENT::MCParticle *parentMCP , mcpVector &MCPsWithVertex );

EVENT::MCParticle* getLinkedMCP( EVENT::ReconstructedParticle *recoParticle , LCRelationNavigator RecoMCParticleNav , LCRelationNavigator MCParticleRecoNav , bool getChargedMCP , bool getNeutralMCP , float &weightPFOtoMCP , float &weightMCPtoPFO );

EVENT::ReconstructedParticle* getLinkedPFO( EVENT::MCParticle *mcParticle , LCRelationNavigator RecoMCParticleNav , LCRelationNavigator MCParticleRecoNav , bool getChargedPFO , bool getNeutralPFO , float &weightPFOtoMCP , float &weightMCPtoPFO );

float getWidestCosAlphaOfDecayProducts( std::vector<EVENT::MCParticle*> mcpVec , TVector3 direction );

float getWidestCosAlphaOfVertices( std::vector<EVENT::MCParticle*> mcpVec , TVector3 direction , std::vector<float> truePrimaryVertex );

float getWidestCosAlphaOfChargedPFOs( std::vector<EVENT::ReconstructedParticle*> pfoVec , TVector3 direction );

void isMCParticleFromSLD( EVENT::MCParticle* parentHadron , EVENT::MCParticle* testMCParticle , bool &MCPisFromSLD );

#endif
