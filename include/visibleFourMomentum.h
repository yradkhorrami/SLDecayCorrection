#ifndef visibleFourMomentum_h_1
#define visibleFourMomentum_h_1


#include "lcio.h"
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include "UTIL/LCRelationNavigator.h"
#include "TLorentzVector.h"
#include <string>

typedef std::vector<EVENT::ReconstructedParticle*> pfoVector;

typedef std::vector<EVENT::MCParticle*> mcpVector;

//EVENT::ReconstructedParticle* getLinkedPFO( EVENT::LCEvent *pLCEvent , EVENT::MCParticle *visibleMCP , std::string recoMCTruthLinkCollection , std::string mcTruthRecoLinkCollection , bool getChargedTLV , bool getNeutralTLV );

int getLeptonFourMomentum( EVENT::LCEvent *pLCEvent , EVENT::MCParticle *SLDLepton  , bool cheatLepton4momentum , TLorentzVector &fourMomentumLepton , TLorentzVector &trueFourMomentumLepton , std::string recoMCTruthLinkCollection , std::string mcTruthRecoLinkCollection );

int getChargedFourMomentum( EVENT::LCEvent *pLCEvent , EVENT::MCParticle *SLDLepton  , bool cheatCharged4momentum , TLorentzVector &fourMomentumCharged , TLorentzVector &trueFourMomentumCharged , float &restCharge , std::string recoMCTruthLinkCollection , std::string mcTruthRecoLinkCollection );

int getLinkedChargedPFOs( EVENT::LCEvent *pLCEvent , EVENT::MCParticle *SLDLepton , EVENT::MCParticle *parentHadron , pfoVector &chargedPFOs , std::string recoMCTruthLinkCollection , std::string mcTruthRecoLinkCollection );

int getLinkedChargedMCPs( EVENT::LCEvent *pLCEvent , EVENT::MCParticle *SLDLepton , EVENT::MCParticle *parentHadron , mcpVector &chargedMCPs );

int getNeutralFourMomentum( EVENT::LCEvent *pLCEvent , EVENT::MCParticle *SLDLepton  , bool cheatNeutral4momentum , TLorentzVector &fourMomentumNeutral , TLorentzVector &trueFourMomentumNeutral , std::string recoMCTruthLinkCollection , std::string mcTruthRecoLinkCollection );

int getLinkedNeutralPFOs( EVENT::LCEvent *pLCEvent , EVENT::MCParticle *SLDLepton , EVENT::MCParticle *parentHadron , pfoVector &neutralPFOs , std::string recoMCTruthLinkCollection , std::string mcTruthRecoLinkCollection );

int getLinkedNeutralMCPs( EVENT::LCEvent *pLCEvent , EVENT::MCParticle *SLDLepton , EVENT::MCParticle *parentHadron , mcpVector &neutralMCPs );

//int getCovMatrixMomenta(EVENT::ReconstructedParticle const *, EVENT::FloatVec &);


#endif
