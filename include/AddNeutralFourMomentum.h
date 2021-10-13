#ifndef AddNeutralFourMomentum_h_1
#define AddNeutralFourMomentum_h_1


#include "lcio.h"
#include "streamlog/streamlog.h"
#include <EVENT/ReconstructedParticle.h>
#include "TVector3.h"
#include "TLorentzVector.h"

int getNeutralFourMomentum( EVENT::ReconstructedParticle *assignedJet , TVector3 pointingVector , double InvMass , TLorentzVector chargedFourMomentum , TLorentzVector &neutralFourMomentum );
std::vector<EVENT::ReconstructedParticle*> sortParticles( std::vector<EVENT::ReconstructedParticle*> sortedParticles , std::vector<EVENT::ReconstructedParticle*> remainedParticles , TVector3 pointingVector );

#endif
