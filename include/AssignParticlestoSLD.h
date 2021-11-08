#ifndef AssignParticlestoSLD_h_1
#define AssignParticlestoSLD_h_1


#include "lcio.h"
#include "streamlog/streamlog.h"
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/Track.h>
#include "EVENT/LCCollection.h"
#include <EVENT/Vertex.h>
#include "TVector3.h"
#include "TLorentzVector.h"

std::vector<EVENT::ReconstructedParticle*> assignNeutralParticles( EVENT::ReconstructedParticle *assignedJet , TVector3 pointingVector , double InvMass , TLorentzVector chargedFourMomentum , TLorentzVector &neutralFourMomentum , double cosOpeningAngle );

std::vector<EVENT::ReconstructedParticle*> assignChargedParticles( EVENT::LCEvent *pLCEvent , EVENT::ReconstructedParticle *assignedJet , std::string inputPrimaryVertex , std::string inputBuildUpVertex , EVENT::ReconstructedParticle *linkedRecoLepton , TVector3 pointingVector , double InvMass , TLorentzVector chargedFourMomentum , TLorentzVector &neutralFourMomentum , double cosOpeningAngle );

std::vector<EVENT::ReconstructedParticle*> sortParticlesInCone( std::vector<EVENT::ReconstructedParticle*> sortedParticles , std::vector<EVENT::ReconstructedParticle*> remainedParticles , TVector3 pointingVector );

#endif
