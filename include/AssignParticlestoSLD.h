#ifndef AssignParticlestoSLD_h_1
#define AssignParticlestoSLD_h_1


#include "lcio.h"
#include "streamlog/streamlog.h"
#include "marlin/Global.h"
#include "marlin/Processor.h"
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/Track.h>
#include <EVENT/LCCollection.h>
#include <EVENT/Vertex.h>
#include <IMPL/TrackImpl.h>
#include "TVector3.h"
#include "TLorentzVector.h"
#include "DDMarlinCED.h"
#include <math.h>
#include <GeometryUtil.h>

typedef std::vector<EVENT::ReconstructedParticle*>	pfoVector;

void assignParticlesToSemiLeptonicDecay( pfoVector &assignedParticles , pfoVector &availableParticles , double invariantMass , TVector3 direction );

bool isParticleInVertex( EVENT::ReconstructedParticle *particle , EVENT::Vertex *vertex );

std::vector<EVENT::ReconstructedParticle*> getParticlesWithAloneTracks( EVENT::ReconstructedParticle *particle , EVENT::ReconstructedParticle *jet , EVENT::Vertex *primaryVertex , std::vector<EVENT::Vertex*> vertexVector );

EVENT::Vertex *getParticleVertex( EVENT::ReconstructedParticle *particle , std::vector<EVENT::Vertex*> vertexVector , bool &foundParticleInVertex );

EVENT::ReconstructedParticle *getJetAssignedToParticle( EVENT::ReconstructedParticle *particle , std::vector<EVENT::ReconstructedParticle*> jetVector , bool &foundParticleInJet );

EVENT::Track *makeTrackParameters( TVector3 momentum , TVector3 point , double charge );

void sortParticles( pfoVector &sortedParticles , pfoVector &unSortedParticles , TVector3 direction );

std::vector<EVENT::Vertex*> getVerticesInJet( EVENT::ReconstructedParticle * assignedJet , std::vector<EVENT::Vertex*> buildUpVertexVector );

EVENT::ReconstructedParticle* getLeadingChargedParticle( EVENT::ReconstructedParticle* Jet );

EVENT::ReconstructedParticle* getLeadingNeutralParticle( EVENT::ReconstructedParticle* Jet );

#endif
