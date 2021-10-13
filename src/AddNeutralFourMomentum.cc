#include "AddNeutralFourMomentum.h"

int getNeutralFourMomentum( EVENT::ReconstructedParticle *assignedJet , TVector3 pointingVector , double InvMass , TLorentzVector chargedFourMomentum , TLorentzVector &neutralFourMomentum )
{
	pointingVector.SetMag( 1.0 );
	int nParticles = ( assignedJet->getParticles() ).size();
	std::vector<EVENT::ReconstructedParticle*> sortedParticles;
	std::vector<EVENT::ReconstructedParticle*> remainedParticles;
	for ( int i_particle = 0 ; i_particle < nParticles ; ++i_particle )
	{
		EVENT::ReconstructedParticle* particle = assignedJet->getParticles()[ i_particle ];
		remainedParticles.push_back( particle );
		streamlog_out(DEBUG1) << "		adding particle " << particle << " to Particles Vector" << std::endl;
	}
	sortedParticles = sortParticles( sortedParticles , remainedParticles , pointingVector );
	streamlog_out(DEBUG1) << "		" << sortedParticles.size() << " particles sorted wrt the angular distance to pointing vector" << std::endl;

	TLorentzVector totalFourMomentum = chargedFourMomentum + neutralFourMomentum;
	int i_par = 0;
	int addedParticles = 0;
	double InvariantMass = chargedFourMomentum.M();
	while ( InvariantMass <= InvMass )
	{
		EVENT::ReconstructedParticle* particle = assignedJet->getParticles()[ i_par ];
		if ( particle->getCharge() == 0.0 )
		{
			neutralFourMomentum += TLorentzVector( particle->getMomentum()[ 0 ] , particle->getMomentum()[ 1 ] , particle->getMomentum()[ 2 ] , particle->getEnergy() );
			InvariantMass = ( chargedFourMomentum + neutralFourMomentum ).M();
			++addedParticles;
		}
		++i_par;
	}
	return addedParticles;
}

std::vector<EVENT::ReconstructedParticle*> sortParticles( std::vector<EVENT::ReconstructedParticle*> sortedParticles , std::vector<EVENT::ReconstructedParticle*> remainedParticles , TVector3 pointingVector )
{
	std::vector<EVENT::ReconstructedParticle*> newSortedParticles;
	std::vector<EVENT::ReconstructedParticle*> newRemainedParticles;
	for ( unsigned int i_particle = 0 ; i_particle < remainedParticles.size() ; ++i_particle )
	{
		newSortedParticles.push_back( sortedParticles[ i_particle ] );
	}
	double closestAngleCos = -1.0;
	unsigned int closestParticleIndex = 0;
	for ( unsigned int i_particle = 0 ; i_particle < remainedParticles.size() ; ++i_particle )
	{
		EVENT::ReconstructedParticle* particle = remainedParticles[ i_particle ];
		TVector3 Momentum( particle->getMomentum()[ 0 ] , particle->getMomentum()[ 1 ] , particle->getMomentum()[ 2 ] );
		Momentum.SetMag( 1.0 );
		if ( Momentum.Dot( pointingVector ) > closestAngleCos )
		{
			closestAngleCos = Momentum.Dot( pointingVector );
			closestParticleIndex = i_particle;
		}
	}
	for ( unsigned int i_particle = 0 ; i_particle < remainedParticles.size() ; ++i_particle )
	{
		if ( i_particle == closestParticleIndex )
		{
			newSortedParticles.push_back( remainedParticles[ i_particle ] );
		}
		else
		{
			newRemainedParticles.push_back( remainedParticles[ i_particle ] );
		}
	}
	newSortedParticles = sortParticles( newSortedParticles , newRemainedParticles , pointingVector );
	return newSortedParticles;
}
