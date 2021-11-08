#include "AssignParticlestoSLD.h"
#include "marlin/Global.h"
#include <marlin/Processor.h>
using namespace lcio ;

std::vector<EVENT::ReconstructedParticle*> assignNeutralParticles( 	EVENT::ReconstructedParticle *assignedJet , TVector3 pointingVector , double InvMass , TLorentzVector chargedFourMomentum ,
									TLorentzVector &neutralFourMomentum , double cosOpeningAngle )
{
	pointingVector.SetMag( 1.0 );
	int nParticles = ( assignedJet->getParticles() ).size();
	std::vector<EVENT::ReconstructedParticle*> assignedNeutralParticles;
	std::vector<EVENT::ReconstructedParticle*> sortedParticles;
	std::vector<EVENT::ReconstructedParticle*> remainedParticles;
	for ( int i_particle = 0 ; i_particle < nParticles ; ++i_particle )
	{
		ReconstructedParticle* particle = assignedJet->getParticles()[ i_particle ];
		if ( particle->getTracks().size() == 0 )
		{
			TVector3 Momentum( particle->getMomentum()[ 0 ] , particle->getMomentum()[ 1 ] , particle->getMomentum()[ 2 ] );
			Momentum.SetMag( 1.0 );
			if ( Momentum.Dot( pointingVector ) >= cosOpeningAngle )
			{
				remainedParticles.push_back( particle );
				streamlog_out(DEBUG1) << "		adding one particle to the cone of particle candidates" << std::endl;
				streamlog_out(DEBUG1) << &particle << std::endl;
			}
		}
	}
	sortedParticles = sortParticlesInCone( sortedParticles , remainedParticles , pointingVector );
	streamlog_out(DEBUG1) << "		" << sortedParticles.size() << " neutral particles sorted wrt the angular distance to the pointing vector" << std::endl;

	TLorentzVector totalFourMomentum = chargedFourMomentum + neutralFourMomentum;
	for ( unsigned int i_par = 0 ; i_par < sortedParticles.size() ; ++i_par )
	{
		ReconstructedParticle* particle = sortedParticles[ i_par ];
		TLorentzVector particleFourMomentum = TLorentzVector( particle->getMomentum()[ 0 ] , particle->getMomentum()[ 1 ] , particle->getMomentum()[ 2 ] , particle->getEnergy() );
		if ( ( totalFourMomentum + particleFourMomentum ).M() <= InvMass )
		{
			totalFourMomentum += particleFourMomentum;
			neutralFourMomentum += particleFourMomentum;
			assignedNeutralParticles.push_back( particle );
			streamlog_out(DEBUG1) << "		added one neutral particle to the semi-leptonic decay" << std::endl;
		}
	}
	return assignedNeutralParticles;
}

std::vector<EVENT::ReconstructedParticle*> assignChargedParticles( 	EVENT::LCEvent *pLCEvent , EVENT::ReconstructedParticle *assignedJet , std::string inputPrimaryVertex ,
									std::string inputBuildUpVertex , EVENT::ReconstructedParticle *linkedRecoLepton , TVector3 pointingVector ,
									double InvMass , TLorentzVector chargedFourMomentum , TLorentzVector &neutralFourMomentum , double cosOpeningAngle )
{
	pointingVector.SetMag( 1.0 );
	int nParticles = ( assignedJet->getParticles() ).size();
	std::vector<EVENT::ReconstructedParticle*> assignedChargedParticles;
	std::vector<EVENT::ReconstructedParticle*> sortedParticles;
	std::vector<EVENT::ReconstructedParticle*> remainedParticles;
	LCCollection *primaryVertexCollection = pLCEvent->getCollection( inputPrimaryVertex );
	Vertex *primaryVtx = dynamic_cast<EVENT::Vertex*>( primaryVertexCollection->getElementAt( 0 ) );
	ReconstructedParticle* primaryParticle = primaryVtx->getAssociatedParticle();
	int n_PrimParticles = primaryParticle->getParticles().size();
	LCCollection *buildUpVertexCollection = pLCEvent->getCollection( inputBuildUpVertex );
	int n_VTX = buildUpVertexCollection->getNumberOfElements();
	for ( int i_particle = 0 ; i_particle < nParticles ; ++i_particle )
	{
		ReconstructedParticle* particle = assignedJet->getParticles()[ i_particle ];
		if ( particle->getTracks().size() != 0 && linkedRecoLepton != particle )
		{
			bool hasTrackInPrimaryVertex = false;
			bool hasTrackInBuildUpVertex = false;
			for ( unsigned int i_track = 0 ; i_track < particle->getTracks().size() ; ++i_track )
			{
				Track* testTrack = particle->getTracks()[ i_track ];
				for ( int i_par = 0 ; i_par < n_PrimParticles ; ++i_par )
				{
					for ( unsigned int i_trk = 0 ; i_trk < primaryParticle->getParticles()[ i_par ]->getTracks().size() ; ++i_trk )
					{
						Track* trk = primaryParticle->getParticles()[ i_par ]->getTracks()[ i_trk ];
						if ( trk == testTrack ) hasTrackInPrimaryVertex = true;
					}

				}
				for ( int i_vtx = 0 ; i_vtx < n_VTX ; ++i_vtx )
				{
					Vertex* buildUpVtx = dynamic_cast<Vertex*>( buildUpVertexCollection->getElementAt( i_vtx ) );
					ReconstructedParticle* buildUpParticle = buildUpVtx->getAssociatedParticle();
					int n_buildUpParticles = buildUpParticle->getParticles().size();
					for ( int i_par = 0 ; i_par < n_buildUpParticles ; ++i_par )
					{
						for ( unsigned int i_trk = 0 ; i_trk < buildUpParticle->getParticles()[ i_par ]->getTracks().size() ; ++i_trk )
						{
							Track* trk = buildUpParticle->getParticles()[ i_par ]->getTracks()[ i_trk ];
							if ( trk == testTrack ) hasTrackInBuildUpVertex = true;
						}

					}
				}
			}
			if ( !hasTrackInPrimaryVertex && !hasTrackInBuildUpVertex )
			{
				TVector3 Momentum( particle->getMomentum()[ 0 ] , particle->getMomentum()[ 1 ] , particle->getMomentum()[ 2 ] );
				Momentum.SetMag( 1.0 );
				if ( Momentum.Dot( pointingVector ) >= cosOpeningAngle )
				{
					remainedParticles.push_back( particle );
					streamlog_out(DEBUG1) << "		adding one particle to the cone of particle candidates" << std::endl;
					streamlog_out(DEBUG1) << &particle << std::endl;
				}
			}
		}
	}
	sortedParticles = sortParticlesInCone( sortedParticles , remainedParticles , pointingVector );
	streamlog_out(DEBUG1) << "		" << sortedParticles.size() << " charged particles sorted wrt the angular distance to the pointing vector" << std::endl;

	TLorentzVector totalFourMomentum = chargedFourMomentum + neutralFourMomentum;
	for ( unsigned int i_par = 0 ; i_par < sortedParticles.size() ; ++i_par )
	{
		ReconstructedParticle* particle = sortedParticles[ i_par ];
		TLorentzVector particleFourMomentum = TLorentzVector( particle->getMomentum()[ 0 ] , particle->getMomentum()[ 1 ] , particle->getMomentum()[ 2 ] , particle->getEnergy() );
		if ( ( totalFourMomentum + particleFourMomentum ).M() <= InvMass )
		{
			totalFourMomentum += particleFourMomentum;
			chargedFourMomentum += particleFourMomentum;
			assignedChargedParticles.push_back( particle );
			streamlog_out(DEBUG1) << "		added one neutral particle to the semi-leptonic decay" << std::endl;
		}
	}
	return assignedChargedParticles;
}

std::vector<EVENT::ReconstructedParticle*> sortParticlesInCone( std::vector<EVENT::ReconstructedParticle*> sortedParticles , std::vector<EVENT::ReconstructedParticle*> remainedParticles , TVector3 pointingVector )
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
	newSortedParticles = sortParticlesInCone( newSortedParticles , newRemainedParticles , pointingVector );
	return newSortedParticles;
}
