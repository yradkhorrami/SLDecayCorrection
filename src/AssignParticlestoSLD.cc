#include "AssignParticlestoSLD.h"

using namespace lcio;
using namespace marlin;

std::vector<EVENT::ReconstructedParticle*> assignNeutralParticles( 	EVENT::ReconstructedParticle *assignedJet , TVector3 pointingVector , double InvMass ,
									TLorentzVector chargedFourMomentum , TLorentzVector &neutralFourMomentum , double cosOpeningAngle )
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
				streamlog_out(DEBUG1) << *particle << std::endl;
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
					streamlog_out(DEBUG1) << *particle << std::endl;
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

bool isParticleInVertex( EVENT::ReconstructedParticle *particle , EVENT::Vertex *vertex )
{
	bool particleIsInVertex = false;
	streamlog_out(DEBUG0) << "" << std::endl;
	streamlog_out(DEBUG1) << "----------------------------------------------------------------------" << std::endl;
	streamlog_out(DEBUG1) << "------------ Looking for particle (" << particle << ") in vertex -------------" << std::endl;
	streamlog_out(DEBUG1) << "----------------------------------------------------------------------" << std::endl;
	streamlog_out(DEBUG0) << "	Looking for particle:" << std::endl;
	streamlog_out(DEBUG0) << *particle << std::endl;
	streamlog_out(DEBUG0) << "	in vertex:" << std::endl;
	streamlog_out(DEBUG0) << *vertex << std::endl;
	ReconstructedParticle* associatedParticle = vertex->getAssociatedParticle();
	int nPar = ( associatedParticle->getParticles() ).size();
	streamlog_out(DEBUG0) << "	Associated Particle of the Vertex: " << std::endl;
	streamlog_out(DEBUG0) << *associatedParticle << std::endl;
	for ( int i_par = 0 ; i_par < nPar ; ++i_par )
	{
		ReconstructedParticle* testParticle = associatedParticle->getParticles()[ i_par ];
		streamlog_out(DEBUG0) << "	Checking particle " << i_par << " (" << testParticle << ") in vertex" << std::endl;
		if ( particle == testParticle )
		{
			particleIsInVertex = true;
			streamlog_out(DEBUG2) << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
			streamlog_out(DEBUG2) << "<<<<<<<<<<<<<<<<<<<<  Found particle in vertex  >>>>>>>>>>>>>>>>>>>>>>" << std::endl;
			streamlog_out(DEBUG2) << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
			streamlog_out(DEBUG2) << *testParticle << std::endl;
			streamlog_out(DEBUG2) << "" << std::endl;
		}
	}
	return particleIsInVertex;
}

std::vector<EVENT::ReconstructedParticle*> getParticlesWithAloneTracks( EVENT::ReconstructedParticle *particle , EVENT::ReconstructedParticle *jet ,
	 								EVENT::Vertex *primaryVertex , std::vector<EVENT::Vertex*> vertexVector )
{
	std::vector<EVENT::ReconstructedParticle*> particlesWithAloneTracks{};
	streamlog_out(DEBUG0) << "----------------------------------------------------------------------" << std::endl;
	streamlog_out(DEBUG0) << "-------- Looking for alone tracks in jet " << jet << " in jets ----------" << std::endl;
	streamlog_out(DEBUG0) << "----------------------------------------------------------------------" << std::endl;
	streamlog_out(DEBUG0) << *jet << std::endl;
	int nParticles = ( jet->getParticles() ).size();
	bool particleIsInPrimaryVertex = false;
	bool particleIsInAVertex = false;
	for ( int i_particle = 0 ; i_particle < nParticles ; ++i_particle )
	{
		ReconstructedParticle* testParticle = jet->getParticles()[ i_particle ];
		if ( ( testParticle->getTracks() ).size() == 0 ) continue;
		if ( testParticle == particle ) continue;
		streamlog_out(DEBUG0) << "	Checking particle " << i_particle << " (" << testParticle << ") in jet" << std::endl;
		streamlog_out(DEBUG0) << "	Looking for particle in primary vertex " << std::endl;
		streamlog_out(DEBUG0) << *primaryVertex << std::endl;
		particleIsInPrimaryVertex = isParticleInVertex( testParticle , primaryVertex );
		for ( unsigned int i_vtx = 0 ; i_vtx < vertexVector.size() ; ++i_vtx )
		{
			if ( !particleIsInAVertex )
			{
				Vertex *vertex = vertexVector[ i_vtx ];
				streamlog_out(DEBUG0) << "	Looking for particle in build up vertex[ " << i_vtx << " ]: " << vertex << std::endl;
				streamlog_out(DEBUG0) << *vertex << std::endl;
				particleIsInAVertex = isParticleInVertex( testParticle , vertex );
			}
		}
		if ( !particleIsInPrimaryVertex && !particleIsInAVertex ) particlesWithAloneTracks.push_back( testParticle );
	}
	return particlesWithAloneTracks;
}

EVENT::Vertex *getParticleVertex( EVENT::ReconstructedParticle *particle , std::vector<EVENT::Vertex*> vertexVector , bool &foundParticleInVertex )
{
	foundParticleInVertex = false;
	EVENT::Vertex *particleVertex = NULL;
	streamlog_out(DEBUG0) << "" << std::endl;
	streamlog_out(DEBUG0) << "----------------------------------------------------------------------" << std::endl;
	streamlog_out(DEBUG0) << "------ Looking for particle (" << particle << ") in vertices -----------------" << std::endl;
	streamlog_out(DEBUG0) << "----------------------------------------------------------------------" << std::endl;
	streamlog_out(DEBUG0) << *particle << std::endl;
	for ( unsigned int i_vtx = 0 ; i_vtx < vertexVector.size() ; ++i_vtx )
	{
		Vertex *vertex = vertexVector[ i_vtx ];
		streamlog_out(DEBUG0) << "	Looking for particle in vertex " << i_vtx << std::endl;
		streamlog_out(DEBUG0) << *vertex << std::endl;
		foundParticleInVertex = isParticleInVertex( particle , vertex );
		if ( foundParticleInVertex ) particleVertex = vertex;
	}
	return particleVertex;
}

EVENT::ReconstructedParticle *getJetAssignedToParticle( EVENT::ReconstructedParticle *particle , std::vector<EVENT::ReconstructedParticle*> jetVector , bool &foundParticleInJet )
{
	foundParticleInJet = false;
	EVENT::ReconstructedParticle *assignedJet = NULL;
	streamlog_out(DEBUG1) << "" << std::endl;
	streamlog_out(DEBUG1) << "----------------------------------------------------------------------" << std::endl;
	streamlog_out(DEBUG1) << "------------- Looking for particle (" << particle << ") in jets --------------" << std::endl;
	streamlog_out(DEBUG1) << "----------------------------------------------------------------------" << std::endl;
	streamlog_out(DEBUG1) << *particle << std::endl;
	for ( unsigned int i_jet = 0 ; i_jet < jetVector.size() ; ++i_jet )
	{
		EVENT::ReconstructedParticle* jet = jetVector[ i_jet ];
		streamlog_out(DEBUG0) << "	Looking for particle in jet " << i_jet << std::endl;
		streamlog_out(DEBUG0) << *jet << std::endl;
		int nParticles = ( jet->getParticles() ).size();
		for ( int i_particle = 0 ; i_particle < nParticles ; ++i_particle )
		{
			ReconstructedParticle* testParticle = jet->getParticles()[ i_particle ];
			streamlog_out(DEBUG0) << "	Checking particle " << i_particle << " (" << testParticle << ") in jet" << std::endl;
			if ( testParticle == particle )
			{
				streamlog_out(DEBUG2) << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
				streamlog_out(DEBUG2) << "<<<<<<<<<<<<<<<<<<<<<<  Found particle in jet  >>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
				streamlog_out(DEBUG2) << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
				streamlog_out(DEBUG2) << *testParticle << std::endl;
				streamlog_out(DEBUG2) << "" << std::endl;
				assignedJet = jet;
				foundParticleInJet = true;
			}
		}
	}
	return assignedJet;
}

std::vector<EVENT::Vertex*> getVerticesInJet( EVENT::ReconstructedParticle* jet , std::vector<EVENT::Vertex*> vertexVector )
{
	std::vector<EVENT::Vertex*> jetVertices{};
	for ( unsigned int i_vtx = 0 ; i_vtx < vertexVector.size() ; ++i_vtx )
	{
		Vertex *vertex = vertexVector[ i_vtx ];
		ReconstructedParticle* associatedParticle = vertex->getAssociatedParticle();
		for ( unsigned int i_par = 0 ; i_par < associatedParticle->getParticles().size() ; ++i_par )
		{
			ReconstructedParticle* vertexParticle = associatedParticle->getParticles()[ i_par ];
			int nParticles = ( jet->getParticles() ).size();
			for ( int i_particle = 0 ; i_particle < nParticles ; ++i_particle )
			{
				ReconstructedParticle* jetParticle = jet->getParticles()[ i_particle ];
				streamlog_out(DEBUG0) << "	Checking particle " << i_particle << " (" << jetParticle << ") in jet" << std::endl;
				if ( vertexParticle == jetParticle )
				{
					streamlog_out(DEBUG2) << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
					streamlog_out(DEBUG2) << "<<<<<<<<<<<<<<<<<<<<<<<  Found Vertex in Jet  >>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
					streamlog_out(DEBUG2) << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
					streamlog_out(DEBUG2) << "---------------------------     Vertex     ---------------------------" << std::endl;
					streamlog_out(DEBUG2) << *vertex << std::endl;
					streamlog_out(DEBUG2) << "---------------     Associated Particle Of Vertex     ----------------" << std::endl;
					streamlog_out(DEBUG2) << *associatedParticle << std::endl;
					streamlog_out(DEBUG2) << "----------------------------     Jet     -----------------------------" << std::endl;
					streamlog_out(DEBUG2) << *jet << std::endl;
					streamlog_out(DEBUG2) << "" << std::endl;
					bool vertexWasAdded = false;
					for ( unsigned int i_v = 0 ; i_v < jetVertices.size() ; ++i_v )
					{
						if ( jetVertices[ i_v ] == vertex ) vertexWasAdded = true;
					}
					if ( !vertexWasAdded ) jetVertices.push_back( vertex );
				}
			}
		}
	}
	return jetVertices;
}

EVENT::Track *makeTrackParameters( TVector3 momentum , TVector3 point , double charge )
{
	TrackImpl *track = new TrackImpl;
	double m_Bfield = MarlinUtil::getBzAtOrigin();
	double Px = momentum.Px();
	double Py = momentum.Py();
	double Pz = momentum.Pz();
	double Xp = point.X();
	double Yp = point.Y();
	double Zp = point.Z();
	const float refPoint[ 3 ] = { 0.0 , 0.0 , 0.0 };

	double Pt = std::sqrt( pow( Px , 2 ) + pow( Py , 2 ) );
	double tanLambda = Pz / Pt;
	double omega = charge * m_Bfield * 3.0e-4 / Pt;
	double phip = atan2( Py , Px );
	double xCenter = Xp + std::sin( phip ) / omega;
	double yCenter = Yp - std::cos( phip ) / omega;
	double d0 = 1.0 / omega - charge * std::sqrt( pow( xCenter , 2 ) + pow( yCenter , 2 ) );
	double phi0 = atan2( xCenter , -1.0 * yCenter );
	while ( phi0 < 0.0 )
	{
		phi0 += 2.0 * M_PI;
	}
	while ( phi0 >= 2.0 * M_PI )
	{
		phi0 -= 2.0 * M_PI;
	}
	double deltaPhi = phip - phi0;
	double zOffset = -1.0 * deltaPhi * tanLambda / std::sqrt( omega );
	double fullTurnsLength = ( zOffset - Zp ) * std::sqrt( omega ) / ( tanLambda * 2.0 * M_PI );
	int nCurls = 0;
	int n1,n2;
        if ( fullTurnsLength >= 0.0 )
	{
		n1 = int( fullTurnsLength );
		n2 = n1 + 1;
        }
        else
	{
		n1 = int( fullTurnsLength ) - 1;
		n2 = n1 + 1;
        }

        if ( fabs( n1 - fullTurnsLength ) < fabs( n2 - fullTurnsLength ) )
	{
		nCurls = n1;
        }
        else
	{
		nCurls = n2;
        }
	double z0 = Zp + ( deltaPhi + nCurls * 2.0 * M_PI ) * tanLambda / std::fabs( omega );
	track->setD0( d0 );
	track->setZ0( z0 );
	track->setPhi( phi0 );
	track->setOmega( omega );
	track->setTanLambda( tanLambda );
	track->setReferencePoint( refPoint );
	return track;

}
