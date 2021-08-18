#include "streamlog/streamlog.h"
#include "marlin/Global.h"
#include "EVENT/Vertex.h"
#include "linkedPFO.h"
#include "flightDirection.h"

using namespace lcio;
using namespace marlin;

int getParentHadronFlightDirection( EVENT::LCEvent *pLCEvent , EVENT::MCParticle *SLDLepton , std::string inputPrimaryVertex , std::string inputBuildUpVertex , std::string inputJetCollection , int vertexinScenario , std::string recoMCTruthLinkCollection , std::string mcTruthRecoLinkCollection )
{
	std::vector<double> primaryVertex;
	std::vector<double> secondayVertex;
	std::vector<double> tritaryVertex;
	streamlog_out(DEBUG0) << "	Look for PFO linked to MCParticle (PDG: " << std::endl;

	return 1;
}

int getPrimaryVertex( EVENT::LCEvent *pLCEvent , EVENT::MCParticle *SLDLepton , std::string inputPrimaryVertex , bool cheatVertices , std::vector<double> &primaryVertex )
{
	int primaryVertexStatus = 0;

	try
	{
		LCCollection *primaryVertexCollection = pLCEvent->getCollection( inputPrimaryVertex );
		streamlog_out(DEBUG2) << "	There is " << primaryVertexCollection->getNumberOfElements() << " Primary Vertex" << std::endl;
	}
	catch (DataNotAvailableException &e)
	{
		streamlog_out(WARNING) << "Could not find the Primary Vertex collection" << std::endl;
		return -1;
	}

	if ( cheatVertices )
	{
		const EVENT::MCParticle *MotherHadron = SLDLepton->getParents()[ 0 ];
		primaryVertex.push_back( MotherHadron->getVertex()[ 0 ] );
		primaryVertex.push_back( MotherHadron->getVertex()[ 1 ] );
		primaryVertex.push_back( MotherHadron->getVertex()[ 2 ] );
		primaryVertexStatus = 1;
		streamlog_out(DEBUG4) << "		true primary Vertex (x,y,z): 		" << primaryVertex[ 0 ] << "	, " << primaryVertex[ 1 ] << "	, " << primaryVertex[ 2 ] << std::endl;
	}
	else
	{
		LCCollection *primaryVertexCollection = pLCEvent->getCollection( inputPrimaryVertex );
		Vertex* primaryVtx = dynamic_cast<Vertex*>( primaryVertexCollection->getElementAt( 0 ) );
		primaryVertex.push_back( primaryVtx->getPosition()[ 0 ] );
		primaryVertex.push_back( primaryVtx->getPosition()[ 1 ] );
		primaryVertex.push_back( primaryVtx->getPosition()[ 2 ] );
		primaryVertexStatus = 2;
		streamlog_out(DEBUG4) << "		reco primary Vertex (x,y,z): 	" << primaryVertex[ 0 ] << "	, " << primaryVertex[ 1 ] << "	, " << primaryVertex[ 2 ] << std::endl;
	}
	return primaryVertexStatus;
}

int getSecondaryVertex( EVENT::LCEvent *pLCEvent , EVENT::MCParticle *SLDLepton , std::string inputPrimaryVertex , std::string inputBuildUpVertex , bool cheatVertices , std::vector<double> &secondayVertex , std::string recoMCTruthLinkCollection , std::string mcTruthRecoLinkCollection )
{
	int secondayVertexStatus = 0;
	try
	{
		LCCollection *buildUpVertexCollection = pLCEvent->getCollection( inputBuildUpVertex );
		streamlog_out(DEBUG2) << "	There are " << buildUpVertexCollection->getNumberOfElements() << " BuildUp Vertices" << std::endl;
	}
	catch (DataNotAvailableException &e)
	{
		streamlog_out(WARNING) << "Could not find the BuildUp Vertex collection" << std::endl;
		secondayVertexStatus = -1;
	}
	if ( cheatVertices )
	{
		const EVENT::MCParticle *MotherHadron = SLDLepton->getParents()[ 0 ];
		secondayVertex.push_back( MotherHadron->getEndpoint()[ 0 ] );
		secondayVertex.push_back( MotherHadron->getEndpoint()[ 1 ] );
		secondayVertex.push_back( MotherHadron->getEndpoint()[ 2 ] );
		secondayVertexStatus = 1;
		streamlog_out(DEBUG4) << "		true secondary Vertex (x,y,z): 		" << secondayVertex[ 0 ] << "	, " << secondayVertex[ 1 ] << "	, " << secondayVertex[ 2 ] << std::endl;
	}
	else
	{
		ReconstructedParticle* linkedRecoLepton = getLinkedPFO( pLCEvent , SLDLepton , recoMCTruthLinkCollection , mcTruthRecoLinkCollection , true , false );
		LCCollection *buildUpVertexCollection = pLCEvent->getCollection( inputBuildUpVertex );
		int n_VTX = buildUpVertexCollection->getNumberOfElements();
		if ( linkedRecoLepton == NULL )
		{
			secondayVertexStatus = 1;
		}
		else
		{
			for ( int i_vtx = 0 ; i_vtx < n_VTX ; ++i_vtx )
			{
				Vertex* secondaryVtx = dynamic_cast<Vertex*>( buildUpVertexCollection->getElementAt( i_vtx ) );
				ReconstructedParticle* vertexRecoParticle = secondaryVtx->getAssociatedParticle();
				int nVertexParticles = ( vertexRecoParticle->getParticles() ).size();
				for ( int i_particle = 0 ; i_particle < nVertexParticles ; ++i_particle )
				{
					ReconstructedParticle* particle = vertexRecoParticle->getParticles()[ i_particle ];
					if ( linkedRecoLepton == particle )
					{
						secondayVertex.push_back( secondaryVtx->getPosition()[ 0 ] );
						secondayVertex.push_back( secondaryVtx->getPosition()[ 1 ] );
						secondayVertex.push_back( secondaryVtx->getPosition()[ 2 ] );
						secondayVertexStatus = 2;
						streamlog_out(DEBUG0) << "	Found Reco Lepton in BuildUp Vertex" << std::endl;
					}
				}
			}
		}
	}
	return secondayVertexStatus;
}
