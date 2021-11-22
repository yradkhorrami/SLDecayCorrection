#ifndef flightDirection_h_1
#define flightDirection_h_1


#include "lcio.h"
#include "marlin/Global.h"
#include "marlin/Processor.h"
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/Vertex.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/VertexImpl.h>
#include "marlinutil/HelixClass.h"
#include "UTIL/LCRelationNavigator.h"
#include "streamlog/streamlog.h"
#include "DDMarlinCED.h"
#include "GeometryUtil.h"

#include "TVector3.h"
#include "TF2.h"

#include "SLDCorrection.h"
#include "AssignParticlestoSLD.h"
//#include "linkedPFO.h"

typedef std::vector<EVENT::MCParticle*>			mcpVector;
typedef std::vector<EVENT::ReconstructedParticle*>	pfoVector;
typedef std::vector<EVENT::Vertex*>			vtxVector;


void getTrueFlightDirection( EVENT::MCParticle *SLDLepton , TVector3 &trueFlightDirection , std::vector<double> &trueStartVertex , std::vector<double> &trueSLDVertex );
int getRecoFlightDirection( EVENT::ReconstructedParticle *linkedRecoLepton , TVector3 &recoFlightDirection , EVENT::Vertex *primaryVertex , EVENT::Vertex *startVertex , vtxVector &SLDVertices , pfoVector &SLDVerticesRP , EVENT::ReconstructedParticle *assignedJet , vtxVector verticesInJet , pfoVector PFOswithAloneTracks , float &helicesDistance , int vertexingScenario );

double intersectTrackLine( EVENT::Track *track , EVENT::Vertex* primaryVertex , TVector3 momentumOfLine , std::vector<double> pointOnLine , std::vector<double> &PCAatTrack , std::vector<double> &PCAatLine );

double intersectTrackTrack( EVENT::Track *track1 , EVENT::Track *track2 , std::vector<double> &PCAatTrack1 , std::vector<double> &PCAatTrack2 );

/*




void cheatTrueFlightDirection( EVENT::MCParticle *SLDLepton , TVector3 &trueFlightDirection , bool m_displayEvent );
int getStartVertexPosition( EVENT::Vertex *startVertex , std::vector<double> &startVertexPosition , bool m_displayEvent );
int getEndVertexPosition( EVENT::Vertex *primaryVtx , std::vector<EVENT::Vertex *> buildUpVertices , EVENT::ReconstructedParticle *linkedRecoLepton , std::vector<double> &endVertexPosition , double helicesDistance , bool m_displayEvent );
//int getPrimaryVertex( EVENT::LCEvent *pLCEvent , EVENT::MCParticle *SLDLepton , std::string inputPrimaryVertex , bool cheatVertices , std::vector<double> &primaryVertex , bool m_displayEvent );
int getParentHadronFlightDirection( EVENT::LCEvent *pLCEvent , EVENT::MCParticle *SLDLepton , TVector3 &trueFlightDirection , TVector3 &recoFlightDirection , std::string inputPrimaryVertex , std::string inputBuildUpVertex , std::string inputJetCollection , int vertexingScenario , std::string recoMCTruthLinkCollection , std::string mcTruthRecoLinkCollection , float &helicesDistance , std::vector<double> &SecondaryVertexPar , bool m_displayEvent );
int getSecondaryVertex( EVENT::LCEvent *pLCEvent , EVENT::MCParticle *SLDLepton , std::string inputPrimaryVertex , std::string inputBuildUpVertex , bool cheatVertices , std::vector<double> &secondayVertex , std::string inputJetCollection , std::string recoMCTruthLinkCollection , std::string mcTruthRecoLinkCollection , float &helicesDistance , std::vector<double> &SecondaryVertexPar , bool m_displayEvent );
int getJetAxis( EVENT::LCEvent *pLCEvent , EVENT::MCParticle *SLDLepton , TVector3 &jetAxis , std::string inputJetCollection , std::string recoMCTruthLinkCollection , std::string mcTruthRecoLinkCollection );
int getLeadingParticleFlightDirection( EVENT::LCEvent *pLCEvent , EVENT::MCParticle *SLDLepton , TVector3 &leadingParticleFlightDirection , std::string inputJetCollection , std::string recoMCTruthLinkCollection , std::string mcTruthRecoLinkCollection );
float intersectHelixHelix( EVENT::ReconstructedParticle *linkedRecoLepton , TVector3 Momentum , std::vector<double> point , std::vector<double> &PCAatLeptonTrack , std::vector<double> &PCAatDownStreamLine );
float intersectHelixLine( EVENT::LCEvent *pLCEvent , EVENT::Track *track , TVector3 Momentum , std::vector<double> point , std::vector<double> &PCAatLeptonTrack , std::vector<double> &PCAatDownStreamLine , std::string inputPrimaryVertex );
int getTrueDownStreamVertex( EVENT::MCParticle *SLDLepton , TVector3 &trueDSVertex , TVector3 &trueDSVertexMomentum );
int getTrueVertex( EVENT::MCParticle *MotherParticle , TVector3 &trueDSVertex ,TVector3 &trueDSVertexMomentum );
*/
void drawMCParticles( EVENT::MCParticle *MotherHadron );
void drawReconstructedParticle( EVENT::ReconstructedParticle *recoParticle , EVENT::Vertex *primaryVertex , int colorCharged , int colorNeutral );

#endif
