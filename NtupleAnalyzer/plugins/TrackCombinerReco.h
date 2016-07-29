//#ifndef Monopoles_TrackCombinerReco_H
//#define Monopoles_TrackCombinerReco_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
//#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
//#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"

#include "DataFormats/Common/interface/Handle.h"

//#include "RecoTracker/DeDx/interface/UnbinnedLikelihoodFit.h"

//#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"

//#include <TLinearFitter.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TFile.h>
#include <TTree.h>

#include <vector>
#include <string>
#include <set>

//typedef std::vector<Trajectory> TrajectoryCollection;
//#endif

