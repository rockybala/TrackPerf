#include "TrackPerf/SeedPerfHistProc.hxx"

#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <EVENT/LCRelation.h>
#include <EVENT/MCParticle.h>
#include <EVENT/Track.h>

#include <AIDA/ITree.h>
#include <marlin/AIDAProcessor.h>
#include <set>

#include "TrackPerf/ResoHists.hxx"
#include "TrackPerf/TrackHists.hxx"
#include "TrackPerf/TruthHists.hxx"

SeedPerfHistProc aSeedPerfHistProc;

SeedPerfHistProc::SeedPerfHistProc() : Processor("SeedPerfHistProc") {
  // modify processor description
  _description =
      "SeedPerfHistProc creates a series of output histograms for track seeding "
      "performance studies.";

  registerInputCollection(LCIO::MCPARTICLE, "MCParticleCollection",
                          "Name of the MCParticle collection", _mcpColName,
                          _mcpColName);

  registerInputCollection(LCIO::TRACK, "SeedCollection",
                          "Name of the seed collection", _seedColName,
                          _seedColName);

  registerInputCollection(
      LCIO::LCRELATION, "MCTrackRelationCollection",
      "Name of LCRelation collection with track to MC matching",
      _seedMatchColName, _seedMatchColName);
}

void SeedPerfHistProc::init() {
  // Print the initial parameters
  printParameters();

  // Create histograms
  AIDA::ITree* tree = marlin::AIDAProcessor::tree(this);
  marlin::AIDAProcessor::histogramFactory(this);

  tree->mkdir("seeding");
  tree->cd("seeding");

  h_number_of_seeds = new TH1F("h_number_of_seeds", ";Number of seeds;Events", 100, 0, 300000);
  h_number_of_matched_seeds = new TH1F("h_number_of_matched_seeds", ";Number of matched seeds;Events", 100, 0, 300000);
  h_number_of_fake_seeds = new TH1F("h_number_of_fake_seeds", ";Number of fake seeds;Events", 100, 0, 300000);
  h_matched_mcparticle_pt = new TH1F("h_matched_mcparticle_pt", ";Particle pT;Particles", 100, 0, 5000);
  h_matched_mcparticle_lambda = new TH1F("h_matched_mcparticle_lambda", ";Particle #lambda; Particles", 100, -3.14, 3.14);
  h_matched_mcparticle_phi = new TH1F("h_matched_mcparticle_phi", ";Particle #phi; Particles", 100, -3.14, 3.14);
  h_matched_mcparticle_eta = new TH1F("h_matched_mcparticle_eta", ";Particle #eta;Events", 100, -5, 5);
  h_matched_mcparticle_theta = new TH1F("h_matched_mcparticle_theta", ";Particle #theta;Events", 36, 0, 180);
  
}

void SeedPerfHistProc::processRunHeader(LCRunHeader* /*run*/) {}

void SeedPerfHistProc::processEvent(LCEvent* evt) {

  // MCParticles
  LCCollection* mcpCol = evt->getCollection(_mcpColName);

  if (mcpCol->getTypeName() != lcio::LCIO::MCPARTICLE) {
    throw EVENT::Exception("Invalid collection type: " + mcpCol->getTypeName());
  }

  std::set<const EVENT::MCParticle*> mcpSet;
  for (uint32_t i = 0; i < mcpCol->getNumberOfElements(); i++) {
    const EVENT::MCParticle* mcp =
        static_cast<const EVENT::MCParticle*>(mcpCol->getElementAt(i));

    if (mcp->getGeneratorStatus() != 1) {
      continue;
    }

    if (mcp->getCharge() == 0) {
      continue;
    }

    if (mcp->isDecayedInTracker()) {
      continue;
    }

    // Tracker acceptance
    const double* mom = mcp->getMomentum();
    double pt = std::sqrt(std::pow(mom[0], 2) + std::pow(mom[1], 2));
    double lambda = std::atan2(mom[2], pt);
    if (fabs(lambda) > 75. / 180 * 3.14) {
      continue;
    }

    mcpSet.insert(mcp);
  }

  // Seeds

  LCCollection* seedCol = evt->getCollection(_seedColName);

  if (seedCol->getTypeName() != lcio::LCIO::TRACK) {
    throw EVENT::Exception("Invalid collection type: " + seedCol->getTypeName());
  }

  std::set<const EVENT::Track*> seedSet;
  for (uint32_t i = 0; i < seedCol->getNumberOfElements(); i++) {
    const EVENT::Track* seed =
        static_cast<const EVENT::Track*>(seedCol->getElementAt(i));

    seedSet.insert(seed);
  }

  h_number_of_seeds->Fill(seedSet.size());
  
  // Loop over track to MC associations to save matched objects
  std::set<const EVENT::MCParticle*> mcpSet_matched;
  std::set<const EVENT::Track*> seedSet_matched;
  LCCollection* seed2mcCol = evt->getCollection(_seedMatchColName);
  if (seed2mcCol->getTypeName() != lcio::LCIO::LCRELATION) {
    throw EVENT::Exception("Invalid collection type: " + seed2mcCol->getTypeName());
  }

  for (int i = 0; i < seed2mcCol->getNumberOfElements(); ++i) {
    const EVENT::LCRelation* rel =
        static_cast<const EVENT::LCRelation*>(seed2mcCol->getElementAt(i));
    const EVENT::MCParticle* mcp =
        static_cast<const EVENT::MCParticle*>(rel->getFrom());
    const EVENT::Track* seed = static_cast<const EVENT::Track*>(rel->getTo());

    if (mcpSet.count(mcp) == 0) {
      continue;
    }  // truth particle not selected

    if (seedSet.find(seed) != seedSet.end()) {
      mcpSet_matched.insert(mcp);
      seedSet_matched.insert(seed);
      
      mcpSet.erase(mcp);
      seedSet.erase(seed);
    }
    
  }

  h_number_of_matched_seeds->Fill(seedSet_matched.size());
  h_number_of_fake_seeds->Fill(seedSet.size());

  for (const auto* particle : mcpSet_matched) {
    const double* mom = particle->getMomentum();
    double pt = std::sqrt(std::pow(mom[0], 2) + std::pow(mom[1], 2));
    h_matched_mcparticle_pt->Fill(pt);

    double lambda = std::atan2(mom[2], pt);
    h_matched_mcparticle_lambda->Fill(lambda);

    double phi = std::atan2(mom[1], mom[0]);
    h_matched_mcparticle_phi->Fill(phi);

    double eta = std::atanh(mom[2] / std::sqrt(std::pow(pt, 2) + std::pow(mom[2], 2)));
    h_matched_mcparticle_eta->Fill(eta);

    double theta = 2 * std::atan(std::exp(-eta));
    h_matched_mcparticle_theta->Fill(theta);
  }
} 


void SeedPerfHistProc::check(LCEvent* /*evt*/) {}
void SeedPerfHistProc::end() {}

