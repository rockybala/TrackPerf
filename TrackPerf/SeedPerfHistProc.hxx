#pragma once

#include <TH1.h>

#include <marlin/Processor.h>

namespace TrackPerf {


} // namespace TrackPerf

class SeedPerfHistProc : public marlin::Processor {
public:
  virtual Processor* newProcessor() { return new SeedPerfHistProc; }

  SeedPerfHistProc(const SeedPerfHistProc&) = delete;
  SeedPerfHistProc& operator=(const SeedPerfHistProc&) = delete;
  SeedPerfHistProc();

  virtual void init();
  virtual void processRunHeader(LCRunHeader* run);
  virtual void processEvent(LCEvent* evt);
  virtual void check(LCEvent* evt);
  virtual void end();

private:
  //! Seed collection
  std::string _seedColName{};

  //! MC Particle Collection
  std::string _mcpColName{};

  //! Seed to MC truth match collection
  std::string _seedMatchColName{};

  //! Histograms
  TH1* h_number_of_seeds;
  TH1* h_number_of_matched_seeds;
  TH1* h_number_of_fake_seeds;
  TH1* h_matched_mcparticle_pt;
  TH1* h_matched_mcparticle_lambda;
  TH1* h_matched_mcparticle_phi;
  TH1* h_matched_mcparticle_eta;
  TH1* h_matched_mcparticle_theta;
  
  
};
