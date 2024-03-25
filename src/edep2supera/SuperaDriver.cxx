
#include "SuperaDriver.h"
#include "Utilities.h"
#include <iostream>
#include <map>

namespace edep2supera {

	void SuperaDriver::Configure(const YAML::Node& cfg)
	{
		_allowed_detectors.clear();
		if (cfg["ActiveDetectors"]) {
			if(cfg["ActiveDetectors"].IsSequence())
				_allowed_detectors = cfg["ActiveDetectors"].as<std::vector<std::string> >();
			else
				_allowed_detectors.push_back(cfg["ActiveDetectors"].as<std::string>());
		}

		_segment_size_max = 0.03;
		if (cfg["MaxSegmentSize"])
			_segment_size_max = cfg["MaxSegmentSize"].as<double>();

		supera::Driver::Configure(cfg);
	}

	/*
	void SuperaDriver::ConfigureFromTest(const std::string& yaml_text)
	{
        auto cfg = YAML::LoadFile(yaml_text);

		_allowed_detectors.clear();
		if (cfg["ActiveDetectors"]) {
			if(cfg["ActiveDetectors"].IsSequence())
				_allowed_detectors = cfg["ActiveDetectors"].as<std::vector<std::string> >();
			else
				_allowed_detectors.push_back(cfg["ActiveDetectors"].as<std::string>());
		}

		_segment_size_max = 0.03;
		if (cfg["MaxSegmentSize"])
			_segment_size_max = cfg["MaxSegmentSize"].as<double>();

		supera::Driver::Configure(cfg);
	}

	void SuperaDriver::ConfigureFromFile(const std::string& yaml_file)
	{
        auto cfg = YAML::LoadFile(yaml_file);

		_allowed_detectors.clear();
		if (cfg["ActiveDetectors"]) {
			if(cfg["ActiveDetectors"].IsSequence())
				_allowed_detectors = cfg["ActiveDetectors"].as<std::vector<std::string> >();
			else
				_allowed_detectors.push_back(cfg["ActiveDetectors"].as<std::string>());
		}

		_segment_size_max = 0.03;
		if (cfg["MaxSegmentSize"])
			_segment_size_max = cfg["MaxSegmentSize"].as<double>();

		supera::Driver::Configure(cfg);
	}
	*/
	supera::EventInput SuperaDriver::ReadEvent(const TG4Event *ev) // returns a supera.Event to be used in SuperaAtomic
	{
		supera::EventInput result;
		LOG_VERBOSE() << "Processing " << ev->Trajectories.size() << " trajectories.\n";
		result.reserve(ev->Trajectories.size());

		// index map from track_id to result vector index 
		// (do not use std::map but a simple vector index mapping)
		_trackid2idx.clear();
		_trackid2idx.reserve(ev->Trajectories.size());

		for (auto const &traj : ev->Trajectories)
		{
			supera::ParticleInput part_input;

			part_input.valid   = true;
			part_input.part    = this->TG4TrajectoryToParticle(traj);
			part_input.part.id = result.size();
			//part_input.part.type = this->InferProcessType(traj,part_input.part);

			LOG_VERBOSE() << "  Track ID " << part_input.part.trackid 
			<< " PDG " << part_input.part.pdg 
			<< " Energy " << part_input.part.energy_init << "\n";
			// Critical check: all track ID must be an integer >=0
			if(traj.GetTrackId() < 0) {
				LOG_FATAL() << "Negative track ID found " << traj.GetTrackId() << "\n";
				throw supera::meatloaf();
			}
			_trackid2idx.resize(traj.GetTrackId()+1,supera::kINVALID_INDEX);
			_trackid2idx[traj.GetTrackId()] = part_input.part.id;
			result.push_back(part_input);
		}

		// Fill parent information (needed before TG4TrajectoryToParticle is called)
		for (size_t i=0; i<result.size(); ++i) {
			auto& part = result[i].part;
			auto const& traj = ev->Trajectories[i];

			if(part.parent_trackid < _trackid2idx.size()) {
				auto const& parent_index = _trackid2idx[part.parent_trackid];
				if(parent_index != supera::kINVALID_INDEX)
					part.parent_pdg = result[parent_index].part.pdg;
			}
			this->SetProcessType(traj,part);
		}

		VoxelizeEvent(ev,result);

		return result;
	}

	std::vector<supera::EDep> 
	SuperaDriver::MakeEDeps(const TG4HitSegment &hit) const {

		double energy = hit.GetEnergyDeposit();
		auto const& start = hit.GetStart();
		auto const& end   = hit.GetStop();

		supera::Point3D pt_start(start.X()/10.,start.Y()/10.,start.Z()/10.);
		supera::Point3D pt_end  (end.X()/10.,end.Y()/10.,end.Z()/10.);
		auto points = SamplePointsFromLine(pt_start,pt_end,_segment_size_max);

		double segment_energy = energy / ((double)(points.size()));
		double segment_size   = pt_start.distance(pt_end) / ((double)(points.size()));
		double segment_dedx   = segment_energy / segment_size;

		std::vector<supera::EDep> result(points.size());
		for(size_t i=0; i<points.size(); ++i) {
			result[i].x = points[i].x;
			result[i].y = points[i].y;
			result[i].z = points[i].z;
			result[i].e = segment_energy;
			result[i].dedx = segment_dedx;
		}

		return result;
	}

	void SuperaDriver::VoxelizeEvent(const TG4Event *ev, 
		supera::EventInput &result) const
	{
		LOG_DEBUG() << "starting voxelization"<< "\n";

		for (auto const &det : ev->SegmentDetectors)
		{

			LOG_DEBUG() << "Accepting list of active regions from config\n";
			if (std::find(_allowed_detectors.begin(), _allowed_detectors.end(), det.first) == _allowed_detectors.end())
			{
				LOG_INFO() << det.first<< "not in acceptable active regions\n";
				continue;
			}

			for (auto const &hit : det.second)
			{

				auto track_id = hit.Contrib.front();

				for(size_t i=0; hit.Contrib.size()>1 && i<hit.Contrib.size(); ++i) {
					auto const& tid = hit.Contrib[i];
					int pdg = -1;
					double energy = -1;
					if(tid < (int)(_trackid2idx.size()) && _trackid2idx[tid] != supera::kINVALID_INDEX) {
						auto const& part = result[_trackid2idx[tid]].part;
						pdg = part.pdg;
						energy = part.energy_init; 
					}
					LOG_WARNING() << "A segment with multiple tracks: ID " << tid 
					<< " PDG " << pdg << " Energy " << energy << "\n";
				}

				if(track_id >= (int)(_trackid2idx.size()) || _trackid2idx[track_id] == supera::kINVALID_INDEX) {
					LOG_ERROR() << "Segment for invalid particle (Track ID " << track_id << " unknown)\n";
					continue;
				}

				auto& part = result[_trackid2idx[track_id]];
				auto& pcloud = part.pcloud;
				auto edeps = MakeEDeps(hit);
				LOG_DEBUG() << "Segment for track " << part.part.trackid << " " << part.part.pdg 
				<< " Energy total " << hit.GetEnergyDeposit() 
				<< " (" << hit.GetSecondaryDeposit() << ")" 
				<< " dE/dX " << edeps[0].dedx
				<< " length " << hit.GetEnergyDeposit() / edeps[0].dedx << "\n";
				if (pcloud.size() == 0) pcloud = edeps;
				else {
					pcloud.reserve(pcloud.size()+edeps.size());
					for(auto& edep : edeps) pcloud.push_back(edep);
				}


			}
		}
	}


	supera::Particle SuperaDriver::TG4TrajectoryToParticle(const TG4Trajectory& edepsim_part)
	{
		supera::Particle result;
		auto const& start = edepsim_part.Points.front().GetPosition();
		auto const& end   = edepsim_part.Points.back().GetPosition();

		// first, last steps, distance travel, energy_deposit, num_voxels

		result.trackid = edepsim_part.GetTrackId(); ///< Geant4 track id
		result.pdg = edepsim_part.GetPDGCode();	 ///< PDG code
		result.px = edepsim_part.GetInitialMomentum().Px();
		result.py = edepsim_part.GetInitialMomentum().Py();
		result.pz = edepsim_part.GetInitialMomentum().Pz();				///< (x,y,z) component of particle's initial momentum
		result.end_px = edepsim_part.Points.back().GetMomentum().X();
		result.end_py = edepsim_part.Points.back().GetMomentum().Y();
		result.end_pz = edepsim_part.Points.back().GetMomentum().Z();							  
		result.vtx = supera::Vertex(start.X() / 10., start.Y() / 10., start.Z() / 10., start.T()); ///< (x,y,z,t) of particle's vertex information
		result.end_pt = supera::Vertex(end.X() / 10., end.Y() / 10., end.Z() / 10., end.T());		///< (x,y,z,t) at which particle disappeared from G4WorldVolume
		//result.process = edepsim_part.Points[0].GetProcess(); 										///< string identifier of the particle's creation process from Geant4
		result.energy_init = edepsim_part.GetInitialMomentum().E();										///< initial energy of the particle
		if (edepsim_part.GetParentId() < -1) {
			LOG_FATAL() << "Parent ID " << edepsim_part.GetParentId() << " is unexpected (cannot be < -1)\n";
			throw supera::meatloaf();
		}

		if(edepsim_part.GetParentId() == -1)
			result.parent_trackid = result.trackid;
		else
			result.parent_trackid = edepsim_part.GetParentId();

		if(result.trackid == supera::kINVALID_TRACKID || result.parent_trackid == supera::kINVALID_TRACKID) {
			LOG_FATAL() << "Unexpected to have an invalid track ID " << edepsim_part.GetTrackId() 
			<< " or parent track ID " << edepsim_part.GetParentId() << "\n";
			throw supera::meatloaf();
		}
		for (size_t i = 0; i < edepsim_part.Points.size() - 1; ++i)
		{
			result.dist_travel += (edepsim_part.Points[i].GetPosition() - edepsim_part.Points[i + 1].GetPosition()).Vect().Mag();
		}
		return result;
	}
/*
	void SuperaDriver::BBox_bounds(const TG4HitSegment &deposition, supera::ParticleInput &PI)
	{
        if (deposition.GetEnergyDeposit()>0)
        {
        //   std::cout << p.trackid << " " << p.energy_deposit << " adding first " << p.first_step.pos.x << " " << p.first_step.pos.y << " " << p.first_step.pos.z << std::endl;
        //   std::cout << p.trackid << " " << p.energy_deposit << " adding last " << p.last_step.pos.x << " " << p.last_step.pos.y << " " << p.last_step.pos.z << std::endl;
		const double epsilon = 1.e-3;
		PI.edep_bottom_left.x = std::min(deposition.GetStart().X() / 10. - epsilon, PI.edep_bottom_left.x);
		PI.edep_bottom_left.y = std::min(deposition.GetStart().Y() / 10. - epsilon, PI.edep_bottom_left.y);
		PI.edep_bottom_left.z = std::min(deposition.GetStart().Z() / 10. - epsilon, PI.edep_bottom_left.z);
		PI.edep_bottom_left.x = std::min(deposition.GetStop().X() / 10. - epsilon, PI.edep_bottom_left.x);
		PI.edep_bottom_left.y = std::min(deposition.GetStop().Y() / 10. - epsilon, PI.edep_bottom_left.y);
		PI.edep_bottom_left.z = std::min(deposition.GetStop().Z() / 10. - epsilon, PI.edep_bottom_left.z);

		PI.edep_top_right.x = std::max(deposition.GetStart().X() / 10. + epsilon, PI.edep_top_right.x);
		PI.edep_top_right.y = std::max(deposition.GetStart().Y() / 10. + epsilon, PI.edep_top_right.y);
		PI.edep_top_right.z = std::max(deposition.GetStart().Z() / 10. + epsilon, PI.edep_top_right.z);
		PI.edep_top_right.x = std::max(deposition.GetStop().X() / 10. + epsilon, PI.edep_top_right.x);
		PI.edep_top_right.y = std::max(deposition.GetStop().Y() / 10. + epsilon, PI.edep_top_right.y);
		PI.edep_top_right.z = std::max(deposition.GetStop().Z() / 10. + epsilon, PI.edep_top_right.z);
		}
	}
*/

	/*
	Main process type enum from G4ProcessType.hh
	enum G4ProcessType
	{
		fNotDefined,
		fTransportation,
		fElectromagnetic,
		fOptical,             
		fHadronic,
		fPhotolepton_hadron,
		fDecay,
		fGeneral,
		fParameterisation,
		fUserDefined,
		fParallel
	};

	Sub process type enum from G4EmProcessSubType
	enum G4EmProcessSubType 
	{ 
		fCoulombScattering = 1, 
		fIonisation = 2, 
		fBremsstrahlung = 3, 
		fPairProdByCharged = 4,
		fAnnihilation = 5, 
		fAnnihilationToMuMu = 6,
		fAnnihilationToHadrons = 7,
		fNuclearStopping = 8,
		fMultipleScattering = 10, 
		fRayleigh = 11,
		fPhotoElectricEffect = 12,
		fComptonScattering = 13,
		fGammaConversion = 14,
		fGammaConversionToMuMu = 15,
		fCerenkov = 21,
		fScintillation = 22,
		fSynchrotronRadiation = 23,
		fTransitionRadiation = 24 
	};

	Sub process type enum from G4HadronicProcessType.hh
	enum G4HadronicProcessType
	{
		fHadronElastic =    111,
		fHadronInelastic =  121,
		fCapture =          131,
		fFission =          141,
		fHadronAtRest =     151,
		fLeptonAtRest =     152,
		fChargeExchange =   161,
		fRadioactiveDecay = 210
	};
	*/

	void
	SuperaDriver::SetProcessType(const TG4Trajectory& edepsim_part, 
		supera::Particle& supera_part)
	{

		auto pdg_code    = supera_part.pdg;
		auto g4type_main = edepsim_part.Points.front().GetProcess();
		auto g4type_sub  = edepsim_part.Points.front().GetSubprocess();

		std::stringstream ss;
		ss << (int)(g4type_main) << "::" << (int)(g4type_sub);

		supera_part.process = ss.str();
		if(pdg_code == 2112 || pdg_code > 1000000000)
			supera_part.type = supera::kNeutron;
		else if(supera_part.trackid == supera_part.parent_trackid) {
			supera_part.type = supera::kPrimary;
		}else if(pdg_code == 22) {
			supera_part.type = supera::kPhoton;
		}else if(std::abs(pdg_code) == 11) {
			if( g4type_main == TG4TrajectoryPoint::G4ProcessType::kProcessElectromagetic ) {
				if( g4type_sub == TG4TrajectoryPoint::G4ProcessSubtype::kSubtypeEMPhotoelectric ) {
					supera_part.type = supera::kPhotoElectron;
				}
				else if( g4type_sub == TG4TrajectoryPoint::G4ProcessSubtype::kSubtypeEMComptonScattering ) {
					supera_part.type = supera::kCompton;
				}
				else if (g4type_sub == TG4TrajectoryPoint::G4ProcessSubtype::kSubtypeEMGammaConversion ||
						 g4type_sub == TG4TrajectoryPoint::G4ProcessSubtype::kSubtypeEMPairProdByCharged)
				{
					supera_part.type = supera::kConversion;
				}else if( g4type_sub == TG4TrajectoryPoint::G4ProcessSubtype::kSubtypeEMIonization ) {
					if( std::abs(supera_part.parent_pdg) == 11 ) {
						supera_part.type = supera::kIonization;
					}else if(std::abs(supera_part.parent_pdg) == 211 || 
						std::abs(supera_part.parent_pdg) == 13 || 
						std::abs(supera_part.parent_pdg) == 2212 ||
						std::abs(supera_part.parent_pdg) == 321) {
						supera_part.type = supera::kDelta;
					}else if(supera_part.parent_pdg == 22){
						supera_part.type = supera::kCompton;
					}else{
						LOG_WARNING() << "UNEXPECTED CASE for IONIZATION " << std::endl
						<< "PDG " << pdg_code 
						<< " TrackId " << edepsim_part.TrackId
						<< " Energy " << supera_part.energy_init 
						<< " Parent PDG " << supera_part.parent_pdg 
						<< " Parent TrackId " << edepsim_part.ParentId
						<< " G4ProcessType " << g4type_main 
						<< " SubProcessType " << g4type_sub
						<< std::endl;
						//throw supera::meatloaf();
						supera_part.type = supera::kIonization;
					}
				}
				else{
					LOG_FATAL() << "UNEXPECTED EM SubType " << std::endl
					<< "PDG " << pdg_code 
					<< " TrackId " << edepsim_part.TrackId
					<< " Energy " << supera_part.energy_init 
					<< " Parent PDG " << supera_part.parent_pdg 
					<< " Parent TrackId " << edepsim_part.ParentId
					<< " G4ProcessType " << g4type_main 
					<< " SubProcessType " << g4type_sub
					<< std::endl;
					throw supera::meatloaf();
				}
			}
			else if( g4type_main == TG4TrajectoryPoint::G4ProcessType::kProcessDecay ) {
				supera_part.type = supera::kDecay;
			}else{
				LOG_WARNING() << "Cannot classify this shower" << std::endl 
				<< "PDG " << pdg_code 
				<< " TrackId " << edepsim_part.TrackId
				<< " Energy " << supera_part.energy_init 
				<< " Parent PDG " << supera_part.parent_pdg 
				<< " Parent TrackId " << edepsim_part.ParentId
				<< " G4ProcessType " << g4type_main 
				<< " SubProcessType " << g4type_sub
				<< std::endl;
				supera_part.type = supera::kOtherShower;
			}
		}
		else {
			supera_part.type = supera::kTrack;
		}

	} 


}
