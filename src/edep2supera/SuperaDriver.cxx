
#include "SuperaDriver.h"
#include "Utilities.h"
#include <iostream>
#include <map>

namespace edep2supera {
	
	void SuperaDriver::Configure(const std::string& name, const std::map<std::string,std::string>& params)
    {
        if(name == "edep2supera") {
            supera::PSet cfg;
            cfg.data = params;
			Configure(cfg);
			SetLogConfig(supera::Logger::parseStringThresh(cfg.get<std::string>("LogLevel","WARNING")));
        }
        else{
            std::string msg = name + " is not known to Supera...";
            throw supera::meatloaf(msg);
        }

    }

	void SuperaDriver::Configure(const supera::PSet &cfg)
	{
		if (cfg.exists("ActiveDetectors") > 0) 
			allowed_detectors = cfg.get<std::vector<std::string>>("ActiveDetectors");
		_segment_size_max = cfg.get<double>("MaxSegmentSize",0.03);
	}

	supera::EventInput SuperaDriver::ReadEvent(const TG4Event *ev) // returns a supera.Event to be used in SuperaAtomic
	{
		supera::EventInput result;
		LOG.VERBOSE() << "Processing " << ev->Trajectories.size() << " trajectories.\n";
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
			part_input.type    = this->InferProcessType(traj);

			LOG.VERBOSE() << "  Track ID " << part_input.part.trackid 
			<< " PDG " << part_input.part.pdg 
			<< " Energy " << part_input.part.energy_init << "\n";
			// Critical check: all track ID must be an integer >=0
			if(traj.GetTrackId() < 0) {
				LOG.FATAL() << "Negative track ID found " << traj.GetTrackId() << "\n";
				throw supera::meatloaf();
			}
			_trackid2idx.resize(traj.GetTrackId()+1,-1);
			_trackid2idx[traj.GetTrackId()] = part_input.part.id;
			result.push_back(part_input);
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
		LOG.DEBUG() << "starting voxelization"<< "\n";

		for (auto const &det : ev->SegmentDetectors)
		{

			LOG.DEBUG() << "Accepting list of active regions from config\n";
			if (std::find(allowed_detectors.begin(), allowed_detectors.end(), det.first) == allowed_detectors.end())
			{
				LOG.INFO() << det.first<< "not in acceptable active regions\n";
				continue;
			}

			for (auto const &hit : det.second)
			{

				auto track_id = hit.Contrib.front();

				for(size_t i=0; hit.Contrib.size()>1 && i<hit.Contrib.size(); ++i) {
					auto const& tid = hit.Contrib[i];
					int pdg = -1;
					double energy = -1;
					if(tid < _trackid2idx.size() && _trackid2idx[tid]>=0) {
						auto const& part = result[_trackid2idx[tid]].part;
						pdg = part.pdg;
						energy = part.energy_init; 
					}
					LOG.WARNING() << "A segment with multiple tracks: ID " << tid 
					<< " PDG " << pdg << " Energy " << energy << "\n";
				}

				if(track_id >=_trackid2idx.size() || _trackid2idx[track_id] < 0) {
					LOG.ERROR() << "Segment for invalid particle (Track ID " << track_id << " unknown)\n";
					continue;
				}

				auto& part = result[_trackid2idx[track_id]];
				auto& pcloud = part.pcloud;
				auto edeps = MakeEDeps(hit);
				LOG.DEBUG() << "Segment for track " << part.part.trackid << " " << part.part.pdg 
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
		result.pz = edepsim_part.GetInitialMomentum().Pz();										///< (x,y,z) component of particle's initial momentum
		result.vtx = supera::Vertex(start.X() / 10., start.Y() / 10., start.Z() / 10., start.T()); ///< (x,y,z,t) of particle's vertex information
		result.end_pt = supera::Vertex(end.X() / 10., end.Y() / 10., end.Z() / 10., end.T());		///< (x,y,z,t) at which particle disappeared from G4WorldVolume
		//result.process = edepsim_part.Points[0].GetProcess(); 										///< string identifier of the particle's creation process from Geant4
		result.energy_init = edepsim_part.GetInitialMomentum().E();										///< initial energy of the particle
		result.parent_trackid = edepsim_part.GetParentId();												///< PDG code of the parent particle


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
	supera::ProcessType SuperaDriver::InferProcessType(const TG4Trajectory& edepsim_part)
	{
		auto pdg_code = edepsim_part.GetPDGCode();
		auto g4type_main = edepsim_part.Points.front().GetProcess();
		auto g4type_sub  = edepsim_part.Points.front().GetSubprocess();


		if(pdg_code == 22) return supera::kPhoton;
		else if(std::abs(pdg_code) == 11) 
		{
			/*
			std::cout << "PDG " << pdg_code << " G4ProcessType " << g4type_main 
			<< " SubProcessType " << g4type_sub
			<< std::endl;
			*/

			/*
			if (g4type_sub==2) return kIonization;
			if (g4type_sub==13) return kPhotoElectron;
			if (g4type_sub==13) return kCompton;
			if (g4type_sub==14) return kConversion;




			std::string prc = mcpart.Process();
			if( prc == "muIoni" || prc == "hIoni" || prc == "muPairProd" )
				grp.type = kDelta;
			else if( prc == "muMinusCaptureAtRest" || prc == "muPlusCaptureAtRest" || prc == "Decay" )
				grp.type = kDecay;
			else if( prc == "compt"  )
				grp.type = kCompton;
			else if( prc == "phot"   )
				grp.type = kPhotoElectron;
			else if( prc == "eIoni"  )
				grp.type = kIonization;
			else if( prc == "conv"   )
				grp.type = kConversion;
			else if( prc == "primary")
				grp.type = kPrimary;
			else
				grp.type = kOtherShower;
			*/

			return supera::kOtherShower;			
		}
		else if(pdg_code == 2112) return supera::kNeutron;
		else 
		{
			
			return supera::kTrack;
		}
	}


}
