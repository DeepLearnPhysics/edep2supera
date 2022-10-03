#ifndef __SUPERADRIVER_H__
#define __SUPERADRIVER_H__

#include "EDepSim/TG4Event.h"
#include "supera/data/Particle.h"
#include "supera/base/PSet.h"
#include "supera/algorithm/Loggable.h"
//#include "LArCVBaseUtilFunc.h"
#include "supera/process/Driver.h"
#include "supera/algorithm/ParticleIndex.h"
//#include "Voxelize.h"
#include "supera/base/meatloaf.h"
#include <algorithm>

namespace edep2supera {

	class SuperaDriver : public supera::Driver, public supera::Loggable
	{
	public:
		SuperaDriver() : _segment_size_max(0.03) {}

		supera::EventInput ReadEvent(const TG4Event *ev);
		void VoxelizeEvent(const TG4Event *ev, supera::EventInput&result ) const;
		supera::Particle TG4TrajectoryToParticle(const TG4Trajectory& edepsim_part);
		supera::ProcessType InferProcessType(const TG4Trajectory& edepsim_part);
		void ExpandBBox(supera::EventInput& result);
		void Configure(const std::string& name, const std::map<std::string,std::string>& params);
		void Configure(const supera::PSet &cfg);
		//void BBox_bounds(const TG4HitSegment &deposition, supera::ParticleInput &PI);
		std::vector<supera::EDep> MakeEDeps(const TG4HitSegment &hit) const;
		std::vector<std::string> allowed_detectors;
	private:
		double _segment_size_max; // in CM
		std::vector<int> _trackid2idx;
	};
	
}
#endif
