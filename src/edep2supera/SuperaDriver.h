#ifndef __SUPERADRIVER_H__
#define __SUPERADRIVER_H__

#include "EDepSim/TG4Event.h"
#include "supera/data/Particle.h"
//#include "supera/data/Neutrino.h"
//#include "LArCVBaseUtilFunc.h"
#include "supera/process/Driver.h"
//#include "supera/algorithm/ParticleIndex.h"
//#include "Voxelize.h"
#include "supera/base/meatloaf.h"
#include <algorithm>

namespace edep2supera {

	class SuperaDriver : public supera::Driver
	{
	public:
		SuperaDriver(std::string name="SuperaDriver") : Driver(name), _segment_size_max(0.03) {}

		supera::EventInput ReadEvent(const TG4Event *ev);
		void VoxelizeEvent(const TG4Event *ev, supera::EventInput&result ) const;
		supera::Particle TG4TrajectoryToParticle(const TG4Trajectory& edepsim_part);
		void SetProcessType(const TG4Trajectory& edepsim_part,
			supera::Particle& supera_part);
		void ExpandBBox(supera::EventInput& result);

		virtual void Configure(const YAML::Node& cfg) override;

		std::vector<supera::EDep> MakeEDeps(const TG4HitSegment &hit) const;

	private:
		double _segment_size_max; // in CM
		std::vector<supera::Index_t> _trackid2idx;
		std::vector<std::string> _allowed_detectors;
	};
	
}
#endif
