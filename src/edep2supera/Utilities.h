#ifndef _EDEP2SUPERA_UTILITIES_H_
#define _EDEP2SUPERA_UTILITIES_H_

#include "supera/base/Point.h"
#include <vector>
namespace supera {
	std::vector<supera::Point3D> 
	SamplePointsFromLine(const supera::Point3D& p0,
		const supera::Point3D& p1,
		const double max_segment_size
		);

}

#endif