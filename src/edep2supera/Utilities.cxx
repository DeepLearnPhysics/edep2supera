#include "Utilities.h"
#include <assert.h>
namespace supera {
	std::vector<supera::Point3D> 
	SamplePointsFromLine(const supera::Point3D& p0,
		const supera::Point3D& p1,
		double max_segment_size) 
	{

		assert(max_segment_size>0);

		double dist = p0.distance(p1);
		size_t num_segment = int(dist / max_segment_size) + 1;
		double len_segment = dist / (double(num_segment));

		auto dir = p0.direction(p1);

		std::vector<supera::Point3D> result;
		result.reserve(num_segment);

		// first point
		result.push_back(p0 + (dir * (len_segment/2)));
		for(size_t i=1; i<num_segment; ++i)
			result.push_back(result.back() + (dir * len_segment));

		return result;
	}
}

