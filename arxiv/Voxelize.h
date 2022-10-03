/**
 * \file Voxelize.h
 *
 * \brief Utilities for segmenting truth energy deposits into 3D voxels
 *
 * @author J. Wolcott <jwolcott@fnal.gov>
 */

#ifndef EDEP2SUPERA_VOXELIZE_H
#define EDEP2SUPERA_VOXELIZE_H

#include <vector>

#include "supera/data/Particle.h"
#include "supera/base/Voxel.h"
#include "supera/data/ImageMeta3D.h"
#include "EDepSim/TG4Event.h"
#include "raybox.h"
#include "supera/algorithm/Loggable.h"

namespace supera
{
  // some forward declarations
  template <typename T>
  class AABBox;

  template <typename T>
  class Vec3;
  typedef Vec3<float> Vec3f;
  typedef Vec3<double> Vec3d;

  /// Split an edep-sim TG4HitSegment true energy deposit into voxels
  ///
  /// \param[in]  hitSegment   The TG4HitSegment to operate on
  /// \param[in]  meta         Metadata about the voxel array
  /// \param[in] particles     Collection of true particles whose recorded energy deposits will be updated.  Pass empty vector to disable
  /// \return                  A vector of voxels
   std::vector<EDep> MakeEDeps(const ::TG4HitSegment &hitSegment,
                                       const ImageMeta3D &meta,
                                       EventInput &particles,
                                       double     &dist_travel);

   // --------------------------------------------------

   /// Split an edep-sim TG4HitSegment true energy deposit into voxels
   ///
   /// \param[in]  hitSegment   The TG4HitSegment to operate on
   /// \param[in]  meta         Metadata about the voxel array
   /// \param[in] particles     Collection of true particles whose recorded energy deposits will be updated.  Pass empty vector to disable
   /// \return                  A vector of voxels
   std::vector<EDep> MakeEDeps(const ::TG4HitSegment &hitSegment,
                               const ImageMeta3D &meta,
                               double &dist_travel);


   /// Split an edep-sim TG4HitSegment true energy deposit into voxels
   ///
   /// \param[in]  hitSegment   The TG4HitSegment to operate on
   /// \param[in]  meta         Metadata about the voxel array
   /// \param[in] particles     Collection of true particles whose recorded energy deposits will be updated.  Pass empty vector to disable
   /// \return                  A vector of voxels
   std::vector<EDep> MakeEDeps(const ::TG4HitSegment &hitSegment,
                               const ImageMeta3D &meta);

  /// Where (if anywhere) does a line segment intersect a given bounding box?
  /// (If the entire line segment is contained, the entry and exit points
  ///  will be set to the start and stop points provided.)
  ///
  /// \param bbox        Bounding box in question
  /// \param startPoint  Start point of the line segment in question
  /// \param stopPoint   Stop point of the line segment in question
  /// \param entryPoint  Computed entry point of the line segment into the box, if any
  /// \param exitPoint   Computed exit point of the line segment from the box, if any
  /// \return            Number of intersections (will be 0, 1, or 2)
  template <typename T>
  char Intersections(const AABBox<T> &bbox,
                     const TVector3 &startPoint,
                     const TVector3 &stopPoint,
                     Vec3<T> &entryPoint,
                     Vec3<T> &exitPoint);
}

#endif //EDEP2SUPERA_VOXELIZE_H
