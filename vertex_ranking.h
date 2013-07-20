#ifndef VERTEX_RANKING_H
#define VERTEX_RANKING_H

#include "vr_dual_plane.h"
#include "vr_std_sort.h"

/// concept VertexRanking
template <class Kernel, class Storage>
class VertexRanking
{
public:
  /// default constructor
  VertexRanking(Storage& storage);

  void insert_input_vertex(const typename Storage::Vertex_index& vertex);
  void rank_vertices();
};

#endif // VERTEX_RANKING_H
