#ifndef VR_STD_SORT_H
#define VR_STD_SORT_H

#include <CGAL/Kernel/global_functions_2.h>

#include "log.h"

template <class Kernel_,
          class Storage_>
class Vr_Std_Sort
{
public:
  /// \name template parameters
  //@{
  typedef Kernel_       Kernel;
  typedef Storage_      Storage;
  //@}
private:
  typedef typename Storage::Vertex_index  Vertex_index;
  typedef std::vector<Vertex_index>       Vertices;

  class CompBySlope
  {
  public:
    /// default constructor
    CompBySlope(const Storage& storage) :
      storage_(storage)
    {
    }

    /// \return true if a should be ranked before b
    bool operator() (const Vertex_index& a, const Vertex_index& b)
    {
      CGAL_precondition(a.assigned() && b.assigned());

      const CGAL::Point_2<Kernel>& point_a = storage_.vertex(a).point,
          point_b = storage_.vertex(b).point;

      if(CGAL::collinear(storage_.query_point(), point_a, point_b))
      {
        if(storage_.query_point().x() == point_a.x())
        {
          return (a < b);
        }

        if(storage_.allow_collinear())
        {
          return CGAL::less_x(point_a, point_b);
        }
        else
        {
          return CGAL::less_x(point_b, point_a);
        }
      }

      CGAL::Line_2<Kernel> line_a(storage_.query_point(), point_a),
          line_b(storage_.query_point(), point_b);
      return (CGAL::compare_slopes(line_a, line_b) == CGAL::SMALLER);
    }
  private:
    const Storage&  storage_;
  };

public:
  /// default constructor
  Vr_Std_Sort(Storage& storage) :
    comparator_(storage), storage_(storage)
  {

  }

  void insert_input_vertex(const Vertex_index& vertex)
  {
    CGAL_precondition(vertex.assigned());
    CGAL_precondition((0 <= vertex)
                      && (vertex < storage_.num_input_vertices()));

    vertices_.push_back(vertex);
  }

  void rank_vertices()
  {
    std::sort(vertices_.begin(), vertices_.end(), comparator_);

    for(typename Vertices::iterator it = vertices_.begin(); it != vertices_.end(); ++it)
    {
      // ignore split vertices
      if(!storage_.vertex(*it).is_split_vertex()
         && (!storage_.is_query_point_on_vertex()
             || (storage_.vertex(*it).point != storage_.query_point())))
      {
        storage_.rank_vertex(*it);
      }
    }
  }

private:
  CompBySlope comparator_;
  Storage&  storage_;
  Vertices  vertices_;
};

#endif // VR_STD_SORT_H
