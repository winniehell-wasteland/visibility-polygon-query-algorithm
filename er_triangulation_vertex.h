#ifndef ER_TRIANGULATION_VERTEX_H
#define ER_TRIANGULATION_VERTEX_H

#include <CGAL/Triangulation_vertex_base_2.h>

/**
 * custom vertex base for the edge rank triangulation
 */
template < typename Storage_, typename GT,
           typename Vb = CGAL::Triangulation_vertex_base_2<GT> >
class Er_Triangulation_vertex
    : public Vb
{
public:
  /// \name template parameters
  //@{
  typedef typename Vb::Face_handle  Face_handle;
  typedef typename Vb::Point        Point;
  typedef Storage_                  Storage;
  //@}

  /// template voodoo to reassign the triangulation data structure
  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS2>::Other    Vb2;
    typedef Er_Triangulation_vertex<Storage, GT, Vb2>   Other;
  };

  /// default constructor
  Er_Triangulation_vertex() :
    Vb(),
    vertex_index_()
  {

  }

  /// construct from point
  Er_Triangulation_vertex(const Point & p) :
    Vb(p),
    vertex_index_()
  {

  }

  /// construct from point and face
  Er_Triangulation_vertex(const Point & p, Face_handle c) :
    Vb(p, c),
    vertex_index_()
  {

  }

  /// construct from face
  Er_Triangulation_vertex(Face_handle c) :
    Vb(c),
    vertex_index_()
  {

  }

  const typename Storage::Vertex_index& vertex_index() const
  {
    return vertex_index_;
  }

  void set_vertex_index(const typename Storage::Vertex_index& index)
  {
    CGAL_precondition(index.assigned());
    CGAL_precondition(!vertex_index_.assigned());

    vertex_index_ = index;
  }
private:
  /// index of the vertex in storage
  typename Storage::Vertex_index  vertex_index_;
};

#endif // ER_TRIANGULATION_VERTEX_H
