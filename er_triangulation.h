#ifndef ER_TRIANGULATION_H
#define ER_TRIANGULATION_H

#include <algorithm>
#include <bitset>
#include <cassert>
#include <deque>
#include <vector>

#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_euclidean_traits_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>

#include "er_triangulation_vertex.h"
#include "visibility_index.h"
#include "visibility_storage.h"

#include "log.h"

#undef debug
#if ENABLE_VP_DEBUG
  #define debug(msg) \
  { \
    std::cout << (msg) << std::endl; \
  }
#else
  #define debug(msg)
#endif

template <class Kernel_,
          class Storage_>
class Er_Triangulation
{
public:
  /// \name template parameters
  //@{
  typedef Kernel_         Kernel;
  typedef Storage_        Storage;
  //@}

private:
  /**
   * contains additional information of triangulation faces
   */
  struct Face_info
  {
    /// default constructor
    Face_info() :
      covered(),
      ranked(false)
    {
    }

    /// copy constructor
    Face_info(const Face_info& other) :
      covered(other.covered),
      ranked(other.ranked)
    {
    }

    /// stores for each incident edge whether it is covered by another one of this face
    std::bitset<3> covered;
    /// flag for the graph ranking
    bool ranked;
  };

  /// \name triangulation defines
  //@{
  typedef CGAL::Triangulation_euclidean_traits_2<Kernel>                             Traits;

  typedef CGAL::Constrained_triangulation_face_base_2<Kernel>                        Face_base;
  typedef CGAL::Triangulation_face_base_with_info_2<Face_info,Kernel,Face_base>      Extended_face;
  typedef Er_Triangulation_vertex<Storage, Kernel>                                   Vertex_base;
  typedef CGAL::Triangulation_data_structure_2<Vertex_base, Extended_face>           Data_structure;

  typedef CGAL::No_intersection_tag                                                  Itag;

  typedef CGAL::Constrained_Delaunay_triangulation_2<Traits, Data_structure, Itag>   Triangulation;

  typedef typename Triangulation::Face_handle                                        Face_handle;
  typedef typename Triangulation::Line_face_circulator                               Line_face_circulator;
  typedef typename Triangulation::size_type                                          Size;
  typedef typename Triangulation::Vertex_handle                                      Vertex_handle;
  //@}

  typedef std::vector<Vertex_handle>                                                 Vertex_map;
  typedef std::pair<Vertex_handle, Vertex_handle>                                    Vertex_pair;
  typedef std::vector<Vertex_pair>                                                   Split_edges;
  typedef std::vector<Vertex_handle>                                                 Split_vertices;

public:
  /**
   * default constructor
   */
  Er_Triangulation(Storage& storage) :
    split_edges_(),
    split_vertices_(),
    storage_(storage),
    triangulation_(),
    vertex_map_()
  {

  }

  /**
   * remove all query vertex dependent data
   */
  void clear()
  {
    CGAL_precondition(triangulation_.is_valid(true));
    CGAL_precondition(query_vertex_ != Vertex_handle());

    // remove split vertices
    CGAL_precondition(split_vertices_.size() == split_edges_.size());
    for(Size i = 0; i < split_vertices_.size(); ++i)
    {
      const Vertex_handle& split_vertex = split_vertices_[i];
      const Vertex_pair& split_edge = split_edges_[i];

      CGAL_precondition(CGAL::Segment_2<Kernel>(split_edge.first->point(), split_edge.second->point()).has_on(split_vertex->point()));

      triangulation_.remove_incident_constraints(split_vertex);
      triangulation_.remove(split_vertex);

      // restore input edge
      triangulation_.insert_constraint(split_edge.first, split_edge.second);
    }

    // remove query vertex
    if(!storage_.is_query_point_on_edge()
       && !storage_.is_query_point_on_vertex())
    {
      CGAL_precondition(!triangulation_.are_there_incident_constraints(query_vertex_));

      triangulation_.remove(query_vertex_);
    }

    query_vertex_ = Vertex_handle();

    split_vertices_.clear();
    split_edges_.clear();

    // reset face info
    for(typename Triangulation::All_faces_iterator face = triangulation_.all_faces_begin(); face != triangulation_.all_faces_end(); ++face)
    {
      face->info().covered.reset();
      face->info().ranked = false;
    }

    CGAL_postcondition(triangulation_.is_valid(true));
  }

#if ENABLE_VP_DEBUG
  /**
   * draw the triangulation vertices and edges
   */
  void draw(CairoImage& img)
  {
    CairoColor red(CGAL::RED), black(CGAL::BLACK);

    for(typename Triangulation::Finite_edges_iterator edge = triangulation_.finite_edges_begin();
        edge != triangulation_.finite_edges_end();
        ++edge)
    {
      img.draw_segment(edge->first->vertex(triangulation_.cw(edge->second))->point(),
               edge->first->vertex(triangulation_.ccw(edge->second))->point(),
               edge->first->is_constrained(edge->second)?red:black);
    }

    storage_.draw_labels(img);

    img.set_color(CGAL::VIOLET);
    for(typename std::vector<Vertex_handle>::iterator split_vertex = split_vertices_.begin(); split_vertex != split_vertices_.end(); ++split_vertex)
    {
      img.draw_point((*split_vertex)->point());
    }
  }
#endif

#if ENABLE_VP_DEBUG
  /**
   * draw the covering graph
   */
  void draw_graph(CairoImage& img)
  {
    CairoColor red(CGAL::RED), black(CGAL::BLACK);

    for(typename Triangulation::Finite_edges_iterator edge = triangulation_.finite_edges_begin();
        edge != triangulation_.finite_edges_end();
        ++edge)
    {
      img.draw_segment(edge->first->vertex(triangulation_.cw(edge->second))->point(),
               edge->first->vertex(triangulation_.ccw(edge->second))->point(),
               edge->first->is_constrained(edge->second)?red:black);
    }

    img.set_color(CGAL::YELLOW);

    for(typename Triangulation::Finite_faces_iterator face = triangulation_.finite_faces_begin();
        face != triangulation_.finite_faces_end();
        ++face)
    {
      for(int i = 0; i < 3; ++i)
      {
        if(face->info().covered.test(i))
        {
          int visible = face->info().covered.test(triangulation_.cw(i))?triangulation_.ccw(i):triangulation_.cw(i);
          CGAL::Point_2<Kernel> visible_midpoint = CGAL::midpoint(face->vertex(triangulation_.cw(visible))->point(), face->vertex(triangulation_.ccw(visible))->point());

          CGAL::Point_2<Kernel> edge_midpoint = CGAL::midpoint(face->vertex(triangulation_.cw(i))->point(), face->vertex(triangulation_.ccw(i))->point());
          img.draw_segment(visible_midpoint, edge_midpoint);
        }
      }
    }
  }
#endif

  /// insert an input vertex in the triangulation
  void insert_input_vertex(const typename Storage::Vertex_index& index)
  {
    CGAL_precondition(index.assigned());

    Vertex_handle triangulation_vertex = triangulation_.insert(storage_.vertex(index).point);
    triangulation_vertex->set_vertex_index(index);

    CGAL_precondition(index == vertex_map_.size());
    vertex_map_.push_back(triangulation_vertex);

#if ENABLE_VP_DEBUG
    bbox_input_.extend((double[2]){CGAL::to_double(storage_.vertex(index).point.x()), CGAL::to_double(storage_.vertex(index).point.y())});
#endif
  }

  void insert_split_vertex(const typename Storage::Edge_index& index)
  {
    CGAL_precondition(index.assigned());

    const typename Storage::Edge_entry& edge_entry = storage_.split_edge(index);
    CGAL_precondition(edge_entry.source.assigned()
                      && edge_entry.target.assigned());
    CGAL_precondition((edge_entry.source < vertex_map_.size())
                      && (edge_entry.target < vertex_map_.size()));

    Vertex_handle source = vertex_map_[edge_entry.source],
      target = vertex_map_[edge_entry.target];

    Face_handle face;
    int edge;

    bool found_edge = triangulation_.is_edge(source, target, face, edge);

    CGAL_precondition(found_edge);

    // ensure the split edge is part of input arrangement
    CGAL_precondition_msg(face->is_constrained(edge), "Split edge is no input edge.");

    if(found_edge)
    {
      Vertex_handle split_vertex = triangulation_.insert(storage_.split_vertex(index).point, Triangulation::EDGE, face, edge);

      split_vertex->set_vertex_index(storage_.num_input_vertices() + index);

      split_vertices_.push_back(split_vertex);
      split_edges_.push_back(Vertex_pair(source, target));

      CGAL_postcondition(split_vertices_.size() == split_edges_.size());
    }

    CGAL_postcondition(triangulation_.is_valid(true));
  }

  /// insert arrangement edge as triangulation constraint
  void insert_input_edge(const CGAL::Point_2<Kernel>& from, const CGAL::Point_2<Kernel>& to)
  {
    triangulation_.insert_constraint(from, to);

    CGAL_postcondition(triangulation_.is_valid(true));
  }

  /// locate query point in input and insert it into triangulation
  void insert_query_point(const CGAL::Point_2<Kernel>& query_point)
  {
    CGAL_precondition_msg((query_vertex_ == Vertex_handle()), "Query vertex should not be set yet.");
    CGAL_precondition(storage_.is_query_point_assigned());

    int li;
    typename Triangulation::Locate_type lt;
    Face_handle face = triangulation_.locate(query_point, lt, li);

    if((lt == Triangulation::EDGE) && face->is_constrained(li))
    {
      CGAL_precondition(storage_.is_query_point_on_edge());

      // store current edge
      split_edges_.push_back(vertex_pair(face, li));
    }
    else if(lt == Triangulation::VERTEX)
    {
      CGAL_precondition(storage_.is_query_point_on_vertex());
    }
    else
    {
      CGAL_precondition(!storage_.is_query_point_on_vertex()
                        && !storage_.is_query_point_on_edge());
    }

    query_vertex_ = triangulation_.insert(query_point, lt, face, li);

    if(storage_.is_query_point_on_vertex())
    {
      CGAL_postcondition_msg(triangulation_.are_there_incident_constraints(query_vertex_), "Query vertex can not be isolated.");
    }
    else if(storage_.is_query_point_on_edge())
    {
      split_vertices_.push_back(query_vertex_);
      CGAL_postcondition(split_vertices_.size() == split_edges_.size());
    }

    CGAL_postcondition(triangulation_.is_valid(true));
  }

  /**
   * calculates the edge ranks based on triangulation by topological sorting
   * Assumption: triangulation has already been calculated
   */
  void rank_edges()
  {
    CGAL_precondition(triangulation_.is_valid(true));
    CGAL_precondition(query_vertex_ != Vertex_handle());

    create_directed_graph();

#if ENABLE_VP_DEBUG
    {
      CairoImage img("covering_graph_" + tostr(visibilty_polygons_computed, 3), bbox_input_.bbox());
      draw_graph(img);
      storage_.draw_labels(img);
    }
#endif

    rank_graph();

    CGAL_postcondition(triangulation_.is_valid(true));
  }

private:
  Vertex_handle query_vertex_;

  Split_edges split_edges_;
  Split_vertices split_vertices_;

  Storage& storage_;
  Triangulation triangulation_;

  Vertex_map vertex_map_;

#if ENABLE_VP_DEBUG
  CGAL::Box_intersection_d::Box_d<double, 2> bbox_input_;
#endif

  /// run create_subgraph() for each incident face of the query vertex
  void create_directed_graph()
  {
    CGAL_precondition(query_vertex_ != Vertex_handle());
    CGAL_precondition(triangulation_.is_valid(true));

    typename Triangulation::Face_circulator begin = triangulation_.incident_faces(query_vertex_), face = begin;

    do
    {
      if(triangulation_.is_infinite(face))
      {
        continue;
      }

      CGAL_precondition(face->info().covered.count() == 0);

      for(int i = 0; i < 3; i++)
      {
        if(face->vertex(i) == query_vertex_)
        {
          face->info().covered.set(i);
        }
      }

      CGAL_postcondition(face->info().covered.count() == 1);
    } while(++face != begin);

    begin = triangulation_.incident_faces(query_vertex_);
    face = begin;

    do
    {
      if(triangulation_.is_infinite(face))
      {
        continue;
      }

      for(int i = 0; i < 3; i++)
      {
        if(triangulation_.is_infinite(face->neighbor(i))
           || (face->neighbor(i)->info().covered.count() > 0))
        {
          continue;
        }

        create_subgraph(face->neighbor(i));
      }
    } while(++face != begin);
  }

  /// BFS over all faces to store how the edges of a face cover each other
  void create_subgraph(const Face_handle& start)
  {
    CGAL_precondition(triangulation_.is_valid(true));
    CGAL_precondition(start != Face_handle());

    std::deque<Face_handle> face_queue;
    face_queue.push_back(start);

    while(!face_queue.empty())
    {
      Face_handle face = face_queue.front();
      face_queue.pop_front();

      CGAL_precondition(!triangulation_.is_infinite(face));

      if(face->info().covered.count() > 0)
      {
        continue;
      }

      // orientation of query_point according to edge, RIGHT_TURN is outside of face
      CGAL::Orientation orientation[3];

      for(int i = 0; i < 3; i++)
      {
        orientation[i] = CGAL::orientation(face->vertex(triangulation_.ccw(i))->point(),
                                           face->vertex(triangulation_.cw(i))->point(),
                                           storage_.query_point());
      }

      for(int i = 0; i < 3; i++)
      {
        if(orientation[i] == CGAL::LEFT_TURN)
        {
          for(int j = 0; j < 3; j++)
          {
            if(i == j)
            {
              continue;
            }

            if(orientation[j] == CGAL::RIGHT_TURN)
            {
              // edge i (partially) covered by edge j
              face->info().covered.set(i);
            }
          }
        }
        else if(orientation[i] == CGAL::COLLINEAR)
        {
          // collinear edges can not cover anything
          face->info().covered.set(i);
        }
      }

      CGAL_precondition(face->info().covered.count() > 0);

      for(int i = 0; i < 3; i++)
      {
        if(triangulation_.is_infinite(face->neighbor(i))
           || (face->neighbor(i)->info().covered.count() > 0))
        {
          continue;
        }

        face_queue.push_back(face->neighbor(i));
      }
    }
  }

  /**
   * find the vertex ranks which belong to the triangulation edge and use the
   * storage to rank the edge
   */
  void rank_edge(const Face_handle& face, int edge)
  {
    // check if edge is part of input
    if(triangulation_.is_infinite(face, edge)
       || !face->is_constrained(edge))
    {
      // can't rank this! hammertime!
      return;
    }

    CGAL_precondition_msg(!face->info().ranked, "Face should not be ranked yet.");

    storage_.rank_edge(face->vertex(triangulation_.ccw(edge))->vertex_index(),
                       face->vertex(triangulation_.cw(edge))->vertex_index());
  }

  /// walk along the graph and rank the edges
  void rank_graph()
  {
    std::deque<Face_handle> face_queue;

    // rank faces around query point
    {
      typename Triangulation::Face_circulator begin = triangulation_.incident_faces(query_vertex_), face = begin;

      do
      {
        if(triangulation_.is_infinite(face))
        {
          continue;
        }

        CGAL_precondition(face->info().covered.count() == 1);
        CGAL_precondition(!face->info().ranked);

        // rank covering edges
        for(int i = 0; i < 3; ++i)
        {
          // check if edge i is covered
          if(face->info().covered.test(i))
          {
            if(!triangulation_.is_infinite(face->neighbor(i)))
            {
              face_queue.push_back(face->neighbor(i));
            }
          }
          // avoid duplicates
          else if(!face->neighbor(i)->info().ranked)
          {
            rank_edge(face, i);
          }
        }

        // rank edges on convex hull
        for(int i = 0; i < 3; ++i)
        {
          if(face->info().covered.test(i) && triangulation_.is_infinite(face->neighbor(i)))
          {
            rank_edge(face, i);
          }
        }

        face->info().ranked = true;
      } while(++face != begin);
    }

    while(!face_queue.empty())
    {
      Face_handle face = face_queue.front();
      face_queue.pop_front();

      CGAL_precondition(!triangulation_.is_infinite(face));

      // avoid ranking twice
      if(face->info().ranked)
      {
        continue;
      }

      debug("  current face: " + face_label(face));

      CGAL_precondition(face->info().covered.count() > 0);
      CGAL_precondition(face->info().covered.count() < 3);

      bool covering_faces_ranked = true;
      for(int i = 0; i < 3; ++i)
      {
        covering_faces_ranked &= face->info().covered.test(i)
            || face->neighbor(i)->info().ranked;
      }

      if(!covering_faces_ranked)
      {
        debug("  covering faces not ranked yet!");
        continue;
      }

      // rank covering edges and queue covered edges
      for(int i = 0; i < 3; ++i)
      {
        // check if edge i is covered
        if(face->info().covered.test(i))
        {
          if(!triangulation_.is_infinite(face->neighbor(i))
             && !face->neighbor(i)->info().ranked)
          {
            face_queue.push_back(face->neighbor(i));
          }
        }
        else
        {
          CGAL_precondition(!collinear_edge(face, i));
          rank_edge(face, i);
        }
      }

      for(int i = 0; i < 3; ++i)
      {
        if(face->info().covered.test(i))
        {
          // rank edges on convex hull
          if(triangulation_.is_infinite(face->neighbor(i)))
          {
            rank_edge(face, i);
          }
          // handle collinear edges last
          else if(collinear_edge(face, i)
                  && face->neighbor(i)->info().ranked)
          {
            rank_edge(face, i);
          }
        }
      }

      face->info().ranked = true;
    }
  }

  bool collinear_edge(const Face_handle& face, int edge) const
  {
    return (face->info().covered.test(edge)
        == face->neighbor(edge)->info().covered.test(
              triangulation_.mirror_edge(std::make_pair(face, edge)).second));
  }

  /// conversion from edge to Vertex_pair
  Vertex_pair vertex_pair(const Face_handle& face, int edge) const
  {
    return Vertex_pair(face->vertex(triangulation_.ccw(edge)), face->vertex(triangulation_.cw(edge)));
  }

#if ENABLE_VP_DEBUG
  /// \name labels
  //@{
  /// helper function to print an edge label
  std::string edge_label(const Face_handle& face, int edge) const
  {
    return edge_label(face->vertex(triangulation_.cw(edge)), face->vertex(triangulation_.ccw(edge)));
  }

  /// helper function to print an edge label
  std::string edge_label(const Vertex_handle& from, const Vertex_handle& to) const
  {
    return vertex_label(from) + "->" + vertex_label(to);
  }

  std::string face_label(const Face_handle& face) const
  {
    if(triangulation_.is_infinite(face))
    {
      return "inf";
    }

    return vertex_label(face->vertex(0)) + ", " + vertex_label(face->vertex(1)) + ", " + vertex_label(face->vertex(2));
  }

  const std::string& vertex_label(const Vertex_handle& vertex) const
  {
    CGAL_precondition(vertex->vertex_index().assigned());

    return storage_.vertex(vertex->vertex_index()).label;
  }

  //@}
#endif

};

#endif // ER_TRIANGULATION_H
