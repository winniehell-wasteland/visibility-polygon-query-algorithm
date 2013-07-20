#ifndef VR_DUAL_PLANE_H
#define VR_DUAL_PLANE_H

#include <map>
#include <utility>
#include <vector>

#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_curve_data_traits_2.h>
#include <CGAL/Arr_linear_traits_2.h>

#include "log.h"

#if ENABLE_VP_VERTEX_RANKING_DEBUG
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arr_observer.h>
#endif

#undef debug
#if ENABLE_VP_VERTEX_RANKING_DEBUG
  #define debug(msg) \
  { \
    std::cout << (msg) << std::endl; \
  }
#else
  #define debug(msg)
#endif

/// arrangement of dual lines for input vertices for efficient ordering by slope
template <class Kernel_,
          class Storage_>
class Vr_Dual_plane
{
public:
  /// \name template parameters
  //@{
  typedef Kernel_        Kernel;
  typedef Storage_       Storage;
  //@}

private:
  typedef typename Storage::Size          Size;
  typedef typename Storage::Vertex_index  Vertex_index;
  typedef std::vector<Vertex_index>       Ranked_vertices;

#if ENABLE_VP_VERTEX_RANKING_DEBUG
  struct Dual_object_data
  {
    Dual_object_data() :
      label("#")
    {

    }

    std::string label;
  };
#endif

  /// \name kernel types
  //@{
  typedef typename Kernel::Equal_x_2    Equal_x_2;
  typedef CGAL::Line_2<Kernel>          Line_2;
  typedef CGAL::Point_2<Kernel>         Point_2;
  //@}

  /// \name dual plane types
  //@{
  typedef CGAL::Arr_linear_traits_2<Kernel>                Linear_traits;
  typedef CGAL::Arr_curve_data_traits_2<Linear_traits,
                                        Vertex_index>      Curve_traits;
  typedef typename Curve_traits::X_monotone_curve_2        Curve_2;
  typedef CGAL::Arr_vertex_base<Point_2>                   Vertex_base;
  typedef CGAL::Arr_halfedge_base<Curve_2>                 Halfedge_base;
#if ENABLE_VP_VERTEX_RANKING_DEBUG
  typedef CGAL::Arr_extended_vertex<Vertex_base,
                                    Dual_object_data>      Extended_vertex;
  typedef CGAL::Arr_extended_halfedge<Halfedge_base,
                                      Dual_object_data>    Extended_halfedge;
  typedef CGAL::Arr_dcel_base<Extended_vertex,
                              Extended_halfedge,
                              CGAL::Arr_face_base>         Dcel;
#else
  typedef CGAL::Arr_dcel_base<Vertex_base,
                              Halfedge_base,
                              CGAL::Arr_face_base>         Dcel;
#endif
  typedef CGAL::Arrangement_2<Curve_traits, Dcel>          Dual_plane;

  typedef typename Dual_plane::
      Ccb_halfedge_const_circulator               Ccb_halfedge_const_circulator;
  typedef typename Dual_plane::Halfedge_handle    Halfedge_handle;
  //@}

public:
  /// default constructor
  Vr_Dual_plane(Storage& storage) :
    dual_plane_(),
    storage_(storage)
#if ENABLE_VP_VERTEX_RANKING_DEBUG
    , bbox_(),
    debug_info_created_(false)
#endif
  {

  }

  /**
    * create a line in the dual plane with equation y = v.x()*x + v.y() for
    * the input vertex with given storage index
    */
  void insert_input_vertex(const Vertex_index& vertex)
  {
    // make sure, it fits
    CGAL_precondition(Size(vertex) == Vertex_index(vertex));
    CGAL_precondition(Size(vertex) == typename Ranked_vertices::size_type(vertex));

    const Point_2& point = storage_.vertex(vertex).point;
    Line_2 line(-point.x(), 1, -point.y());
    CGAL::insert(dual_plane_, Curve_2(line, vertex));

    CGAL_postcondition(dual_plane_.is_valid());
  }

  void rank_vertices()
  {
    CGAL_precondition(dual_plane_.is_valid());
    CGAL_precondition(!dual_plane_.is_empty());
    CGAL_precondition(dual_plane_.number_of_vertices() > 0);

    debug("RANK VERTICES");

#if ENABLE_VP_VERTEX_RANKING_DEBUG
    init_debug_info();
#endif

    Line_2 query_line(-storage_.query_point().x(), 1, -storage_.query_point().y());

    // edge that is currently cut by sweep line and the edge that led into current face
    Halfedge_handle sweep_edge, entering_edge;

  #if ENABLE_VP_VERTEX_RANKING_DEBUG
    {
      CairoImage img("dual_" + tostr(visibilty_polygons_computed, 3),
                     bbox_, DUAL_SCALE_X, DUAL_SCALE_Y);

      img.draw_dual_plane(dual_plane_, bbox_);
      img.draw_dual_line(query_line, bbox_, CGAL::GREEN);
    }
  #endif

    find_first_edge(query_line, sweep_edge);

    bool done = false;
    Ranked_vertices sorted_vertices;

#if ENABLE_VP_VERTEX_RANKING_DEBUG
    print_face(sweep_edge->face());
#endif

    entering_edge = sweep_edge;

    while(!done)
    {
      if(!sweep_edge->is_fictitious())
      {
        bool edge_lr = (sweep_edge->direction() == CGAL::ARR_LEFT_TO_RIGHT);

        debug("    " + sweep_edge->data().label + ": " + (edge_lr?"LR":"RL"));

        Line_2 line = sweep_edge->curve().supporting_line();
        bool line_below = (CGAL::compare_slopes(line, query_line)
                           == CGAL::SMALLER);

        bool intersects = false;

        // line parallel to query line
        if(CGAL::parallel(line, query_line))
        {
          intersects = false;
        }
        else if(sweep_edge->source()->is_at_open_boundary())
        {
          const Point_2& target = sweep_edge->target()->point();

          if(query_line.has_on(target))
          {
            debug("T=O");
          }
          else if(query_line.has_on_negative_side(target))
          {
            debug("T=R");
          }
          else
          {
            debug("T=L");
          }

          intersects = !query_line.has_on_positive_side(target)
              && (edge_lr == line_below);
        }
        else if(sweep_edge->target()->is_at_open_boundary())
        {
          const Point_2& source = sweep_edge->source()->point();

          if(query_line.has_on(source))
          {
            debug("S=O");
          }
          else if(query_line.has_on_negative_side(source))
          {
            debug("S=R");
          }
          else
          {
            debug("S=L");
          }

          intersects = !query_line.has_on_negative_side(source)
              && (edge_lr == line_below);

        }
        else if(query_line.has_on_positive_side(sweep_edge->source()->point()))
        {
          intersects = !query_line.has_on_positive_side(sweep_edge->target()->point());
        }
        else if(query_line.has_on(sweep_edge->source()->point()))
        {
          intersects = query_line.has_on_negative_side(sweep_edge->target()->point());
        }
        else
        {
          intersects = false;
        }

        // do not cross query line
        CGAL_precondition((line != query_line) || !intersects);

        // found edge intersecting query line
        if(intersects)
        {
          sorted_vertices.push_back(sweep_edge->curve().data());

          debug("  leaving through " + sweep_edge->data().label);

          sweep_edge = sweep_edge->twin();
          entering_edge = sweep_edge;

  #if ENABLE_VP_VERTEX_RANKING_DEBUG
          print_face(sweep_edge->face());
  #endif
        }
        else
        {
          debug("  skipping " + sweep_edge->data().label);
        }
      }

      sweep_edge = sweep_edge->prev();
      done = (sweep_edge == entering_edge);
    }

    for(Size i = 0; i < sorted_vertices.size(); ++i)
    {
      if(!storage_.allow_collinear())
      {
        Size collinear_end = i;

        if(collinear_end < sorted_vertices.size() - 1)
        {
          debug("Testing "+storage_.vertex(sorted_vertices[i]).label
                +" and "+storage_.vertex(sorted_vertices[collinear_end+1]).label);
        }

        const CGAL::Point_2<Kernel>& current_point =
            storage_.vertex(sorted_vertices[i]).point;

        while((collinear_end < sorted_vertices.size() - 1)
              && CGAL::collinear(storage_.query_point(), current_point,
                                 storage_.vertex(sorted_vertices[collinear_end+1]).point)
              && CGAL::Segment_2<Kernel>(storage_.query_point(), current_point)
                  .has_on(storage_.vertex(sorted_vertices[collinear_end+1]).point)
              )
        {
          debug("  collinear");
          ++collinear_end;

          if(collinear_end < sorted_vertices.size() - 1)
          {
            debug("Testing "+storage_.vertex(sorted_vertices[i]).label
                  +" and "+storage_.vertex(sorted_vertices[collinear_end+1]).label);
          }
        }

        Size k = collinear_end;
        while(i < k)
        {
          debug("Swap "+storage_.vertex(sorted_vertices[i]).label
                +" and "+storage_.vertex(sorted_vertices[k]).label);

          std::swap(sorted_vertices[i], sorted_vertices[k]);
          storage_.rank_vertex(sorted_vertices[i]);

          ++i;
          --k;
        }

        if(i != collinear_end)
        {
          debug("");
        }

        while(i < collinear_end)
        {
          storage_.rank_vertex(sorted_vertices[i]);
          ++i;
        }
      }

      storage_.rank_vertex(sorted_vertices[i]);
    }

    rank_parallels();
  }
private:
  Dual_plane        dual_plane_;
  Storage&          storage_;
#if ENABLE_VP_VERTEX_RANKING_DEBUG
  CGAL::Bbox_2      bbox_;
  bool              debug_info_created_;
#endif

  /// walk around the arrangement until query_line is intersected
  void find_first_edge(const Line_2& query_line, Halfedge_handle& first_edge)
  {
    // ensure there is only one face inside the fictitious face
    CGAL_precondition(std::distance(dual_plane_.fictitious_face()->holes_begin(),
                                    dual_plane_.fictitious_face()->holes_end()) == 1);

    // current fictitious edges
    Ccb_halfedge_const_circulator
        begin = *dual_plane_.fictitious_face()->holes_begin(),
        boundary_edge = begin;

    do
    {
      CGAL_precondition(boundary_edge->is_fictitious());

      // left and right neighbors of fictitious edge
      Ccb_halfedge_const_circulator left = boundary_edge->twin(), right = left;

      while(left->is_fictitious())
      {
        left = left->next();
        CGAL_postcondition(left != right);
      }

      while(right->is_fictitious())
      {
        right = right->prev();
        CGAL_postcondition(right != left);
      }

      debug("checking " + left->data().label + " and " + right->data().label);

      CGAL_precondition(!left->is_fictitious()
                        && left->source()->is_at_open_boundary());
      CGAL_precondition(!right->is_fictitious()
                        && right->target()->is_at_open_boundary());

      const Line_2& left_line = left->curve().supporting_line();
      const Line_2& right_line = right->curve().supporting_line();

      Point_2 left_target = left->target()->point() + query_line.direction().vector();

      // true if left edge ist above query line
      bool left_above = (left_line.has_on_negative_side(left_target) && (left->direction() == CGAL::ARR_RIGHT_TO_LEFT))
          || (left_line.has_on_positive_side(left_target) && (left->direction() == CGAL::ARR_LEFT_TO_RIGHT))
          || (left_line.has_on(left_target) && (left_line.y_at_x(0) > query_line.y_at_x(0)));

      if(left_above)
      {
        if(right_line == query_line)
        {
          first_edge = dual_plane_.non_const_handle(right->prev());
          debug("starting with " + first_edge->data().label + " on query line "
                + storage_.vertex(first_edge->curve().data()).label);

          return;
        }
        else
        {
          Point_2 right_source = right->source()->point() + query_line.direction().vector();

          // true if right edge ist below query line
          bool right_below = (right_line.has_on_positive_side(right_source) && (right->direction() == CGAL::ARR_LEFT_TO_RIGHT))
              || (right_line.has_on_negative_side(right_source) && (right->direction() == CGAL::ARR_RIGHT_TO_LEFT))
              || (right_line.has_on(right_source) && (right_line.y_at_x(0) <= query_line.y_at_x(0)));

          if(right_below)
          {
            first_edge = dual_plane_.non_const_handle(right);
            debug("starting with " + first_edge->data().label + " on "
                  + storage_.vertex(first_edge->curve().data()).label);

            return;
          }
        }
      }
    } while(++boundary_edge != begin);

    CGAL_error_msg("No edge intersects query line");
  }

  /**
   * find all vertices with same x-coordinate as query point - i.e. those which
   * correspond to a dual line that is parallel to the query line
   */
  void rank_parallels()
  {
    Equal_x_2 equal_x;

    for(Size i = 0; i < storage_.num_input_vertices(); ++i)
    {
      const Point_2& point = storage_.vertex(i).point;

      if(!storage_.vertex(i).is_split_vertex()
         && equal_x(storage_.query_point(), point)
         && (!storage_.is_query_point_on_vertex()
             || (point != storage_.query_point())))
      {
        CGAL_precondition(CGAL::less_y(point, storage_.query_point()));

        storage_.rank_vertex(i);
      }
    }
  }

#if ENABLE_VP_VERTEX_RANKING_DEBUG
  void init_debug_info()
  {
    if(debug_info_created_)
    {
      return;
    }

    debug_info_created_ = true;

    CGAL::Box_intersection_d::Box_d<double, 2> tmp_bbox;

    // label vertices and extend bbox
    for(typename Dual_plane::Vertex_iterator vertex = dual_plane_.vertices_begin(); vertex != dual_plane_.vertices_end(); ++vertex)
    {
      tmp_bbox.extend((double[2]){CGAL::to_double(vertex->point().x()), CGAL::to_double(vertex->point().y())});

      std::string first, current, label;

      typename Dual_plane::Halfedge_around_vertex_const_circulator begin = vertex->incident_halfedges(),
          edge = begin;

      // label each vertex after its incident lines
      do
      {
        CGAL_precondition(edge->target()->point() == vertex->point());

        if(!edge->target()->is_at_open_boundary())
        {
          current = storage_.vertex(edge->curve().data()).label;

          // stop when half-way around to prevent duplicates
          if(current == first)
          {
            break;
          }

          if(first.length() == 0)
          {
            first = current;
            label = current;
          }
          else
          {
            label = label + ", " + current;
          }
        }
      } while(++edge != begin);

      vertex->data().label = label;
    }

    bbox_ = tmp_bbox.bbox();

    // label edges after there vertices
    for(typename Dual_plane::Halfedge_iterator edge = dual_plane_.halfedges_begin(); edge != dual_plane_.halfedges_end(); ++edge)
    {
      if(edge->is_fictitious())
      {
        continue;
      }

      std::string label = "(";

      if(edge->source()->is_at_open_boundary())
      {
        label += "#";
        label += storage_.vertex(edge->curve().data()).label;
      }
      else
      {
        label += edge->source()->data().label;
      }

      label += ") -> (";

      if(edge->target()->is_at_open_boundary())
      {
        label += storage_.vertex(edge->curve().data()).label;
        label += "#";
      }
      else
      {
        label += edge->target()->data().label;
      }

      label += ")";

      edge->data().label = label;
    }
  }

  void print_face(const typename Dual_plane::Face_const_handle& face) const
  {
    Ccb_halfedge_const_circulator begin = face->outer_ccb(),
    edge = begin;

    std::cout << "face: ";

    do
    {
      std::cout << edge->data().label << ", ";
    } while(++edge != begin);

    std::cout << std::endl;
  }
#endif
};

#endif // VR_DUAL_PLANE_H
