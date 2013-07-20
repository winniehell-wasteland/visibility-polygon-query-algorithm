#ifndef VISIBILITY_POLYGON_H
#define VISIBILITY_POLYGON_H

#include <vector>

#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/box_intersection_d.h>

// TODO: check what the hell is going on there
// needs to be after triangulation stuff
#include <CGAL/Arr_point_location/Arr_lm_grid_generator.h>
#include <CGAL/Arr_landmarks_point_location.h>

#include "edge_ranking.h"
#include "vertex_ranking.h"
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
          class InputContainer = CGAL::Arrangement_2<CGAL::Arr_linear_traits_2<Kernel> >,
          class OutputContainer  = CGAL::Arrangement_2<CGAL::Arr_linear_traits_2<Kernel> > >
class Visibility_polygon_generator
{
public:
  /// \name template parameters
  //@{
  typedef InputContainer   Input_container;
  typedef Kernel_          Kernel;
  typedef OutputContainer  Output_container;
  //@}

  typedef typename Input_container::Halfedge_const_handle  Halfedge_const_handle;

  /// \name storage definitions
  //@{
  typedef Visibility_storage<Kernel>         Storage;
  typedef typename Storage::Edge_index       Edge_index;
  typedef typename Edge_index::Index_type    Edge_index_type;
  typedef typename Storage::Vertex_index     Vertex_index;
  typedef typename Vertex_index::Index_type  Vertex_index_type;
  typedef typename Storage::Sector_index     Sector_index;
  typedef typename Sector_index::Index_type  Sector_index_type;
  typedef typename Storage::Edge_entry       Edge_entry;
  typedef typename Storage::Vertex_entry     Vertex_entry;
  //@}

  //typedef Er_Std_Sort<Kernel, Storage>       Edge_ranking;
  typedef Er_Triangulation<Kernel, Storage>  Edge_ranking;

  typedef Vr_Dual_plane<Kernel, Storage>  Vertex_ranking;
  //typedef Vr_Std_Sort<Kernel, Storage>  Vertex_ranking;

  /// default constructor
  Visibility_polygon_generator(const Input_container& input, bool allow_collinear = false) :
    input_(input),
    storage_(input_.number_of_edges(), input_.number_of_vertices(), allow_collinear),
    edge_ranking_(storage_),
    vertex_ranking_(storage_)
  {
    CGAL_precondition(!input_.is_empty());
    CGAL_precondition(input_.number_of_unbounded_faces() == 1);

    // TODO: handle isolated vertices
    CGAL_precondition(input_.number_of_isolated_vertices() == 0);

    // TODO: plane sweep
    typedef std::vector<CGAL::Point_2<Kernel> >  Input_points;
    Input_points input_points;

    for(typename Input_container::Vertex_const_iterator vertex = input_.vertices_begin(); vertex != input_.vertices_end(); ++vertex)
    {
      input_points.push_back(vertex->point());
    }

    std::sort(input_points.begin(), input_points.end(), typename Kernel::Less_y_2());

    // insert vertices into edge ranking and vertex ranking
    for(typename Input_points::const_iterator point = input_points.begin(); point != input_points.end(); ++point)
    {
      Vertex_index_type index = storage_.insert_input_vertex(*point);

      debugmsg( "Insert point at " << *point
#if ENABLE_VP_DEBUG
      << " - " << storage_.vertex(index).label
#endif
     );

      edge_ranking_.insert_input_vertex(index);
      vertex_ranking_.insert_input_vertex(index);

#if ENABLE_VP_DEBUG
      bbox_input_.extend((double[2]){CGAL::to_double(point->x()), CGAL::to_double(point->y())});
#endif
    }

    CGAL_postcondition(storage_.num_input_vertices() == input_.number_of_vertices());

    // insert edges into edge ranking
    for(typename Input_container::Edge_const_iterator edge = input_.edges_begin(); edge != input_.edges_end(); ++edge)
    {
      CGAL_precondition(!edge->is_fictitious());
      CGAL_precondition(!edge->source()->is_at_open_boundary()
                        && !edge->target()->is_at_open_boundary());

      storage_.insert_input_edge(edge->source()->point(), edge->target()->point());
      edge_ranking_.insert_input_edge(edge->source()->point(), edge->target()->point());
    }

#if ENABLE_VP_DEBUG
    {
      CairoImage img("triangulation", bbox_input_.bbox());
      edge_ranking_.draw(img);
    }
#endif
  }

  void query(const CGAL::Point_2<Kernel>& query_point, Output_container& output, Halfedge_const_handle incident_edge)
  {
    debugmsg( "Query point at (" << query_point << ")");
    if(incident_edge != input_.halfedges_end())
    {
      debugmsg( "  edge at (" << incident_edge->source()->point() << ")->("
                << incident_edge->target()->point() << ")" );
    }

#if ENABLE_VP_DEBUG
  ++visibilty_polygons_computed;
#endif

    CGAL_precondition_msg(incident_edge != Halfedge_const_handle(), "Incident edge can not be empty.");
    CGAL_precondition_msg(!incident_edge->is_fictitious() || (incident_edge == input_.halfedges_end()), "Incident edge can not be fictitious.");

    storage_.insert_query_point(query_point);

  #if ENABLE_VP_DEBUG
    {
      CairoImage img("input_" + tostr(visibilty_polygons_computed, 3), bbox_input_.bbox());
      storage_.draw_input(img);

      if(storage_.is_query_point_on_edge()
        || storage_.is_query_point_on_vertex())
      {
        CGAL_precondition(incident_edge != input_.halfedges_end());
        CGAL_precondition(!incident_edge->source()->is_at_open_boundary() && !incident_edge->target()->is_at_open_boundary());

        img.set_color(CGAL::YELLOW);
        img.draw_edge<Input_container>(*incident_edge);
        img.draw_point(incident_edge->source()->point());

        if(storage_.is_query_point_on_vertex())
        {
          img.draw_edge<Input_container>(*incident_edge->next());
        }
      }
    }
  #endif

    if(storage_.is_query_point_on_edge()
      || storage_.is_query_point_on_vertex())
    {
      CGAL_precondition(incident_edge != input_.halfedges_end());
    }
    else
    {
      CGAL_precondition(storage_.is_query_point_assigned());
      CGAL_precondition(incident_edge == input_.halfedges_end());
    }

    edge_ranking_.insert_query_point(query_point);

    // insert split vertices
    for(Vertex_index_type i = 0; i < storage_.num_split_vertices(); ++i)
    {
      edge_ranking_.insert_split_vertex(i);
    }

#if ENABLE_VP_DEBUG
    {
      CairoImage img("split_triangulation_" + tostr(visibilty_polygons_computed, 3), bbox_input_.bbox());
      edge_ranking_.draw(img);
    }
#endif

    vertex_ranking_.rank_vertices();
    debug(tostr(storage_.num_vertex_ranks()) + " vertices ranked");
    debug("");

  #if ENABLE_VP_DEBUG
    {
      CairoImage img("ranked_vertices_" + tostr(visibilty_polygons_computed, 3), bbox_input_.bbox());
      img.draw_arrangement(input_, CGAL::YELLOW, CGAL::BLACK);
      storage_.draw_vertex_ranks(img, query_point);
      img.draw_point(query_point, CGAL::GREEN);
    }
  #endif

    CGAL_postcondition_msg(storage_.all_vertices_ranked(), "All vertices must be ranked.");

    edge_ranking_.rank_edges();
    debug(tostr(storage_.num_edge_ranks()) + " edges ranked");
    debug("");

#if ENABLE_VP_DEBUG
    {
      CairoImage img("ranked_edges_" + tostr(visibilty_polygons_computed, 3), bbox_input_.bbox());
      storage_.draw_edge_ranks(img);
      img.draw_point(query_point, CGAL::BLUE);
    }
#endif

    CGAL_postcondition_msg(storage_.all_edges_ranked(), "All edges must be ranked.");

    {
      Sector_index skip_start, skip_end;

      CGAL_precondition_msg(!skip_start.assigned() && !skip_end.assigned(), "Skip range should initially be not assigned.");

      if(storage_.is_query_point_on_vertex())
      {
        CGAL_precondition(incident_edge->target()->point() == query_point);

        Vertex_index_type index_start = storage_.find_vertex(incident_edge->source()->point()),
            index_end = storage_.find_vertex(incident_edge->next()->target()->point());

        skip_start = storage_.sector_after(storage_.vertex(index_start).rank);
        skip_end = storage_.sector_after(storage_.vertex(index_end).rank);

        CGAL_postcondition(skip_start.assigned() && skip_end.assigned());
      }
      // query point lies on edge
      else if(storage_.is_query_point_on_edge())
      {
        Vertex_index_type index_start = storage_.find_vertex(incident_edge->source()->point()),
            index_end = storage_.find_vertex(incident_edge->target()->point());

        skip_start = storage_.sector_after(storage_.vertex(index_start).rank);
        skip_end = storage_.sector_after(storage_.vertex(index_end).rank);

        CGAL_postcondition(skip_start.assigned() && skip_end.assigned());
      }
      else
      {
        CGAL_precondition(incident_edge == input_.halfedges_end());
      }

      join_edges(output, skip_start, skip_end);
      debug("");
    }

      debug("VP has " + tostr(output.number_of_vertices() - output.number_of_vertices_at_infinity()) + " vertices");
      debug("\n");

      {
        CGAL_precondition(output.number_of_vertices() > 0);

        typename Output_container::Vertex_iterator
            vertex = output.vertices_begin(),
            first = vertex, next = vertex;
        ++next;

        while(next != output.vertices_end())
        {
          debug("  adding segment");
          CGAL_precondition_msg(vertex->point() != next->point(), "Segment is degenerated.");
          CGAL::insert_non_intersecting_curve(output, CGAL::Segment_2<Kernel>(vertex->point(), next->point()));

          vertex = next;
          ++next;
        }

        debug("  adding segment");
        CGAL_precondition_msg(vertex->point() != first->point(), "Segment is degenerated.");
        CGAL::insert_non_intersecting_curve(output, CGAL::Segment_2<Kernel>(vertex->point(), first->point()));
      }

  #if ENABLE_VP_DEBUG
    {
      CairoImage img("visibility_polygon_" + tostr(visibilty_polygons_computed, 3), bbox_input_.bbox());
      img.draw_arrangement(input_);

      bool contains_query_point = false;

      img.set_color(CGAL::RED);
      for(typename Output_container::Vertex_const_iterator vertex = output.vertices_begin();
          vertex != output.vertices_end();
          ++vertex)
      {
        if(vertex->point() == query_point)
        {
          contains_query_point = true;
        }
        else
        {
          img.draw_point(vertex->point());
        }
      }

      for(typename Output_container::Edge_const_iterator edge = output.edges_begin();
          edge != output.edges_end();
          ++edge)
      {
        img.draw_segment(edge->source()->point(), edge->target()->point());
      }

      if(contains_query_point)
      {
        img.draw_point(query_point, CGAL::YELLOW);
      }
      else
      {
        img.draw_point(query_point, CGAL::GREEN);
      }

      storage_.draw_labels(img);
    }
  #endif

    CGAL_postcondition( output.number_of_unbounded_faces() == 1 );
    CGAL_postcondition( output.number_of_faces() == 2 );

    edge_ranking_.clear();

#if ENABLE_VP_DEBUG
  {
    CairoImage img("reverted_triangulation_" + tostr(visibilty_polygons_computed, 3), bbox_input_.bbox());
    edge_ranking_.draw(img);
  }
#endif

    storage_.clear();
  }

private:
  /// store input polygon/arrangement
  const Input_container& input_;
  /// information for joining visible edges
  Storage storage_;
  /// triangulated input
  Edge_ranking edge_ranking_;
  /// ranking for arrangement vertices
  Vertex_ranking  vertex_ranking_;
#if ENABLE_VP_DEBUG
  CGAL::Box_intersection_d::Box_d<double, 2> bbox_input_;
#endif

  // can't copy this class! (not thread safe)
  Visibility_polygon_generator(const Visibility_polygon_generator& other) :
    input_(other.input_)
  {
    assert(0);
  }

  void add_to_output(Output_container& output, const CGAL::Point_2<Kernel>& point)
  {
    CGAL_precondition_msg(output.number_of_faces() == 1, "Visibility polygon can not have more than one face.");
    if((output.number_of_vertices() > 0)
       && ((--output.vertices_end())->point() == point))
    {
      debug("Duplicate vertex detected.");
    }

    output.insert_in_face_interior(point, output.unbounded_face());
  }

  void join_edges(Output_container& output, const Sector_index& skip_start, const Sector_index& skip_end)
  {
    // ensure either both or no skip parameter is set
    CGAL_precondition(skip_start.assigned() == skip_end.assigned());

    debug("JOIN EDGES");

    if(skip_start.assigned())
    {
      debug("skipping from " + storage_.sector_label(skip_start) + " to " + storage_.sector_label(skip_end));
    }

    // check sector mapping
    for(Vertex_index_type vertex = 0; vertex < storage_.num_vertex_ranks(); ++vertex)
    {
      //debug(storage_.ranked_vertex(vertex).label + " " + storage_.sector_label(storage_.sector_after(vertex)));
    }

    // loop over all edges with ascending rank (i.e. from visible to not visible)
    for(Edge_index_type edge_index = 0; edge_index < storage_.num_edge_ranks(); ++edge_index)
    {
      const Edge_entry& edge = storage_.ranked_edge(edge_index);

      CGAL_precondition(edge.source.assigned() && edge.target.assigned());

      debug("  current edge: " + storage_.ranked_vertex(edge.source).label
            + "->" + storage_.ranked_vertex(edge.target).label + " ("
            + tostr(edge_index) + ")");

      Sector_index_type sector = storage_.sector_after(edge.source);

#if ENABLE_VP_DEBUG
      std::cout << "  current set: ";
      storage_.dump_set(sector);
#endif

      sector = storage_.find_last(sector);

      Vertex_index_type end = edge.target;

      if(!storage_.allow_collinear())
      {
        while((end < storage_.num_vertex_ranks() - 1)
              && CGAL::collinear(storage_.query_point(),
                                 storage_.ranked_vertex(end).point,
                                 storage_.ranked_vertex(end+1).point))
        {
          ++end;
        }
      }

      while(sector < storage_.sector_after(end))
      {
        debug("  sector " + storage_.sector_label(sector) + " covered");

        CGAL_precondition(!storage_.sector(sector).visible_edge.assigned());
        storage_.sector(sector).visible_edge = edge_index;
        storage_.link_set(sector);
        sector = storage_.find_last(sector);

#if ENABLE_VP_DEBUG
      std::cout << "  new set: ";
      storage_.dump_set(sector);
#endif
      }

  #if ENABLE_VP_DEBUG
      std::cout << std::endl;
  #endif

    }

#if ENABLE_VP_DEBUG
    {
      CairoImage img("sectors_" + tostr(visibilty_polygons_computed, 3), bbox_input_.bbox());
      img.draw_arrangement(input_);
      storage_.draw_sectors(img);
      img.draw_point(storage_.query_point(), CGAL::GREEN);
      storage_.draw_labels(img);
    }
#endif

#if ENABLE_VP_STORAGE_DEBUG
    storage_.dump_sectors();
#endif

    // ensure all sectors are covered
    CGAL_postcondition_msg(storage_.count_sets(skip_end, skip_start) == 1,
                           "At least one sector has no visible edge.");

    debug("Calculating visible edge segments...");

    Edge_index first_edge, last_edge;
    Vertex_index first_vertex, last_vertex;

    for(Sector_index_type sector = 0; sector < storage_.num_sectors(); ++sector)
    {
#if ENABLE_VP_DEBUG
      std::cout << storage_.sector_label(sector);
#endif

      if(storage_.in_range(sector, skip_start, skip_end))
      {
#if ENABLE_VP_DEBUG
        std::cout << ":\tskip" << std::endl;
#endif

        if((skip_start < skip_end) || !storage_.in_range(sector, 0, skip_end))
        {

          if(!storage_.same_vertex(first_vertex, storage_.vertex_before(skip_start))
             && !storage_.same_vertex(last_vertex, storage_.vertex_before(skip_start)))
          {
            last_vertex = storage_.vertex_before(skip_start);

            if(storage_.ranked_vertex(last_vertex).is_input_vertex())
            {
              debug("  adding " + storage_.ranked_vertex(last_vertex).label + " (skip first)");
              add_to_output(output, storage_.ranked_vertex(last_vertex).point);
            }

            if(!first_vertex.assigned())
            {
              first_vertex = last_vertex;
            }
          }

          if(storage_.is_query_point_on_vertex())
          {
            debug("  adding query point");
            add_to_output(output, storage_.query_point());
          }
        }

        if((skip_start > skip_end) && !storage_.in_range(sector, 0, skip_end))
        {
          // done
          break;
        }

        last_edge.reset();

        CGAL_precondition(!storage_.same_vertex(last_vertex, storage_.vertex_before(skip_end)));
        last_vertex = storage_.vertex_before(skip_end);

        if(storage_.ranked_vertex(last_vertex).is_input_vertex())
        {
          debug("  adding " + storage_.ranked_vertex(last_vertex).label + " (skip second)");
          add_to_output(output, storage_.ranked_vertex(last_vertex).point);

          if(!first_vertex.assigned())
          {
            first_vertex = last_vertex;
          }
        }

        if(sector < skip_end)
        {
          sector = skip_end - 1;
        }

        continue;
      }

      const Edge_index& visible_edge = storage_.sector(sector).visible_edge;
      CGAL_precondition_msg(visible_edge.assigned(), "No visible edge assigned to sector.");

      Vertex_index from = storage_.ranked_edge(visible_edge).source,
          to = storage_.ranked_edge(visible_edge).target;
      CGAL_precondition_msg(from.assigned() && to.assigned(), "Both vertices of visible edge need to be assigned.");

      Sector_index_type last_covered = sector;

      while((last_covered+1 < storage_.num_sectors())
            && storage_.sector(last_covered+1).visible_edge.assigned()
            && (storage_.sector(last_covered+1).visible_edge == visible_edge))
      {
        ++last_covered;
      }

#if ENABLE_VP_DEBUG
      if(last_covered > sector)
      {
        std::cout << " - " << storage_.sector_label(last_covered);
      }

      std::cout << ":\t" << storage_.ranked_vertex(from).label
                << "->" << storage_.ranked_vertex(to).label
                << std::endl;
#endif

      // edge source is visible
      if(sector == storage_.sector_after(from))
      {
        const Vertex_entry& from_vertex = storage_.ranked_vertex(from);

        if(from_vertex.is_input_vertex()
           && !storage_.same_vertex(first_vertex, from)
           && !storage_.same_vertex(last_vertex, from))
        {
          debug("  adding " + from_vertex.label + " (source)");
          add_to_output(output, from_vertex.point);

          last_vertex = from;

          if(!first_vertex.assigned())
          {
            first_vertex = from;
          }
        }
      }
      else if(last_vertex.assigned())
      {
        CGAL_precondition(!last_edge.assigned() || (last_edge < visible_edge));

        split_edge(output, visible_edge, last_vertex, last_vertex);
      }

      // edge target is visible
      if(storage_.sector_after(to) <= last_covered + 1)
      {
        const Vertex_entry& to_vertex = storage_.ranked_vertex(to);

        CGAL_precondition(!to_vertex.is_split_vertex()
                          || storage_.is_second_entry(to));

        if(to_vertex.is_input_vertex()
           && !storage_.same_vertex(first_vertex, to)
           && !storage_.same_vertex(last_vertex, to))
        {
          debug("  adding " + to_vertex.label + " (target)");
          add_to_output(output, to_vertex.point);

          last_vertex = to;

          if(!first_vertex.assigned())
          {
            first_vertex = to;
          }
        }
      }
      else if(last_covered + 1 < storage_.num_sectors())
      {
        const Edge_index& next_edge = storage_.sector(last_covered + 1).visible_edge;

        split_edge(output, visible_edge, storage_.ranked_edge(next_edge).source, last_vertex);
      }

      last_edge = visible_edge;
      if(!first_edge.assigned())
      {
        first_edge = last_edge;
      }

      sector = last_covered;
    }

    CGAL_precondition(last_vertex.assigned() && first_vertex.assigned());
    CGAL_precondition(last_edge.assigned() && first_edge.assigned());

    Vertex_index first_source = storage_.ranked_edge(first_edge).source,
        last_target = storage_.ranked_edge(last_edge).target;

    // handle overlap
    if(!storage_.same_vertex(last_target, first_source))
    {
      debug("first vertex: "+storage_.ranked_vertex(first_vertex).label);
      debug("first source: "+storage_.ranked_vertex(first_source).label);
      debug("last vertex: "+storage_.ranked_vertex(last_vertex).label);
      debug("last target: "+storage_.ranked_vertex(last_target).label);

      if(!storage_.same_vertex(first_source, first_vertex)
         && !storage_.same_vertex(first_source, last_vertex)
         && storage_.ranked_vertex(first_source).is_split_vertex())
      {
        const Vertex_entry& vertex = storage_.ranked_vertex(first_source);

        debug("  adding " + vertex.label + " (first source)");

        //CGAL_precondition(storage_.same_vertex(last_vertex, last_target));
        //CGAL_precondition(vertex.is_split_vertex());

        add_to_output(output, vertex.point);
      }
      else if(!storage_.same_vertex(last_target, first_vertex)
              && !storage_.same_vertex(last_target, last_vertex)
              && storage_.ranked_vertex(last_target).is_split_vertex())
      {
        //CGAL_precondition(storage_.same_vertex(first_vertex, first_source)
        //                  || storage_.is_query_point_on_vertex());

        const Vertex_entry& vertex = storage_.ranked_vertex(last_target);

        debug("  adding " + vertex.label + " (last target)");

        CGAL_precondition(vertex.is_split_vertex());
        add_to_output(output, vertex.point);
      }
    }
  }

  void split_edge(Output_container& output,
                  const Edge_index& edge_index,
                  const Vertex_index& split_pos,
                  Vertex_index& last_vertex)
  {
    const Edge_entry& edge_entry = storage_.ranked_edge(edge_index);
    const Vertex_index& from = edge_entry.source, to = edge_entry.target;

    CGAL_precondition(from.assigned() && to.assigned());
    CGAL_precondition(split_pos.assigned());

    if(storage_.same_vertex(from, split_pos)
       || storage_.same_vertex(to, split_pos))
    {
      // nothing to split - just a yellow lemon tree
      return;
    }

    debug("  splitting " + storage_.ranked_vertex(edge_entry.source).label
          + "->" + storage_.ranked_vertex(edge_entry.target).label
          + " at " + storage_.ranked_vertex(split_pos).label);

    CGAL::Line_2<Kernel> edge(storage_.ranked_vertex(from).point, storage_.ranked_vertex(to).point);
    CGAL::Ray_2<Kernel> split_ray(storage_.query_point(), storage_.ranked_vertex(split_pos).point);

    CGAL_precondition(!split_ray.is_degenerate());
    CGAL_precondition(!CGAL::parallel(edge, split_ray.supporting_line()));

    if(split_ray.has_on(storage_.ranked_vertex(from).point))
    {
      const Vertex_entry& vertex = storage_.ranked_vertex(from);

      if(!storage_.same_vertex(last_vertex, from))
      {
        debug("  adding " + vertex.label + " (split source)");
        add_to_output(output, vertex.point);

        last_vertex = from;
      }
    }
    else if(split_ray.has_on(storage_.ranked_vertex(to).point))
    {
      CGAL_precondition(!storage_.same_vertex(last_vertex, to));

      debug("  adding " + storage_.ranked_vertex(to).label + " (split target)");
      add_to_output(output, storage_.ranked_vertex(to).point);

      last_vertex = to;
    }
    else
    {
      debug("  adding ... (split)");

      CGAL::Object split_point = CGAL::intersection(edge, split_ray);
      CGAL_precondition_msg(split_point.is<CGAL::Point_2<Kernel> >(), "Intersection is not a point.");

      CGAL_precondition(!last_vertex.assigned()
                        || (storage_.ranked_vertex(last_vertex).point
                            != CGAL::object_cast<CGAL::Point_2<Kernel> >(split_point)));
      add_to_output(output, CGAL::object_cast<CGAL::Point_2<Kernel> >(split_point));
    }
  }
};

#endif // VISIBILITY_POLYGON_H
