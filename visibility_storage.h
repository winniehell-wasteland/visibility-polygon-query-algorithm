#ifndef VISIBILITY_STORAGE_H
#define VISIBILITY_STORAGE_H

#include <vector>
#include <limits>

#include <CGAL/Cartesian_d.h>
#include <CGAL/predicates_d.h>

#include "visibility_index.h"
#include "visibility_set.h"

#include "log.h"

#undef debug
#if ENABLE_VP_STORAGE_DEBUG
  #define debug(msg) \
  { \
    std::cout << (msg) << std::endl; \
  }
#else
  #define debug(msg)
#endif

/**
 * storage for visibility polygon calculation
 *
 * \tparam Kernel_
 * \tparam Size_
 */
template <class Kernel_,
          typename Size_ = std::size_t>
class Visibility_storage
{
public:
  /// \name template parameters
  //@{
  typedef Kernel_  Kernel;
  typedef Size_    Size;
  //@}

  /// \name indices
  //@{
  typedef Visibility_index<Size>  Edge_index;
  typedef Visibility_index<Size>  Vertex_index;
  typedef Visibility_index<Size>  Sector_index;
  //@}

  /**
   * Edge stored by the vertex indices of its end points. Note that the
   */
  struct Edge_entry
  {
    /// default constructor
    Edge_entry() :
      source(),
      target()
    {

    }

    /// \return true if no data has been assigned to this edge entry
    bool empty() const
    {
      return !source.assigned() && !target.assigned();
    }

    /// set entry back to initial state
    void reset()
    {
      source.reset();
      target.reset();
    }

    /**
     * \name vertices
     * end vertices of this edge such that visible face is to the left of
     * the line from source to target
     */
    //@{
    Vertex_index source, target;
    //@}
  };

  // TODO: can we use boost::disjoint_sets here?
  /// \private
  typedef Visibility_set<Size>  Set;

  struct Sector_entry
  {
    /// default constructor
    Sector_entry() :
      set(),
      visible_edge()
    {

    }

    /// \return true if no data has been assigned to this visibility entry
    bool empty() const
    {
      return !visible_edge.assigned();
    }

    /// set this entry belongs to
    Set set;

    /**
     * stores which edge covers this entry
     */
    Edge_index visible_edge;
  };

  struct Vertex_entry
  {
    enum Vertex_type
    {
      /// not assigned
      UNDEFINED = 0,
      /// input vertex
      INPUT,
      /// input vertex on split ray
      FAKE_SPLIT,
      /// split vertex
      SPLIT
    };

    /// default constructor
    Vertex_entry() :
      point(CGAL::ORIGIN),
      rank(),
      type(UNDEFINED)
#if ENABLE_VP_DEBUG
    , label("#")
#endif
    {

    }

    /// \return true if this entry represents an input vertex
    bool is_input_vertex() const
    {
      return (type == INPUT) || (type == FAKE_SPLIT);
    }

    /// \return true if this entry represents a split vertex
    bool is_split_vertex() const
    {
      return (type == SPLIT) || (type == FAKE_SPLIT);
    }

    CGAL::Point_2<Kernel>  point;
    Vertex_index           rank;
    Vertex_type            type;

#if ENABLE_VP_DEBUG
    std::string            label;
#endif
  };

private:
  /// values for position of query point in input arrangement
  enum Query_point_pos
  {
    UNDEFINED = 0,
    ON_VERTEX,
    ON_EDGE,
    IN_FACE
  };

  /// \name container classes
  //@{
  typedef std::vector<Edge_entry>                Edge_entries;
  typedef std::vector<Edge_entry>                Edge_ranks;
  typedef std::vector<Sector_entry>              Sector_entries;
  typedef std::vector<Vertex_entry>              Vertex_entries;
  typedef std::map<CGAL::Point_2<Kernel>, Size>  Vertex_map;
  typedef std::vector<Vertex_index>              Vertex_ranks;
  //@}
public:
  // TODO: can we use less space?
  /**
   * default constructor
   * \param num_edges number of edges to store
   * \param num_vertices number of vertices to store
   */
  Visibility_storage(const Size& edges,
                     const Size& vertices,
                     bool allow_collinear) :
    allow_collinear_(allow_collinear),
    // we need additional space for split edges
    edge_entries_(2*edges),
    edge_ranks_(2*edges),
    sector_entries_(vertices + 2),
    // we need additional space for split vertices
    vertex_entries_(vertices + edges),
    vertex_map_(),
    // we need additional space for split vertices
    vertex_ranks_(vertices + 2*edges),
    num_skipped_edges_(0),
    num_split_vertices_(0),
    num_split_vertex_entries_(0),
    num_input_edges_(0),
    num_input_vertices_(0),
    edge_rank_(0),
    vertex_rank_left_(0),
    vertex_rank_right_(0),
    query_point_(CGAL::ORIGIN),
    query_point_pos_(UNDEFINED)
  {
    CGAL_precondition(edges > 0);
    CGAL_precondition(vertices > 0);

    // ensure that every data structure is large enough
    Size max_value = std::numeric_limits<Size>::max();

    // ignore if it is not used
    (void) max_value;

    // 2*edges fit
    CGAL_precondition(edges <= max_value/2);
    // vertices + 2 fit
    CGAL_precondition((2 <= max_value) && (vertices <= max_value - 2));
    // vertices + 2*edges fit
    CGAL_precondition(vertices/2 <= max_value/2 - edges);

    CGAL_precondition(edge_entries_.max_size()/2 >= edges);
    CGAL_precondition(edge_ranks_.max_size()/2 >= edges);

    CGAL_precondition(vertex_entries_.max_size() >= vertices);
    CGAL_precondition(vertex_entries_.max_size() - vertices >= edges);

    CGAL_precondition(vertex_map_.max_size() >= vertices);

    CGAL_precondition(vertex_ranks_.max_size() >= vertices);
    CGAL_precondition(vertex_ranks_.max_size()/2 - vertices >= edges);

    CGAL_precondition(sector_entries_.max_size() >= 2);
    CGAL_precondition(sector_entries_.max_size() - 2 >= vertices);

    // check if every thing is alrighty with the counts
    CGAL_precondition(num_edge_ranks() == 0);
    CGAL_precondition(num_input_edges() == 0);
    CGAL_precondition(num_sectors() == vertices + 1);
    CGAL_precondition(num_split_vertices() == 0);
    CGAL_precondition(num_vertex_ranks() == 0);
    CGAL_precondition(num_input_vertices() == 0);

    // initialize sets
    for(Size sector = 0; sector < sector_entries_.size(); ++sector)
    {
      sector_entries_[sector].set.last = sector;
    }
  }

  /// \name validation methods
  //@{
  /// \return true if all edges are ranked (i.e. the last used rank is the number of edges)
  bool all_edges_ranked() const
  {
    return (edge_rank_ == num_edge_ranks());
  }

  /// \return true if all vertices are ranked (i.e. the last used rank is the number of vertices)
  bool all_vertices_ranked() const
  {
    return (vertex_rank_right_ == num_vertex_ranks() - vertex_shift());
  }

  /**
   * counts the sets from start (inclusively) to end (exclusively)
   * or - if start > end - the sets until end and after start
   */
  Size count_sets(Sector_index start, Sector_index end) const
  {
    if(!start.assigned())
    {
      start = 0;
    }

    if(!end.assigned())
    {
      end = num_sectors();
    }

    CGAL_precondition((0 <= start) && (end <= num_sectors()));

    if(start > end)
    {
      Size sets = 1;

      if(end > 0)
      {
        sets = count_sets(0, end);
      }

      sets += count_sets(start, num_sectors()) - 1;

      return sets;
    }
    else
    {
      Size sets = 1, index = set(find_set_index(start)).last;

      while(index < end)
      {
        debug(sector_label(index) + " not covered");
        ++sets;
        index = set(find_set_index(index+1)).last;
      }

      return sets;
    }
  }
  //@}

  /// empty all data
  void clear()
  {
    CGAL_precondition(all_edges_ranked());
    CGAL_precondition(all_vertices_ranked());

    // remove edge ranks
    for(Size i = 0; i < num_edge_ranks(); ++i)
    {
      edge_ranks_[i].reset();
    }
    edge_rank_ = 0;

    // remove vertex ranks
    for(Size i = 0; i < num_input_vertices_; ++i)
    {
      CGAL_precondition(vertex_entries_[i].is_input_vertex());
      vertex_entries_[i].rank.reset();
      vertex_entries_[i].type = Vertex_entry::INPUT;
    }

    for(Size i = 0; i < num_vertex_ranks(); ++i)
    {
      vertex_ranks_[i].reset();
    }
    vertex_rank_left_ = 0;
    vertex_rank_right_ = 0;

    // remove split edge entries
    for(Size i = 0; i < num_split_vertices_; ++i)
    {
      edge_entries_[num_input_edges_ + i].reset();
    }

    // remove split vertex entries
    for(Size i = 0; i < num_split_vertices_; ++i)
    {
      vertex_entries_[num_input_vertices_ + i] = Vertex_entry();
    }
    num_split_vertices_ = 0;
    num_split_vertex_entries_ = 0;

    num_skipped_edges_ = 0;

    // reset sector entries
    for(Size i = 0; i < sector_entries_.size(); ++i)
    {
      sector_entries_[i] = Sector_entry();
      sector_entries_[i].set.last = i;
    }

    query_point_pos_ = UNDEFINED;
  }

#if ENABLE_VP_DEBUG
  /// \name drawing methods
  //@{

  void draw_edge_ranks(CairoImage& img) const
  {
    CGAL_precondition(num_edge_ranks() > 0);

    CairoColor label_color(CGAL::BLACK);

    for(Size rank = 0; rank < num_edge_ranks(); ++rank)
    {
      const Edge_entry& edge_entry = ranked_edge(rank);

      if(!edge_entry.source.assigned() || !edge_entry.target.assigned())
      {
        continue;
      }

      CGAL::Point_2<Kernel> source = ranked_vertex(edge_entry.source).point,
          target = ranked_vertex(edge_entry.target).point;

      double rank_ratio = 1;

      if(num_edge_ranks() > 1);
      {
        CGAL_precondition(rank < num_edge_ranks());
        rank_ratio = num_edge_ranks() - rank - 1;
        rank_ratio = rank_ratio/(num_edge_ranks() - 1);
      }

      CairoColor color(0, rank_ratio, 0);

      img.draw_segment(source, target, color);
      img.label(CGAL::midpoint(source, target), tostr(rank), label_color);
    }
  }

  void draw_input(CairoImage& img) const
  {
    img.set_color(CGAL::BLACK);
    for(Size i = 0; i < num_input_edges(); ++i)
    {
      const Edge_entry& edge_entry = edge(i);

      CGAL_precondition(edge_entry.source.assigned()
                        && edge_entry.target.assigned());

      img.draw_segment(vertex(edge_entry.source).point,
                       vertex(edge_entry.target).point);
    }

    for(Size i = 0; i < num_input_vertices(); ++i)
    {
      img.draw_point(vertex(i).point);
    }

    img.draw_point(query_point(), CGAL::GREEN);
    draw_labels(img);
  }

  void draw_labels(CairoImage& img) const
  {
    for(Size index = 0; index < num_input_vertices(); ++index)
    {
      img.label(vertex(index).point, vertex(index).label, CGAL::BLACK);
    }
  }

  void draw_sectors(CairoImage& img) const
  {
    img.set_color(CGAL::RED);

    for(Size i = 0; i < num_sectors(); ++i)
    {
      if(!sector(i).visible_edge.assigned())
      {
        continue;
      }

      const Edge_entry& vis_edge = ranked_edge(sector(i).visible_edge);

      CGAL_precondition(vis_edge.target.assigned() && vis_edge.source.assigned());
      img.draw_segment(ranked_vertex(vis_edge.source).point, ranked_vertex(vis_edge.target).point);
    }
  }

  void draw_vertex_ranks(CairoImage& img, const CGAL::Point_2<Kernel>& query_point) const
  {
    CGAL::Aff_transformation_2<Kernel> next_line(CGAL::TRANSLATION, CGAL::Vector_2<Kernel>(0, -img.font_size()));

    for(Size rank = 0; rank < num_vertex_ranks(); ++rank)
    {
      const Vertex_entry& vertex_entry = ranked_vertex(rank);

      bool right_side = CGAL::less_x(query_point, vertex_entry.point);

      CGAL::Point_2<Kernel> label_pos(vertex_entry.point);

      // second label for split vertex
      if(vertex_entry.is_split_vertex() && is_second_entry(rank))
      {
        label_pos = next_line(label_pos);
        right_side = true;
      }

      double rank_ratio = 1;

      if(right_side && (num_vertex_ranks() > vertex_rank_left_ + 1))
      {
        CGAL_precondition(rank >= vertex_rank_left_);
        rank_ratio = (rank - vertex_rank_left_);
        rank_ratio = rank_ratio / (num_vertex_ranks() - vertex_rank_left_ - 1);
      }
      else if(!right_side && (vertex_rank_left_ > 1))
      {
        rank_ratio = rank;
        rank_ratio = rank_ratio / (vertex_rank_left_ - 1);
      }

      CGAL_precondition((0 <= rank_ratio) && (rank_ratio <= 1));
      CairoColor color(rank_ratio, 0, right_side?1:0);

      img.draw_point(vertex_entry.point, color);
      img.label(label_pos, tostr(rank), color);
    }
  }

  //@}
#endif

  /// \name set operations
  //@{
  /// return the last entry of the set which contains the entry with given index
  const Sector_index& find_last(const Sector_index& index)
  {
    CGAL_precondition(index.assigned());
    CGAL_precondition(in_range(index, 0, num_sectors() + 1));

    debug("FIND LAST of " + sector_label(index));
    return set(find_set_index(index)).last;
  }

  /// find the root of a set tree (i.e. the set the entry with given index belongs to)
  const Sector_index& find_set_index(const Sector_index& index) const
  {
    CGAL_precondition(index.assigned());
    CGAL_precondition(in_range(index, 0, num_sectors() + 1));

    // set is root of the tree
    if(!set(index).parent.assigned())
    {
      return index;
    }
    // find the root
    else
    {
      // avoid endless loop
      CGAL_precondition(set(index).parent != index);

      typename Set::Index_type root = set(index).parent;

      // path compression
      while(set(root).parent.assigned())
      {
        root = set(root).parent;
      }

      set(index).parent = root;

      return set(index).parent;
    }
  }

  /**
   * join set with given last entry and the following set
   * \return false if already at the last entry
   */
  void link_set(const Sector_index& last)
  {
    CGAL_precondition(last.assigned());
    CGAL_precondition(in_range(last, 0, num_sectors() + 1));

    debug("LINK SET " + sector_label(last) + " (" + tostr(last) +  ")");
    debug("  with " + sector_label(last+1) + " (" + tostr(last+1) + ")");

    Size first = find_set_index(last);
    CGAL_precondition(set(first).last == last);

    Size second = find_set_index(last+1);
    CGAL_precondition(set(first).last < set(second).last);

    // height of the set tree increses by one
    if(set(first).rank == set(second).rank)
    {
      set(first).rank++;
    }

    // attach the smaller set tree to the larger one
    if(set(first).rank > set(second).rank)
    {
      set(second).parent = first;
      set(first).last = set(second).last;

      CGAL_postcondition(set(second).parent.assigned());
      CGAL_postcondition(set(second).parent == first);
    }
    else
    {
      set(first).parent = second;

      CGAL_postcondition(set(first).parent.assigned());
      CGAL_postcondition(set(first).parent == second);
    }

    CGAL_postcondition(find_set_index(first) == find_set_index(second));
  }
  //@}

  /// \name member counts
  //@{

  /// \return number of total edges ranks (including split edges)
  Size num_edge_ranks() const
  {
    CGAL_precondition(num_split_vertices_ <= num_input_edges());
    CGAL_precondition(num_skipped_edges_ <= num_input_edges());

    return num_input_edges() + num_split_vertices_ - num_skipped_edges_;
  }

  /// \return number of input edges
  Size num_input_edges() const
  {
    CGAL_precondition(num_input_edges_ <= edge_entries_.size());

    return num_input_edges_;
  }

  /// \return number of input vertices
  Size num_input_vertices() const
  {
    CGAL_precondition(num_input_vertices_ < vertex_entries_.size());

    return num_input_vertices_;
  }

  /// \return number of visibility sectors
  Size num_sectors() const
  {
    CGAL_precondition(sector_entries_.size() > 1);

    return sector_entries_.size() - 1 - (is_query_point_on_vertex()?1:0);
  }

  /// \return number of split vertices
  Size num_split_vertices() const
  {
    CGAL_precondition(num_split_vertices_ <= num_input_edges());

    return num_split_vertices_;
  }

  /// \return total number of vertices to be ranked
  Size num_vertex_ranks() const
  {
    CGAL_precondition(num_input_vertices_ + num_split_vertex_entries_ <= vertex_ranks_.size());

    return num_input_vertices_ + num_split_vertex_entries_  - (is_query_point_on_vertex()?1:0);
  }

  //@}


  /// \name ranking
  //@{

  /**
   * rank the edge given by the indices of its vertices
   */
  void rank_edge(Vertex_index source, Vertex_index target)
  {
    if(!source.assigned() || !target.assigned())
    {
      CGAL_precondition(is_query_point_on_edge());
      CGAL_precondition(source.assigned() || target.assigned());

      return;
    }

    debug("RANK EDGE " + vertex(source).label + "->" + vertex(target).label);

    CGAL_precondition(edge_rank_ < num_edge_ranks());
    CGAL_precondition(vertex(source).is_input_vertex()
                      || vertex(target).is_input_vertex());

    // swap source to be split vertex
    if(vertex(target).is_split_vertex())
    {
      std::swap(source, target);
    }

    Vertex_index source_rank = vertex(source).rank,
        target_rank = vertex(target).rank;

    // found edge incident to query point
    if(!source_rank.assigned() || !target_rank.assigned())
    {
      CGAL_precondition(is_query_point_on_vertex());

      skip_edge();
      return;
    }

    // check for split vertex
    if(vertex(source).is_split_vertex())
    {
      CGAL_precondition(source_rank < vertex_shift());
      CGAL_precondition(vertex(target).is_input_vertex());

      // edge target is on right side
      if(target_rank >= vertex_rank_left_)
      {
        // use second split vertex entry
        source_rank = second_entry(source_rank);
        CGAL_postcondition(source_rank < num_vertex_ranks());
      }
    }

    // use smaller rank as source
    if(source_rank > target_rank)
    {
      std::swap(source_rank, target_rank);
    }

    edge_ranks_[edge_rank_].source = source_rank;
    edge_ranks_[edge_rank_].target = target_rank;

    ++edge_rank_;
  }

  void rank_vertex(const Size& index)
  {
    debug("RANK VERTEX " + vertex(index).label);

    CGAL_precondition(!vertex(index).is_split_vertex());
    CGAL_precondition(!is_query_point_on_vertex()
                      || (vertex(index).point != query_point()));

    bool right_of_query_point = (CGAL::compare_lexicographically(vertex(index).point, query_point_) == CGAL::LARGER);

    if(vertex_rank_right_ < vertex_rank_left_)
    {
      CGAL_precondition((vertex_rank_right_ == 0) && !right_of_query_point);
    }
    else if(vertex_rank_right_ == vertex_rank_left_)
    {
      CGAL_precondition(right_of_query_point);
    }

    Size& rank = right_of_query_point?vertex_rank_right_:vertex_rank_left_;

    CGAL_precondition(rank < num_vertex_ranks());
    CGAL_precondition(!vertex_ranks_[rank].assigned());
    vertex_ranks_[rank] = index;
    vertex(index).rank = rank;

    ++rank;
  }
  //@}

  /// \name insertion methods
  //@{

  /**
   * Creates vertex entry for given point. Inserted points are assumed to be
   * ordered by y-coordinate
   */
  Size insert_input_vertex(const CGAL::Point_2<Kernel>& point)
  {
    CGAL_precondition((0 <= num_input_vertices_)
                      && (num_input_vertices_ < vertex_entries_.size()));
    CGAL_precondition((0 == num_input_vertices())
                      || !CGAL::less_y(point, vertex(num_input_vertices_ - 1).point));

    vertex_entries_[num_input_vertices_].point = point;
    vertex_entries_[num_input_vertices_].type = Vertex_entry::INPUT;
#if ENABLE_VP_DEBUG
    vertex_entries_[num_input_vertices_].label = next_vertex_label();
#endif

    vertex_map_.insert(std::make_pair(point, num_input_vertices_));

    return num_input_vertices_++;
  }

  /// create edge entry for given end points
  void insert_input_edge(const CGAL::Point_2<Kernel>& source,
                         const CGAL::Point_2<Kernel>& target)
  {
    CGAL_precondition((0 <= num_input_edges_)
                      && (num_input_edges_ < edge_entries_.size()));

    edge_entries_[num_input_edges_].source = find_vertex(source);
    edge_entries_[num_input_edges_].target = find_vertex(target);

    ++num_input_edges_;
  }

  //@}

  /// \name member access
  //@{
  bool allow_collinear() const
  {
    return allow_collinear_;
  }

  /// \return vertex index for given point
  Vertex_index find_vertex(const CGAL::Point_2<Kernel>& point) const
  {
    typename Vertex_map::const_iterator it = vertex_map_.find(point);

    // not found
    if(it == vertex_map_.end())
    {
      return Vertex_index();
    }
    else
    {
      return it->second;
    }
  }

  /**
   * \param index sector to check
   * \param start start of range (inclusively)
   * \param end end of range (exclusively)
   * \return true if sector index is in the given range
   */
  bool in_range(const Sector_index& index, const Sector_index& start, const Sector_index& end) const
  {
    CGAL_precondition(start.assigned() == end.assigned());

    if(!start.assigned())
    {
      return false;
    }

    if(start <= end)
    {
      return (start <= index) && (index < end);
    }
    else
    {
      return (index < end) || (start <= index);
    }
  }
  //@}

  /// \name query point access
  //@{
  const CGAL::Point_2<Kernel>& query_point() const
  {
    CGAL_precondition(is_query_point_assigned());

    return query_point_;
  }

  void insert_query_point(const CGAL::Point_2<Kernel>& query_point)
  {
    CGAL_precondition(!is_query_point_assigned());

    query_point_ = query_point;

    find_fake_split_vertices();
    find_split_edges();

    // neither on edge nor on vertex
    if(!is_query_point_on_edge()
       && !is_query_point_on_vertex())
    {
      query_point_pos_ = IN_FACE;
    }

    CGAL_postcondition(is_query_point_assigned());

    // rank split vertices
    for(Size i = 0; i < num_input_vertices() + num_split_vertices(); ++i)
    {
      if(vertex(i).is_split_vertex())
      {
        rank_split_vertex(i);
      }
    }
  }

  bool is_query_point_assigned() const
  {
    return (query_point_pos_ != UNDEFINED);
  }

  bool is_query_point_on_edge() const
  {
    return (query_point_pos_ == ON_EDGE);
  }

  bool is_query_point_on_vertex() const
  {
    return (query_point_pos_ == ON_VERTEX);
  }
  //@}

  /// \name entry access
  //@{

  /// \return reference to edge entry with given index
  Edge_entry& edge(const Size& index)
  {
    CGAL_precondition((0 <= index)
                      && (index < num_input_edges_ + num_split_vertices_));
    CGAL_precondition(index < edge_entries_.size());

    return edge_entries_[index];
  }

  /// \return const reference to edge entry with given index
  const Edge_entry& edge(const Size& index) const
  {
    CGAL_precondition((0 <= index)
                      && (index < num_input_edges_ + num_split_vertices_));
    CGAL_precondition(index < edge_entries_.size());

    return edge_entries_[index];
  }

  /// \return const reference to edge entry with given rank
  const Edge_entry& ranked_edge(const Size& rank) const
  {
    CGAL_precondition((0 <= rank) && (rank < num_edge_ranks()));

    return edge_ranks_[rank];
  }

  /// \return const reference to vertex entry with given rank
  const Vertex_entry& ranked_vertex(const Size& rank) const
  {
    CGAL_precondition((0 <= rank) && (rank < num_vertex_ranks()));
    CGAL_precondition(vertex_ranks_[rank].assigned());

    return vertex(vertex_ranks_[rank]);
  }

  /// \return reference to sector entry with given index
  Sector_entry& sector(const Size& index)
  {
    CGAL_precondition(in_range(index, 0, num_sectors()));

    return sector_entries_[index];
  }

  /// \return const reference to sector entry with given index
  const Sector_entry& sector(const Size& index) const
  {
    CGAL_precondition(in_range(index, 0, num_sectors()));

    return sector_entries_[index];
  }

  /// \return const reference to split edge entry with given index
  const Edge_entry& split_edge(const Size& index) const
  {
    CGAL_precondition((0 <= index) && (index < num_split_vertices_));

    return edge(num_input_edges_ + index);
  }

  /// \return const reference to split vertex entry with given index
  const Vertex_entry& split_vertex(const Size& index) const
  {
    CGAL_precondition((0 <= index) && (index < num_split_vertices_));

    return vertex(num_input_vertices_ + index);
  }

  /// \return reference to vertex entry with given index
  Vertex_entry& vertex(const Size& index)
  {
    CGAL_precondition((0 <= index)
                      && (index < num_input_vertices_ + num_split_vertices_));
    CGAL_precondition(index < vertex_entries_.size());

    return vertex_entries_[index];
  }

  /// \return const reference to vertex entry with given index
  const Vertex_entry& vertex(const Size& index) const
  {
    CGAL_precondition((0 <= index)
                      && (index < num_input_vertices_ + num_split_vertices_));
    CGAL_precondition(index < vertex_entries_.size());

    return vertex_entries_[index];
  }
  //@}

#if ENABLE_VP_DEBUG
  /// \name labels
  //@{

  /// \return label for sector with given index
  std::string sector_label(const Size& index) const
  {
    CGAL_precondition(in_range(index, 0, num_sectors() + 1));

    Size from, to;

    if(index == 0)
    {
      from = 0;
      to = vertex_shift();
    }
    else if(index == num_sectors() - 1)
    {
      CGAL_precondition(num_input_vertices() > 1);

      from = num_input_vertices() + num_split_vertices() - 2;
      to = (num_split_vertex_entries_ == 0)? 0:num_vertex_ranks() - 1;
    }
    else if(index == num_sectors())
    {
      return "#";
    }
    else
    {
      from = index + vertex_shift() - 1;
      to = index + vertex_shift();
    }

    return ranked_vertex(from).label + "-" + ranked_vertex(to).label;
  }
  //@}
#endif

  /// \name index calculations
  //@{

  /// map vertex rank to sector index
  Size sector_after(const Vertex_index& vertex_rank) const
  {
    CGAL_precondition(vertex_rank.assigned());

    if(ranked_vertex(vertex_rank).is_split_vertex())
    {
      // first entry of split vertex
      if(vertex_rank < vertex_shift())
      {
        return 0;
      }
      // second entry of split vertex
      else
      {
        return num_sectors();
      }
    }
    else
    {
      // map first input vertex to 1
      return vertex_rank - vertex_shift() + 1;
    }
  }

  /// map sector index to vertex rank
  Size vertex_before(const Sector_index& sector_index) const
  {
    CGAL_precondition(sector_index.assigned());

    if(sector_index == 0)
    {
      return 0;
    }
    else if(sector_index == num_sectors())
    {
      return num_vertex_ranks() - 1;
    }
    else
    {
      CGAL_precondition(sector_index + vertex_shift() - 1 < num_vertex_ranks());
      return sector_index + vertex_shift() - 1;
    }
  }

  bool is_second_entry(const Size& rank) const
  {
    CGAL_precondition(ranked_vertex(rank).is_split_vertex());

    return (rank >= second_entries_start());
  }

  bool same_vertex(const Vertex_index& rank1, const Vertex_index& rank2) const
  {
    if(!rank1.assigned() || !rank2.assigned())
    {
      return false;
    }

    if((ranked_vertex(rank1).is_split_vertex()?second_entry(rank1):rank1)
       == (ranked_vertex(rank2).is_split_vertex()?second_entry(rank2):rank2))
    {
      return true;
    }

    return false;
  }

  Vertex_index second_entry(const Vertex_index& rank) const
  {
    CGAL_precondition(ranked_vertex(rank).is_split_vertex());

    if(is_second_entry(rank))
    {
      return rank;
    }

    CGAL_precondition(rank < vertex_shift());
    CGAL_precondition(second_entries_start() + rank < num_vertex_ranks());

    return second_entries_start() + rank;
  }

  //@}

#if ENABLE_VP_DEBUG
  /// prints all input edges
  void dump_input_edges()
  {
    debug("Edges:");

    for(Size i = 0; i < num_input_edges(); ++i)
    {
      debug("  " + tostr(i) + ": " + tostr(edge(i).source) + " -> "
            + tostr(edge(i).target));
    }
  }

  /// prints all ranked edges
  void dump_ranked_edges()
  {
    debug("Ranked edges:");

    for(Size i = 0; i < num_edge_ranks(); ++i)
    {
      debug("  " + tostr(i) + ": " + tostr(ranked_edge(i).source) + " -> "
            + tostr(ranked_edge(i).target));
    }
  }

  /// prints all sectors along with their visibile edge
  void dump_sectors()
  {
    debug("Sectors:");

    for(Size i = 0; i < num_sectors(); ++i)
    {
      std::string edge_label;

      if(sector(i).visible_edge.assigned())
      {
        const Edge_entry& edge_entry = ranked_edge(sector(i).visible_edge);
        edge_label = ranked_vertex(edge_entry.source).label + "->"
            + ranked_vertex(edge_entry.target).label;
      }

      debug("  " + sector_label(i) + ": " + edge_label);
    }
  }

  /// prints all elements of the set with given index
  void dump_set(const Sector_index& index) const
  {
    CGAL_precondition(index.assigned());

    Size set_index = find_set_index(index), left_most = set_index;

    // find the first entry that belongs to the set
    while((left_most > 0)
          && set(left_most-1).parent.assigned()
          && (set(left_most-1).parent == set_index))
    {
      --left_most;
    }

    while(left_most <= set(set_index).last)
    {
      std::cout << sector_label(left_most) << ", ";
      ++left_most;
    }

    std::cout << std::endl;
  }

  /// prints all input vertices
  void dump_input_vertices()
  {
    debug("Vertices:");

    for(Size i = 0; i < num_input_vertices(); ++i)
    {
      debug("  " + tostr(i) + ": " + vertex(i).label);
    }
  }

#endif

private:
  bool allow_collinear_;
  /// \name containers
  //@{
  /// edge entries for input edges
  Edge_entries edge_entries_;
  /// edges (including split edges) ordered by edge rank
  Edge_ranks edge_ranks_;

  /// entries for visibility sectors
  Sector_entries sector_entries_;

  /// associated data of vertices (including split vertices)
  Vertex_entries  vertex_entries_;
  /// map from point to vertex index (only input vertices)
  Vertex_map      vertex_map_;
  /// vertex indices ordered by rank (including two entries per split vertex)
  Vertex_ranks    vertex_ranks_;
  //@}

  /// \name counts
  //@{
  /// current number of edges ignored by the ranking
  Size num_skipped_edges_;
  /// current number of split vertices/edges
  Size num_split_vertices_;
  /// current number of additional ranks for split vertices
  Size num_split_vertex_entries_;
  /// current number of input edges
  Size num_input_edges_;
  /// current number of input vertices
  Size num_input_vertices_;
  //@}

  /// \name current ranks
  //@{
  /// current edge rank
  Size edge_rank_;
  /// current vertex rank left of query point
  Size vertex_rank_left_;
  /// current vertex rank right of query point
  Size vertex_rank_right_;
  //@}

  /// query point
  CGAL::Point_2<Kernel> query_point_;
  /// position of query point \see Query_point_pos
  Query_point_pos query_point_pos_;

  /// \name set entry access
  //@{
  const Set& set(const Size& index) const
  {
    CGAL_precondition(in_range(index, 0, num_sectors() + 1));

    return sector_entries_[index].set;
  }

  Set& set(const Size& index)
  {
    CGAL_precondition(in_range(index, 0, num_sectors() + 1));

    return sector_entries_[index].set;
  }
  //@}

  /**
   * \return the difference between current input vertex rank and the rank of
   * same vertex in absence of split vertices
   */
  Size vertex_shift() const
  {
    if((num_split_vertex_entries_ == 0) && (num_split_vertices_ == 0))
    {
      return 0;
    }

    CGAL_precondition(num_split_vertex_entries_ > num_split_vertices_);

    return num_split_vertex_entries_ - num_split_vertices_;
  }

#if ENABLE_VP_DEBUG
  /// \name labels
  //@{

  std::string edge_label(const Edge_entry& entry) const
  {
    CGAL_precondition(entry.source.assigned() && entry.target.assigned());

    return vertex(entry.source).label + "->" + vertex(entry.target).label;
  }

  std::string next_vertex_label() const
  {
    std::string label;
    Size index = num_input_vertices_;

    do
    {
      label.insert(label.begin(), 'A' + (index % 26));
      index = index/26;
    } while(index > 0);

    return label;
  }

  //@}
#endif

  /// increase the number of ranks to the left of query point
  void add_left_rank()
  {
    CGAL_precondition((0 <= vertex_rank_right_)
                      && (vertex_rank_right_ < num_input_vertices() + num_split_vertices()));

    ++vertex_rank_right_;
  }

  /**
   * finds all input vertices on upward ray of query point
   */
  void find_fake_split_vertices()
  {
    CGAL_precondition(num_split_vertex_entries_ == 0);

    for(Size i = 0; i < num_input_vertices(); ++i)
    {
      CGAL::Comparison_result comp_x = CGAL::compare_x(vertex(i).point, query_point_);

      if(comp_x == CGAL::EQUAL)
      {
        CGAL::Comparison_result comp_y = CGAL::compare_y(vertex(i).point, query_point_);

        // found vertex above query point
        if(comp_y == CGAL::LARGER)
        {
          CGAL_precondition(vertex(i).is_input_vertex());

          vertex(i).type = Vertex_entry::FAKE_SPLIT;

          CGAL_precondition(num_vertex_ranks() < vertex_ranks_.size());
          num_split_vertex_entries_ += 1;
        }
        // found vertex below query point
        else if(comp_y ==CGAL::SMALLER)
        {
          add_left_rank();
        }
        // found query point
        else
        {
          CGAL_assertion(comp_y == CGAL::EQUAL);

          CGAL_precondition(!is_query_point_assigned());

          query_point_pos_ = ON_VERTEX;
        }
      }
      // found vertex left of query point
      else if(comp_x == CGAL::SMALLER)
      {
        add_left_rank();
      }
    }
  }

  /**
   * finds all edges intersecting the upward ray from query point and inserts
   * split vertices
   */
  void find_split_edges()
  {
    CGAL_precondition(num_split_vertices_ == 0);

    for(Size i = 0; i < num_input_edges(); ++i)
    {
      const Edge_entry& edge_entry = edge(i);

      CGAL_precondition(edge_entry.source.assigned()
                        && edge_entry.target.assigned());
      CGAL::Segment_2<Kernel> segment(vertex(edge_entry.source).point,
                                      vertex(edge_entry.target).point);

      if(segment.has_on(query_point_))
      {
        if((segment.source() != query_point_)
           && (segment.target() != query_point_))
        {
          CGAL_precondition(!is_query_point_assigned());

          query_point_pos_ = ON_EDGE;
          skip_edge();
        }

        continue;
      }

      // don't split vertical edges
      if(segment.is_vertical())
      {
        continue;
      }

      // check if edge vertices lie on opposite sides of vertical ray
      if((CGAL::compare_x(segment.min(), query_point_) == CGAL::SMALLER)
         && (CGAL::compare_x(segment.max(), query_point_) == CGAL::LARGER))
      {
        // calculate intersection point with ray
        CGAL::Point_2<Kernel> intersection(query_point_.x(), segment.supporting_line().y_at_x(query_point_.x()));

        if(intersection.y() > query_point_.y())
        {
          debug("  found split edge " + edge_label(edge_entry));

          CGAL_precondition((intersection != segment.min()) && (intersection != segment.max()));

          insert_split_vertex(intersection, i);
        }
      }
    }
  }

  /// inserts entry for split vertex with given coordinates and split edge
  void insert_split_vertex(const CGAL::Point_2<Kernel>& point,
                           const Size& split_edge)
  {
    CGAL_precondition((0 <= num_split_vertices_)
                      && (num_split_vertices_ < num_input_edges_));
    CGAL_precondition(num_input_vertices_ + num_split_vertices_ < vertex_entries_.size());

    vertex_entries_[num_input_vertices_ + num_split_vertices_].point = point;
    vertex_entries_[num_input_vertices_ + num_split_vertices_].type = Vertex_entry::SPLIT;

#if ENABLE_VP_DEBUG
    CGAL_precondition(edge(split_edge).source.assigned()
                      && edge(split_edge).target.assigned());

    std::string label = vertex(edge(split_edge).source).label
        + vertex(edge(split_edge).target).label;
    std::transform(label.begin(), label.end(), label.begin(), tolower);

    vertex_entries_[num_input_vertices_ + num_split_vertices_].label = label;
#endif

    CGAL_precondition(num_input_edges_ + num_split_vertices_ < edge_entries_.size());

    edge_entries_[num_input_edges_ + num_split_vertices_] = edge(split_edge);

    CGAL_precondition(num_vertex_ranks() + 1 < vertex_ranks_.size());
    num_split_vertex_entries_ += 2;

    CGAL_precondition(num_split_vertices_ < num_input_edges());
    num_split_vertices_ += 1;
  }

  /**
   * ranks a split vertex with given index
   */
  void rank_split_vertex(const Size& index)
  {
    debug("RANK SPLIT VERTEX " + vertex(index).label);

    add_left_rank();

    CGAL_precondition(vertex_rank_left_ < vertex_rank_right_);
    CGAL_precondition(vertex_rank_left_ < num_vertex_ranks());

    Size rank = vertex_rank_left_;

    vertex(index).rank = rank;

    // first entry
    CGAL_precondition(rank < num_vertex_ranks());
    CGAL_precondition(!vertex_ranks_[rank].assigned());
    vertex_ranks_[rank] = index;

    rank = second_entry(rank);

    // second entry
    CGAL_precondition(rank < num_vertex_ranks());
    CGAL_precondition(!vertex_ranks_[rank].assigned());
    vertex_ranks_[rank] = index;

    ++vertex_rank_left_;
  }

  /// increase counter of edges which will not be ranked
  void skip_edge()
  {
    debug("  skip edge");

    CGAL_precondition(num_skipped_edges_ < num_input_edges());

    ++num_skipped_edges_;
  }

  /// \return second rank of first split vertex
  Size second_entries_start() const
  {
    CGAL_precondition(num_input_vertices() >= (is_query_point_on_vertex()?1:0));

    return num_input_vertices() + num_split_vertices()
        - (is_query_point_on_vertex()?1:0);
  }
};

#endif // VISIBILITY_STORAGE_H
