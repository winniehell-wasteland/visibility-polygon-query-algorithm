#ifndef CAIRO_IMAGE_H
#define CAIRO_IMAGE_H

#include "cgal_defs.h"

#include <cairo.h>
#include <cairo-pdf.h>
#include <cairo-svg.h>

#include <CGAL/box_intersection_d.h>

#include "cairo_color.h"

class CairoImage
{
public:
  static const Kernel::RT EPS_FRACTION()
  {
    return Kernel::RT("1/1000000000000000");
  }

  enum file_format
  {
    PDF,
    SVG,
    FILE_FORMAT = SVG
  };

  cairo_surface_t* cairo_surface;
  cairo_t* cairo_ref;
  double xscale_, yscale_;

  /// \name sizes
  //@{
  double font_size() const
  {
    return 2*point_width();
  }

  double line_width() const
  {
    return 0.25;
  }

  double padding() const
  {
    return 4*point_width();
  }

  double point_width() const
  {
    return 2*line_width();
  }
  //@}

  /**
   * initializes cairo handles
   */
  CairoImage(const std::string& filename, const CGAL::Bbox_2& bbox, double xscale = 1.0, double yscale = 1.0) :
    xscale_(xscale), yscale_(yscale)
  {
    double width = xscale*(bbox.xmax() - bbox.xmin()) + 2*padding(),
        height = yscale*(bbox.ymax() - bbox.ymin()) + 2*padding();

    switch(FILE_FORMAT)
    {
      case PDF:
      {
        cairo_surface = static_cast<cairo_surface_t*>( cairo_pdf_surface_create((filename + "." + extension()).c_str(), width, height) );
        break;
      }
      case SVG:
      {
        cairo_surface = static_cast<cairo_surface_t*>( cairo_svg_surface_create((filename + "." + extension()).c_str(), width, height) );
        break;
      }
      default: assert(0);
    }
    cairo_ref = cairo_create(cairo_surface);
    cairo_set_line_width(cairo_ref, line_width());

    cairo_matrix_t matrix;

    cairo_matrix_init(&matrix, 1, 0, 0, -1, -xscale*bbox.xmin() + padding(), yscale*bbox.ymax() + padding());
    cairo_transform(cairo_ref, &matrix);
  }

  /**
   * destroys cairo handles
   */
  ~CairoImage()
  {
    cairo_destroy (cairo_ref);
    cairo_surface_destroy (cairo_surface);
  }

  /// \name draw points
  //@{
  void draw_point( const CGAL::Point_2<Kernel>&  point )
  {
    cairo_new_path(cairo_ref);
    cairo_arc (cairo_ref, xscale_*CGAL::to_double(point.x()), yscale_*CGAL::to_double(point.y()), point_width() / 2, 0, 2 * M_PI);

    cairo_fill(cairo_ref);
  }

  void draw_point( const CGAL::Point_2<Kernel>&  point, const CairoColor& color )
  {
    set_color(color);
    draw_point(point);
  }

  template <class Container>
  void draw_vertex( const typename Container::Vertex& vertex)
  {
    CGAL_precondition_msg(!vertex.is_at_open_boundary(), "Can not draw infinite vertex!");

    draw_point(vertex.point());
  }

  template <class Container>
  void draw_vertex( const typename Container::Vertex& vertex,  const CairoColor& color )
  {
    set_color(color);
    draw_vertex<Container>(vertex);
  }

  //@}

  /**
   * \name draw segments
   */
  //@{
  void draw_segment( const CGAL::Point_2<Kernel>& source, const CGAL::Point_2<Kernel>& target)
  {
    cairo_move_to(cairo_ref, xscale_*CGAL::to_double(source.x()), yscale_*CGAL::to_double(source.y()));
    cairo_line_to(cairo_ref, xscale_*CGAL::to_double(target.x()), yscale_*CGAL::to_double(target.y()));

    cairo_stroke( cairo_ref );
  }

  void draw_segment( const CGAL::Point_2<Kernel>& source, const CGAL::Point_2<Kernel>& target, const CairoColor& color )
  {
    set_color(color);
    draw_segment(source, target);
  }

  template <class Container>
  void draw_edge( const typename Container::Vertex& source, const typename Container::Vertex& target )
  {
    CGAL_precondition_msg(!source.is_at_open_boundary() && !target.is_at_open_boundary(), "Can not draw infinite vertex!");

    draw_segment(source.point(), target.point());
  }

  template <class Container>
  void draw_edge( const typename Container::Vertex& source, const typename Container::Vertex& target, const CairoColor& color )
  {
    set_color(color);
    draw_edge<Container>(source, target);
  }

  template <class Container>
  void draw_edge( const typename Container::Halfedge& edge )
  {
    CGAL_precondition_msg(!edge.is_fictitious(), "Can not draw infinite edge!");

    draw_edge<Container>(*edge.source(), *edge.target());
  }

  template <class Container>
  void draw_edge( const typename Container::Halfedge& edge, const CairoColor& color )
  {
    set_color(color);
    draw_edge<Container>(edge);
  }

  //@}

  /**
   * \name draw triangles
   */
  //@{
  void draw_triangle(const CGAL::Triangle_2<Kernel>& triangle)
  {
    cairo_new_path(cairo_ref);
    cairo_move_to(cairo_ref, xscale_*CGAL::to_double(triangle.vertex(0).x()), yscale_*CGAL::to_double(triangle.vertex(0).y()));
    cairo_line_to(cairo_ref, xscale_*CGAL::to_double(triangle.vertex(1).x()), yscale_*CGAL::to_double(triangle.vertex(1).y()));
    cairo_line_to(cairo_ref, xscale_*CGAL::to_double(triangle.vertex(2).x()), yscale_*CGAL::to_double(triangle.vertex(2).y()));
    cairo_line_to(cairo_ref, xscale_*CGAL::to_double(triangle.vertex(0).x()), yscale_*CGAL::to_double(triangle.vertex(0).y()));

    cairo_fill(cairo_ref);
  }

  void draw_triangle(const CGAL::Triangle_2<Kernel>& triangle, const CairoColor& color)
  {
    set_color(color);
    draw_triangle(triangle);
  }

  template <class Triangulation>
  void draw_triangle(const Triangulation& triangulation,
                     typename Triangulation::Face_const_handle face)
  {
    CGAL_precondition_msg(!triangulation.is_infinite(face), "Can not draw infinite face!");

    draw_triangle(triangulation.triangle(face));
  }

  template <class Triangulation>
  void draw_triangle(const Triangulation& triangulation,
                     typename Triangulation::Face_const_handle face,
                     const CairoColor& color)
  {
    set_color(color);
    draw_triangle<Triangulation>(triangulation, face);
  }

  //@}

  template <class Arrangement>
  void draw_arrangement(const Arrangement& arrangement, const CairoColor& vertex_color, const CairoColor& edge_color)
  {
    set_color(edge_color);
    for(typename Arrangement::Halfedge_const_iterator edge = arrangement.edges_begin(); edge != arrangement.edges_end(); ++edge)
    {
      draw_edge<Arrangement>(*edge);
    }

    set_color(vertex_color);
    for(typename Arrangement::Vertex_const_iterator vertex = arrangement.vertices_begin(); vertex != arrangement.vertices_end(); ++vertex)
    {
      draw_vertex<Arrangement>(*vertex);
    }
  }

  template <class Arrangement>
  void draw_arrangement(const Arrangement& arrangement)
  {
    draw_arrangement(arrangement, CairoColor(CGAL::BLACK), CairoColor(CGAL::BLACK));
  }

  void draw_dual_line(const CGAL::Line_2<Kernel>& line,
                      const CGAL::Bbox_2& bbox)
  {
    CGAL::Point_2<Kernel> p1(bbox.xmin(), line.y_at_x(bbox.xmin())),
        p2(bbox.xmax(), line.y_at_x(bbox.xmax()));

    draw_dual_line(line, p1, p2, bbox);
  }

  void draw_dual_line(const CGAL::Line_2<Kernel>& line,
                      const CGAL::Bbox_2& bbox,
                      const CairoColor& color)
  {
    set_color(color);
    draw_dual_line(line, bbox);
  }

  template <class DualPlane>
  void draw_dual_plane(const DualPlane& dual_plane,
                 const CGAL::Bbox_2& bbox,
                 const CairoColor& color)
  {
    set_color(color);

    for(typename DualPlane::Edge_const_iterator edge = dual_plane.edges_begin(); edge != dual_plane.edges_end(); ++edge)
    {
      if(edge->is_fictitious())
      {
        continue;
      }

      CGAL::Point_2<Kernel> from, to;

      if(edge->source()->is_at_open_boundary())
      {
        if(edge->direction() == CGAL::ARR_LEFT_TO_RIGHT)
        {
          from = CGAL::Point_2<Kernel>(bbox.xmin(), edge->curve().supporting_line().y_at_x(bbox.xmin()));
        }
        else
        {
          from = CGAL::Point_2<Kernel>(bbox.xmax(), edge->curve().supporting_line().y_at_x(bbox.xmax()));
        }
      }
      else
      {
        from = edge->source()->point();
        label(from, edge->source()->data().label);
      }

      if(edge->target()->is_at_open_boundary())
      {
        if(edge->direction() == CGAL::ARR_LEFT_TO_RIGHT)
        {
          to = CGAL::Point_2<Kernel>(bbox.xmax(), edge->curve().supporting_line().y_at_x(bbox.xmax()));
        }
        else
        {
          to = CGAL::Point_2<Kernel>(bbox.xmin(), edge->curve().supporting_line().y_at_x(bbox.xmin()));
        }
      }
      else
      {
        to = edge->target()->point();
        label(to, edge->target()->data().label);
      }

      draw_dual_line(edge->curve().supporting_line(), from, to, bbox);
    }
  }

  template <class DualPlane>
  void draw_dual_plane(const DualPlane& dual_plane,
                 const CGAL::Bbox_2& bbox)
  {
    draw_dual_plane(dual_plane, bbox, CairoColor(CGAL::BLACK));
  }

  inline const char* extension() const
  {
    switch(FILE_FORMAT)
    {
      case PDF: return "pdf";
      case SVG: return "svg";
      default: return 0;
    }
  }

  void label( const CGAL::Point_2<Kernel>& position, const std::string& label_text )
  {
    cairo_select_font_face(cairo_ref, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
    cairo_set_font_size(cairo_ref, font_size());

    cairo_move_to(cairo_ref, xscale_*CGAL::to_double(position.x()) + point_width()/2, yscale_*CGAL::to_double(position.y()) + point_width()/2);

    // do not scale text
    cairo_save(cairo_ref);
    cairo_identity_matrix(cairo_ref);
    cairo_show_text(cairo_ref, label_text.c_str());
    cairo_restore(cairo_ref);
  }

  void label( const CGAL::Point_2<Kernel>& position, const std::string& label_text, const CairoColor& color )
  {
    set_color(color);
    label(position, label_text);
  }

  void set_color(const CairoColor& color)
  {
    cairo_set_source( cairo_ref, color );
  }
private:
  /// draw dual line between two points
  void draw_dual_line(const Kernel::Line_2& line,
                      CGAL::Point_2<Kernel> from,
                      CGAL::Point_2<Kernel> to,
                      const CGAL::Box_intersection_d::Box_d<double, 2>& bbox)
  {
    if(to < from)
    {
      std::swap(from, to);
    }

    if(from.y() < bbox.min_coord(1))
    {
      from = CGAL::Point_2<Kernel>(line.x_at_y(bbox.min_coord(1)), bbox.min_coord(1));

      CGAL_precondition_msg(from.x() - typename Kernel::RT(bbox.max_coord(0)) <= EPS_FRACTION(), "First point is out of bounding box.");
    }

    if(to.y() > bbox.max_coord(1))
    {
      to = CGAL::Point_2<Kernel>(line.x_at_y(bbox.max_coord(1)), bbox.max_coord(1));

      CGAL_precondition_msg(to.x() - typename Kernel::RT(bbox.max_coord(0)) <= EPS_FRACTION(), "Second point is out of bounding box.");
    }

    draw_segment(from, to);
  }
};

#endif // CAIRO_IMAGE_H
