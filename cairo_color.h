#ifndef CAIRO_COLOR_H
#define CAIRO_COLOR_H

#include <cairo.h>

#include <CGAL/IO/Color.h>

class CairoColor
{
public:
	/// default constructor
	CairoColor() :
		pattern_(0)
	{
	}

	CairoColor(const CGAL::Color& other) :
		pattern_(0)
	{
		*this = other;
	}

	CairoColor(double red,
						 double green,
						 double blue,
						 double alpha = 1.0) :
		pattern_(cairo_pattern_create_rgba(red, green, blue, alpha))
	{
	}

	~CairoColor()
	{
		if(pattern_ != 0)
		{
			cairo_pattern_destroy(pattern_);
		}
	}

	CairoColor& operator=(const CairoColor& other)
	{
		if(pattern_ != 0)
		{
			cairo_pattern_destroy(pattern_);
		}

		pattern_ = cairo_pattern_reference(other.pattern_);

		return *this;
	}

	CairoColor& operator=(const CGAL::Color& other)
	{
		if(pattern_ != 0)
		{
			cairo_pattern_destroy(pattern_);
		}

		pattern_ = cairo_pattern_create_rgba(static_cast<double>(other.red())/255, static_cast<double>(other.green())/255, static_cast<double>(other.blue()/255), static_cast<double>(other.alpha())/120);

		return *this;
	}

	operator cairo_pattern_t* () const
	{
		return pattern_;
	}

private:
	cairo_pattern_t* pattern_;
};

#endif // CAIRO_COLOR_H
