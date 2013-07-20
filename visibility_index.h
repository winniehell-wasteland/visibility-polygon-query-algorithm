#ifndef VISIBILITY_INDEX_H
#define VISIBILITY_INDEX_H

template <typename IndexType>
class Visibility_index
{
public:
	/// \name template parameters
	//@{
	typedef IndexType Index_type;
	//@}

	/// default constructor
	Visibility_index() :
		assigned_(false),
		value_()
	{

	}

	/// copy constructor
	Visibility_index(const Visibility_index& other) :
		assigned_(other.assigned_),
		value_(other.value_)
	{

	}

	/// wrapping constructor
	Visibility_index(const Index_type& other) :
		assigned_(true),
		value_(other)
	{

	}

	/// assignment operator
	const Visibility_index& operator =(const Visibility_index& other)
	{
		assigned_ = other.assigned();

		if(other.assigned())
		{
			value_ = other.value_;
		}

		return *this;
	}

	/// wrapping assignment operator
	const Visibility_index& operator =(const Index_type& other)
	{
		assigned_ = true;
		value_ = other;

		return *this;
	}

	/// unwrapping operator
	operator Index_type() const
	{
		assert(assigned_);
		return value_;
	}

	/// \return true if a value has been assigned
	bool assigned() const
	{
		return assigned_;
	}

	/// virtually unassigns value
	void reset()
	{
		assigned_ = false;
	}

private:
	/// true if a value has been assigned
	bool assigned_;
	/// index value
	Index_type value_;
};

#endif // VISIBILITY_INDEX_H
