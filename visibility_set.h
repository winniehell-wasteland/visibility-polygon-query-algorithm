#ifndef VISIBILITY_SET_H
#define VISIBILITY_SET_H

#include "visibility_index.h"

template <typename IndexType>
class Visibility_set
{
public:
  /// \name template parameters
  //@{
  typedef IndexType Index_type;
  //@}

  /// default constructor
  Visibility_set() :
    last(),
    parent(),
    rank(0)
  {

  }

  /// \return true if no data has been assigned to this set
  bool empty() const
  {
    return !parent.assigned()
        && (rank == 0)
        && !last.assigned();
  }

  /// last set in subtrees (i.e. the one with highest index)
  Visibility_index<Index_type> last;
  // TODO: can we use pointer for parent?
  /// root of the set tree or node on the path to the root
  mutable Visibility_index<Index_type> parent;
  /// approximate height of the subtrees
  Index_type rank;
};

#endif // VISIBILITY_SET_H
