#ifndef ER_STD_SORT_H
#define ER_STD_SORT_H

template <class Kernel_,
          class Storage_,
          class InputContainer>
class Er_Std_Sort
{
public:
  /// \name template parameters
  //@{
  typedef InputContainer  Input_container;
  typedef Kernel_         Kernel;
  typedef Storage_        Storage;
  //@}

  /// default constructor
  Er_Std_Sort(Storage& storage) :
    storage_(storage)
  {
    // TODO:
    assert(0);
  }

private:
  Storage& storage_;
};

#endif // ER_STD_SORT_H
