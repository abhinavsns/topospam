// tools.hpp
// some utility functions

#ifndef TOOLS_HPP_
#define TOOLS_HPP_

#include <cmath>
#include <memory>
#include <sstream>
#include <iostream>
#include <limits>
#include <vector>
#include <map>

/** Modulo function that works correctly with negative values */
template<class T>
inline T modu(const T& num, const T& div)
{
  if(num < 0) return div + num%div;
  else return num%div;
}

/** Modulo function that works correctly with negative values */
inline double modu(double num, double div)
{
  if(num < 0) return div + std::fmod(num, div);
  else return std::fmod(num, div);
}

/** Unsigned difference
 *
 * This function can be used with unsigned types with no fear of loosing the
 * sign.
 * */
template<class T>
inline T diff(const T& a, const T& b)
{
  return a>b ? a-b : b-a;
}

/** Check that two doubles are equal
 *
 * This somewhat complicated expression avoids false negative from finite
 * precision of the floating point arithmetic. Precision threshold can be
 * set using the error_factor parameter.
 * */
inline bool check_equal(double a, double b, double error_factor=1.)
{
  return a==b ||
    std::abs(a-b)<std::abs(std::min(a,b))*std::numeric_limits<double>::epsilon()*
                  error_factor;
}

/** Wrap point around periodic boundaries
 *
 * This function is useful when dealing with periodic boundary conditions,
 * simply returns the minimum of x%L and L-x%L
 * */
template<class T, class U>
inline T wrap(const T& x, const U& L)
{
  return std::min(x%L, L-x%L);
}

/** Specialization for double */
inline double wrap(double x, double L)
{
  const auto y = modu(x, L);
  if(abs(y)<abs(L-y))
    return y;
  else
    return L-y;
}

/** Set if smaller
 *
 * This function sets the first variable to be equal to the second if it is
 * smaller.
 * */
template<class T, class U>
inline void set_if_smaller(T& dst, const U& src)
{
  if(src<dst) dst = src;
}

/** Set if bigger
 *
 * This function sets the first variable to be equal to the second if it is
 * bigger.
 * */
template<class T, class U>
inline void set_if_bigger(T& dst, const U& src)
{
  if(src>dst) dst = src;
}

namespace detail
{
  /** Convert to strig and catenate arguments */
  template<class Head>
  void inline_str_add_args(std::ostream& stream, Head&& head)
  {
    stream << std::forward<Head>(head);
  }
  /** Convert to strig and catenate arguments */
  template<class Head, class... Tail>
  void inline_str_add_args(std::ostream& stream, Head&& head, Tail&&... tail)
  {
    stream << std::forward<Head>(head);
    inline_str_add_args(stream, std::forward<Tail>(tail)...);
  }
} // namespace detail

/** Convert any number of arguments to string and catenate
 *
 * It does pretty much what is advertised. Look at the code if you want to learn
 * some pretty neat modern C++.
 * */
template<class... Args>
std::string inline_str(Args&&... args)
{
  std::stringstream s;
  detail::inline_str_add_args(s, std::forward<Args>(args)...);
  return s.str();
}

/** Convert iterable to string of the form {a,b,c,...} */
template<class T>
std::string vec2str(const T& iterable)
{
  std::stringstream s;
  s << '{';
  for(auto it = begin(iterable);;)
  {
    s << *it;
    if(++it==end(iterable)) break;
    s << ',';
  }
  s << '}';

  return s.str();
}

/** Split string */
inline std::vector<std::string> split(const std::string& s, char sep=' ')
{
  std::vector<std::string> words;
  for(size_t p=0, q=0; p!=s.npos; p=q)
    words.push_back(s.substr(p+(p!=0), (q=s.find(sep, p+1))-p-(p!=0)));
  return words;
}

// =============================================================================
// shared_const_resource

/** Share const objects
 *
 * This small utility allows to share constant objects using smart pointers. The
 * objects are referenced by their constructor parameters (using aggregate init)
 * and can not be modified after their construction. The function get() either
 * returns a smart pointer on the resource if it has already been allocated, or
 * constructs it using aggregate initialization.
 * */
template<typename T, typename ...Params>
struct shared_const_resource
{
  /** The shared_data associates a single .... */
  static std::map<std::tuple<Params...>, std::weak_ptr<const T>> shared;

  /** Return shared_ptr on resource
   *
   * If the resource corresponding to the parameters is not available, it is
   * readily constructed from the parameters using aggregate initialization.
   * */
  std::shared_ptr<const T> get(Params... params)
  {
    if(shared.count({params...})==0)
    {
      std::shared_ptr<const T> ptr(new const T {params...});
      shared[{params...}] = ptr;
      return ptr;
    }
    else
      return shared[{params...}].lock();
  }
};

// declaration of the static map
template<typename T, typename ...Params>
std::map<std::tuple<Params...>, std::weak_ptr<const T>>
shared_const_resource<T, Params...>::shared;

#endif//TOOLS_HPP_
