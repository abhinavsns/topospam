#ifndef ERROR_MSG_HPP_
#define ERROR_MSG_HPP_

/** Error messages made easy
  *
  * This class makes it easy to construct simple error messages. Usage:
  *
  *   error_msg(A, B, C, ...)
  *
  * automatically converts arguments A, B, C, ... to string and catenate. The resulting message can be obtained with what().
  */
template<typename Tag>
class inline_msg : std::exception
{
  /** The stream to store the msg */
  std::string msg;
  /** Convert to strig and catenate arguments */
  template<class Head>
  void add_args(std::ostream& stream, Head&& head)
  {
    stream << std::forward<Head>(head);
  }
  /** Convert to strig and catenate arguments */
  template<class Head, class... Tail>
  void add_args(std::ostream& stream, Head&& head, Tail&&... tail)
  {
    stream << std::forward<Head>(head);
    add_args(stream, std::forward<Tail>(tail)...);
  }

public:
  /** Convert to strig and catenate arguments */
  template<class... Args>
  inline_msg(Args&&... args)
  {
    std::stringstream s;
    add_args(s, std::forward<Args>(args)...);
    msg = s.str();
  }

  inline_msg(inline_msg& e) = default;
  inline_msg(inline_msg&& e) = default;

  /** Returns msg */
  const char* what() const throw ()
  { return msg.c_str(); }
};

/** Dummy exception (doing nothing) */
class dummy_exception
{};

/** Error message */
using error_msg = inline_msg<struct Error>;
/** Warning message */
using warning_msg = inline_msg<struct Warning>;

#endif//ERROR_MSG_HPP_
