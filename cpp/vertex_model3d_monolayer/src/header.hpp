#ifndef HEADER_HPP_
#define HEADER_HPP_

// defines project-wide header list to be precompiled
// this allows us to reduce the compile time considerably

#include <iostream>
#include <fstream>
#include <iomanip>
#include <filesystem>
#include <string>
#include <exception>
#include <sstream>
#include <utility>
#include <limits>
#include <random>
#include <map>
#include <vector>
#include <complex>
#include <array>
#include <functional>
#include <memory>
#include <type_traits>
#include "../accessories/include/json.hpp"
#include <boost/program_options.hpp>

// =============================================================================
// Constants

/** display width (change it directly here) */
constexpr unsigned width = 70;
/** An infmaous constant */
// constexpr double Pi = 3.14159265358979323846;

#endif//HEADER_HPP_
