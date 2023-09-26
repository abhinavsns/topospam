#ifndef OPTIONS_HPP_
#define OPTIONS_HPP_

namespace opt = boost::program_options;
//using option_list = std::vector<opt::options_description>;

/** Declare and parse all program options
  *
  * This function takes care of the option parsing. This tend to be a little
  * intricate due to the different ways the user can interact with both the
  * command line and the runcard file. Moreover, we have the complication that
  * the different models generally defines options that can only be parsed after
  * the current model is known.
  * */

std::vector<opt::options_description> GetOptions();


void ParseProgramOptions(int ac, char **av);

/** Process program options
 *
 * This functions makes small adjustments and checks to the inputed program
 * options obtained from ParseProgramOptions and is called just afterwards. Note
 * that most options do not need any processing.
 * */
void ProcessProgramOptions();

/** Print simulation parameters
  *
  * Simply prints all set parameters with a (somewhat) nice formatting. */
void PrintParameters();

void write_parameters_to_jsonfile();

#endif//OPTIONS_HPP_
