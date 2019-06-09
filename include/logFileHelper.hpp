#ifndef LOG_FILE_HELPER_HPP
#define LOG_FILE_HELPER_HPP
#include "defs.hpp"
#include "utils.hpp"
#include <stdio.h>
#include <tuple>

// Helper constants
char const * const simNames[SIM_COUNT] = { "mibin", "tebin", "gtebin", "params" };
char const * const extensions[2] = { ".bin", ".log" };

// Open files for writing
void openFiles(bool runSims[SIM_COUNT], FILE* logFiles[SIM_COUNT], FILE* binFiles[SIM_COUNT], char const*const outputFileBase);
// Close files
void endFiles(bool runSims[SIM_COUNT], FILE* logFiles[SIM_COUNT], FILE* binFiles[SIM_COUNT], char const*const outputFileBase);

// Write header information in files
void writeFileHeaders(FILE* logFiles[SIM_COUNT], FILE* binFiles[SIM_COUNT],
					  GeneratorParams const& gp, SimParams const& sp,
					  char const*const outputFileBase);

// Read header from binary file
std::tuple<GeneratorParams, SimParams> readBinaryHeader(FILE* file);

// Read flock state contents from binary file
void readBinaryContents(FILE* file, FlockState & fs);

// Write header to binary file
void writeBinaryHeader(FILE* file, GeneratorParams const& gp, SimParams const& sp);

// Write flock state contents to binary file
void writeBinaryContents(FILE* file, FlockState const& fs);

// Calculate order param and write to binary file
void calc_and_write_order_param(size_t N, double const * const h, double rot, bool sims[SIM_COUNT], FILE* binFiles[SIM_COUNT]);

// Overwrite initial state file with new state for next iteration
void reseed(char const*const initialState, GeneratorParams const& gp, SimParams const& sp, FlockState const& fs);
#endif
