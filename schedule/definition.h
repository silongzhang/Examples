#pragma once

#include"header.h"

constexpr int NMAX = 208;							// NMAX > the number of nodes.

constexpr int InvalidIndex = -1;

constexpr auto CENTI = 1e-2;
constexpr auto MILLI = 1e-3;
constexpr auto PPM = 1e-6;
constexpr auto PPB = 1e-9;
constexpr auto InfinityPos = 1e8;
constexpr auto InfinityNeg = -InfinityPos;

constexpr auto CENTI_34 = 0.34;
constexpr auto CENTI_51 = 0.51;
constexpr auto CENTI_67 = 0.67;

enum class SolutionStatus { Unkown, Infeasible, Feasible, Optimal };
