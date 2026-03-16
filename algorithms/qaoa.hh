#pragma once

#include "algorithms/vqa_qaoa.hh"

// QaoaOptions and QaoaResult are aliases for the canonical VqaQaoa types.
using QaoaOptions = VqaQaoaOptions;
using QaoaResult  = VqaQaoaResult;

QaoaResult run_qaoa_qubo(int n,
                         const std::vector<double>& q,
                         const QaoaOptions& options);
