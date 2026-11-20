/*
 * Optional hook that lets a front end (the shell's EXPORT recorder) watch
 * intermediate states inside library algorithms without the algorithms
 * knowing about any front end. Zero cost when no observer is installed.
 */
#pragma once

#include "state.hh"
#include <functional>
#include <string>

namespace qsim {

/** kind is a short machine tag ("init", "h_all", "oracle", "diffusion"); detail is human text. */
using StepObserver = std::function<void(const std::string& kind, const std::string& detail, const State&)>;

inline StepObserver& step_observer_storage()
{
  static StepObserver observer;
  return observer;
}

inline void set_step_observer(StepObserver observer) { step_observer_storage() = std::move(observer); }
inline void clear_step_observer() { step_observer_storage() = nullptr; }

inline void notify_step(const std::string& kind, const std::string& detail, const State& state)
{
  if (step_observer_storage()) step_observer_storage()(kind, detail, state);
}

} // namespace qsim
