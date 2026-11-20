/*
 * Optional hook that lets a front end (the shell's EXPORT recorder) watch
 * intermediate states inside library algorithms without the algorithms
 * knowing about any front end. Zero cost when no observer is installed.
 */
#pragma once

#include "state.hh"
#include <functional>
#include <iomanip>
#include <sstream>
#include <string>

namespace qsim {

struct StepInfo {
  /** "init" (fresh |0..0>), "gate" (name is a shell command such as "CP 1 0 1.57"), or "block". */
  std::string kind;
  /** For "gate": the equivalent shell command. For "block": a descriptive name. */
  std::string name;
  /** Block only: label short enough for a narrow box in the viewer. */
  std::string short_label;
  /** Block only: extra human text (e.g. "iteration 1/2"). */
  std::string detail;
  /** Block only: qubit span; -1 means the whole register. */
  int first_qubit = -1;
  int last_qubit = -1;
};

using StepObserver = std::function<void(const StepInfo&, const State&)>;

inline StepObserver& step_observer_storage()
{
  static StepObserver observer;
  return observer;
}

inline void set_step_observer(StepObserver observer) { step_observer_storage() = std::move(observer); }
inline void clear_step_observer() { step_observer_storage() = nullptr; }

inline void notify_step(const StepInfo& info, const State& state)
{
  if (step_observer_storage()) step_observer_storage()(info, state);
}

inline void notify_init(const State& state)
{
  if (!step_observer_storage()) return;
  StepInfo info;
  info.kind = "init";
  step_observer_storage()(info, state);
}

inline void notify_gate(const std::string& command, const State& state)
{
  if (!step_observer_storage()) return;
  StepInfo info;
  info.kind = "gate";
  info.name = command;
  step_observer_storage()(info, state);
}

/** "NAME q" (e.g. "H 2"). */
inline void notify_gate1(const char* name, int q, const State& state)
{
  if (!step_observer_storage()) return;
  notify_gate(std::string(name) + " " + std::to_string(q), state);
}

/** "NAME a b" (e.g. "CX 0 1"). */
inline void notify_gate2(const char* name, int a, int b, const State& state)
{
  if (!step_observer_storage()) return;
  notify_gate(std::string(name) + " " + std::to_string(a) + " " + std::to_string(b), state);
}

/** Parametric single-qubit gate: "NAME q angle". */
inline void notify_gate_angle(const char* name, int q, double angle, const State& state)
{
  if (!step_observer_storage()) return;
  std::ostringstream cmd;
  cmd << std::setprecision(9) << name << " " << q << " " << angle;
  notify_gate(cmd.str(), state);
}

/** "MEASURE q c". */
inline void notify_measure(int q, int c, const State& state) { notify_gate2("MEASURE", q, c, state); }

/**
 * Mutes the observer for the lifetime of the object. Used around the inner loops of
 * variational algorithms so only the final, optimized circuit is reported.
 */
class StepObserverPause {
 public:
  StepObserverPause() : saved_(step_observer_storage()) { step_observer_storage() = nullptr; }
  ~StepObserverPause() { step_observer_storage() = saved_; }
  StepObserverPause(const StepObserverPause&) = delete;
  StepObserverPause& operator=(const StepObserverPause&) = delete;
 private:
  StepObserver saved_;
};

inline void notify_block(const std::string& name, const std::string& short_label, const std::string& detail,
                         const State& state, int first_qubit = -1, int last_qubit = -1)
{
  if (!step_observer_storage()) return;
  StepInfo info;
  info.kind = "block";
  info.name = name;
  info.short_label = short_label;
  info.detail = detail;
  info.first_qubit = first_qubit;
  info.last_qubit = last_qubit;
  step_observer_storage()(info, state);
}

} // namespace qsim
