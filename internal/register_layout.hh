#pragma once

struct QubitRange {
  int start;
  int end;

  int size() const {
    return (end >= start) ? (end - start + 1) : 0;
  }
};

struct RegisterLayout {
  QubitRange target;
  QubitRange control;
};

inline RegisterLayout make_shor_layout(int n_t, int n_c)
{
  // Target register is lower bits; control register follows.
  RegisterLayout layout;
  layout.target.start = 0;
  layout.target.end = n_t - 1;
  layout.control.start = n_t;
  layout.control.end = n_t + n_c - 1;
  return layout;
}
